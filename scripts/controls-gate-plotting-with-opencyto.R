library(openCyto)
library(ggplot2)
library(ggcyto)
library(flowCore)
library(flowClust)

rm(list = ls())

# ------------------------------
# 1. Custom gate: peak density
# ------------------------------
peak_density_gate <- function(fr, pp_res, channels, adjust = 1.5, drop_frac = 0.05, ...) {
  vals <- exprs(fr)[, channels]
  
  dens <- density(vals, adjust = adjust)
  peak_idx <- which.max(dens$y)
  
  left_idx <- peak_idx
  while (left_idx > 1 && dens$y[left_idx] > drop_frac * max(dens$y)) {
    left_idx <- left_idx - 1
  }
  right_idx <- peak_idx
  while (right_idx < length(dens$y) && dens$y[right_idx] > drop_frac * max(dens$y)) {
    right_idx <- right_idx + 1
  }
  
  min_val <- dens$x[left_idx]
  max_val <- dens$x[right_idx]
  
  gate <- rectangleGate(.gate = structure(list(c(min_val, max_val)), names = channels))
  gate@filterId <- "peak_density_gate"
  return(gate)
}

register_plugins(fun = peak_density_gate, methodName = "peak_density_gate")

# ------------------------------
# 2. Gating for DG samples
# ------------------------------
gate_dg <- function(dg_files, prominence_threshold = 0.3) {
  # --- Step 1: decide method for each file ---
  mapping <- data.frame(file = dg_files, method = NA, stringsAsFactors = FALSE)
  
  for (i in seq_along(dg_files)) {
    f <- dg_files[i]
    fr <- read.FCS(f, truncate_max_range = FALSE, alter.names = TRUE)
    vals <- exprs(fr)[, "FL1.H"]
    
    dens <- density(vals)
    peak_idx <- which(diff(sign(diff(dens$y))) == -2) + 1
    peak_heights <- dens$y[peak_idx]
    
    method <- "flowClust"
    if (length(peak_idx) >= 2) {
      ord <- order(dens$x[peak_idx])
      left_peak  <- peak_heights[ord[1]]
      right_peak <- peak_heights[ord[length(ord)]]
      if (left_peak / right_peak > prominence_threshold) {
        method <- "mindensity"
      }
    }
    mapping$method[i] <- method
  }
  print(mapping)  # diagnostic
  
  # --- Step 2: split files ---
  mindensity_files <- mapping$file[mapping$method == "mindensity"]
  flowclust_files  <- mapping$file[mapping$method == "flowClust"]
  
  # GatingSet for mindensity
  gs_mindensity <- NULL
  if (length(mindensity_files) > 0) {
    cs <- load_cytoset_from_fcs(mindensity_files)
    gs_mindensity <- GatingSet(cs)
    lgcl <- estimateLogicle(gs_mindensity[[1]], c("FSC.H", "FL1.H"))
    gs_mindensity <- transform(gs_mindensity, lgcl)
    
    gs_add_gating_method(
      gs_mindensity,
      alias = "nonDebris",
      pop = "+",
      parent = "root",
      dims = "FL1.H",
      gating_method = "mindensity"
    )
    recompute(gs_mindensity)
    
    gs_add_gating_method(
      gs_mindensity,
      alias = "beads",
      pop = "+",
      parent = "nonDebris",
      dims = "FSC.H,FL1.H",
      gating_method = "flowClust.2d",
      gating_args = list(K = 1)
    )
    recompute(gs_mindensity)
  }
  
  # GatingSet for flowClust
  gs_flowclust <- NULL
  if (length(flowclust_files) > 0) {
    cs <- load_cytoset_from_fcs(flowclust_files)
    gs_flowclust <- GatingSet(cs)
    lgcl <- estimateLogicle(gs_flowclust[[1]], c("FSC.H", "FL1.H"))
    gs_flowclust <- transform(gs_flowclust, lgcl)
    
    gs_add_gating_method(
      gs_flowclust,
      alias = "nonDebris",
      pop = "+",
      parent = "root",
      dims = "FL1.H",
      gating_method = "flowClust",
      gating_args = list(K = 2)
    )
    recompute(gs_flowclust)
    
    gs_add_gating_method(
      gs_flowclust,
      alias = "beads",
      pop = "+",
      parent = "nonDebris",
      dims = "FSC.H,FL1.H",
      gating_method = "flowClust.2d",
      gating_args = list(K = 1)
    )
    recompute(gs_flowclust)
  }
  
  return(list(mapping = mapping,
              gs_mindensity = gs_mindensity,
              gs_flowclust = gs_flowclust))
}


# ------------------------------
# 3. Gating for MMB samples
# ------------------------------
gate_mmb <- function(mmb_files) {
  cs_mmb <- load_cytoset_from_fcs(mmb_files)
  gs_mmb <- GatingSet(cs_mmb)
  lgcl <- estimateLogicle(gs_mmb[[1]], c("FSC.H", "FL1.H"))
  gs_mmb <- transform(gs_mmb, lgcl)
  
  # Step 1: keep central bead peak
  gs_add_gating_method(
    gs_mmb,
    alias = "beads_only",
    pop = "+",
    parent = "root",
    dims = "FSC.H",
    gating_method = "peak_density_gate",
    gating_args = list(adjust = 1.5, drop_frac = 0.05)
  )
  recompute(gs_mmb)
  
  # Step 2: cluster beads within that peak
  gs_add_gating_method(
    gs_mmb,
    alias = "beads",
    pop = "+",
    parent = "beads_only",
    dims = "FSC.H,FL1.H",
    gating_method = "flowClust.2d",
    gating_args = list(K = 1)
  )
  recompute(gs_mmb)
  
  return(gs_mmb)
}

# ------------------------------
# 4. Plot gating results
# ------------------------------
plot_gating_results <- function(gs_list, gate_name = "beads") {
  for (nm in names(gs_list)) {
    if (!is.null(gs_list[[nm]])) {
      message("Plotting: ", nm)
      p <- ggcyto(gs_list[[nm]], aes(x = FL1.H, y = FSC.H)) +
        geom_hex(bins = 500) +
        geom_gate(gate_name) +
        facet_wrap(~ name, scales = "free") +
        ggtitle(paste("Gating results -", nm))
      print(p)
    }
  }
}

# ------------------------------
# 5. Plot 1D histograms
# ------------------------------
plot_1d_histogram <- function(gs, sn, channel, title = NULL) {
  fr <- gh_pop_get_data(gs[[sn]], "root", returnType = "flowFrame")
  autoplot(fr, channel) + ggtitle(title %||% paste(channel, "distribution for", sn))
}

# ------------------------------
# 6. Main function
# ------------------------------
process_controls <- function() {
  fcsPath <- "data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA"
  fcsFiles <- list.files(fcsPath, pattern = "_1\\.fcs$", full.names = TRUE)
  
  print(fcsFiles)
  
  dg_files  <- fcsFiles[grepl("^.*/DG", fcsFiles)]
  mmb_files <- setdiff(fcsFiles, dg_files)
  
  res_dg  <- gate_dg(dg_files)
  res_mmb <- gate_mmb(mmb_files)
  
  gs_all <- list(
    "MMB samples"     = res_mmb,
    "DG mindensity"   = res_dg$gs_mindensity,
    "DG flowClust"    = res_dg$gs_flowclust
  )
  
  p <- plot_gating_results(gs_all, gate_name = "beads")
  
  
  return(gs_all)
}

# ------------------------------
# Run
# ------------------------------
res <- process_controls()
print(res)
# Example: 1D histogram
print(plot_1d_histogram(, "8MMB_1.fcs", "FSC.H", "FSC.H for 8MMB_1.fcs"))
