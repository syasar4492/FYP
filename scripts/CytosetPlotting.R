library(flowWorkspace)
library(flowCore)
library(flowViz)
library(openCyto)

rm(list = ls())
read_files <- function(path, pattern){
  cs <- load_cytoset_from_fcs(path = path, pattern = pattern, alter.names = TRUE, truncate_max_range = TRUE)
  tmp <- tempfile()
  save_cytoset(cs, tmp)
  cs <- load_cytoset(tmp)
  sn <- sampleNames(cs)
  for (s in sn){
    cf <- cs[[s]]
    cat("Cytoframe", s, "has", dim(cf)[1], "events", "\n")
  }
  return(tmp)
}

apply_transform <- function(obj, channels = c("FSC-H", "FL1-H")) {
  trans_list <- estimateLogicle(obj[[1]], channels)
  return(transform(obj, trans_list))
}


plot_cytoset <- function(cs, xchannel = "FL1-H", ychannel = "FSC-H") {
  # compute ranges from first cytoframe
  xlim <- range(exprs(cs[[1]])[, xchannel])
  ylim <- range(exprs(cs[[1]])[, ychannel])
  
  # convert to flowSet for flowViz::xyplot
  fs <- cytoset_to_flowSet(cs)
  
  # build formula dynamically
  form <- as.formula(paste0("`", ychannel, "` ~ `", xchannel, "`"))
  
  xyplot(form, fs,
         smooth = FALSE,
         nbin = 500,
         xlim = xlim,
         ylim = ylim,
         xlab = paste(xchannel, "(logicle)"),
         ylab = paste(ychannel, "(logicle)"),
         scales = list(
           x = list(alternating = 1),  # only bottom
           y = list(alternating = 1)   # only left
         ))
}


#histogram
plot_1d_density <- function(cs, channel = 'FL1.H', gs = NULL, gate_alias = "nonDebris", per_gate = NULL, global_gate = NULL) {
  fs <- cytoset_to_flowSet(cs)
  f <- as.formula(paste("~", channel))
  p <- densityplot(f, fs, overlap = 0)
  
  if (!is.null(gs)) {
    gate_obj <- gs_pop_get_gate(gs, gate_alias)[[1]]
    # if (!is.null(per_gate)) {
    #   # Store per_gate in local variable to avoid scoping issues
    #   per_gate_vals <- as.numeric(per_gate)
    #   p <- p + latticeExtra::layer(
    #     panel.abline(v = per_gate_vals, col = "red", lty = 2)
    #   )
    # }
    
    if (!is.null(global_gate)) {
      global_gate_val <- as.numeric(global_gate)  
      p <- p + latticeExtra::layer(
        panel.abline(v = global_gate_val, col = "blue", lty = 3)
      )
    }
  }
  
  return(p)
}


#gating dg beads
gating_dg <- function(cs_dg, yChannel = "FL1-H") {
  # Create GatingSet
  gs_dg <- GatingSet(cs_dg)
  fs_all <- gs_pop_get_data(gs_dg)
  sample_ids <- sampleNames(fs_all)
  
  # Store per-sample cutoffs
  gate_cutoffs <- data.frame(
    sample = character(),
    cutoff = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sid in sample_ids) {
    fr <- fs_all[[sid]]
    
    # Compute mindensity gate for this sample
    gate <- mindensity(fr, channel = yChannel)
    
    # Add gate to the GatingSet
    gs_pop_add(gs_dg[[sid]], gate, parent = "root", name = "nonDebris")
    
    # Extract the cutoff value
    cutoff <- gate@min  # mindensity rectangleGate defines a min boundary along the channel
    
    # Store in data frame
    gate_cutoffs <- rbind(gate_cutoffs, data.frame(sample = sid, cutoff = cutoff))
    
    # Print info
    message("Sample ", sid, " cutoff: ", signif(cutoff, 5))
  }
  
  # Recompute the GatingSet
  recompute(gs_dg)
  
  # Global cutoff (minimum of all sample cutoffs)
  global_cutoff <- min(gate_cutoffs$cutoff)
  message("\nGlobal gate cutoff: ", signif(global_cutoff, 5))
  
  # Store global value as attribute
  attr(gs_dg, "global_cutoff") <- global_cutoff
  
  return(list(gs = gs_dg, gate_cutoffs = gate_cutoffs, global_cutoff = global_cutoff))
}


gating_mmb <- function(cs_mmb, yChannel = "FSC-H", smooth_ratio = 1, drop_frac = 0.05) {
  gs_mmb <- GatingSet(cs_mmb)
  fs_all <- gs_pop_get_data(gs_mmb)
  sample_ids <- sampleNames(fs_all)
  
  gate_ranges <- data.frame(
    sample = character(),
    cutoff = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (sid in sample_ids) {
    fr <- fs_all[[sid]]
    vals <- as.vector(exprs(fr)[, yChannel])
    dens <- density(vals, adjust = smooth_ratio)
    
    peak_idx <- which.max(dens$y)
    
    idx_left <- peak_idx
    while (idx_left > 1 && dens$y[idx_left] > drop_frac * max(dens$y)) idx_left <- idx_left - 1
    idx_right <- peak_idx
    while (idx_right < length(dens$y) && dens$y[idx_right] > drop_frac * max(dens$y)) idx_right <- idx_right + 1
    
    min_val <- dens$x[idx_left]
    max_val <- dens$x[idx_right]
    
    gate <- rectangleGate(
      filterId = "nonDebris",
      .gate = structure(list(c(min_val, max_val)), names = yChannel))
    
    
    gs_pop_add(gs_mmb[[sid]], gate, parent = "root", name = "nonDebris")
    
    gate_ranges <- rbind(gate_ranges, data.frame(sample = sid, min_val = min_val, max_val = max_val))
    message("Sample ", sid, " gate: ", signif(min_val, 5), " - ", signif(max_val, 5))
  }
  
  recompute(gs_mmb)
  
  global_cutoff <- min(gate_ranges$min_val)
  global_max <- max(gate_ranges$max_val)
  message("\nGlobal gate range: ", signif(global_cutoff, 5), " - ", signif(global_max, 5))
  
  # Attach as attributes
  attr(gs_mmb, "global_gate_range") <- c(global_cutoff, global_max)
  attr(gs_mmb, "gate_ranges") <- gate_ranges
  
  # Return both gate info and GatingSet
  return(list(gs = gs_mmb, gate_ranges = gate_ranges, global_range = c(global_cutoff, global_max)))
}


# Ungated visualization function
visualize_controls <- function(pattern, path = "data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA") {
  tmp <- read_files(path = path, pattern = pattern)
  cs <- load_cytoset(tmp)
  cs_t <- realize_view(cs)
  cs_t <- apply_transform(cs_t)
  
  # Determine channel based on pattern
  if (pattern == "DG") {
    channel <- "FL1.H"
  } else {
    channel <- "FSC.H"
  }
  
  # Generate plots
  p_2d <- plot_cytoset(cs_t)
  p_1d <- plot_1d_density(cs_t, channel = channel)
  
  return(list(
    pattern = pattern,
    cytoset = cs_t,
    plot_2d = p_2d,
    plot_1d = p_1d
  ))
}

# Usage example:
# Visualize DG controls
dg_vis <- visualize_controls("DG")
print(dg_vis$plot_2d)
print(dg_vis$plot_1d)

mmb_vis <- visualize_controls("MMB")
print(mmb_vis$plot_2d)
print(mmb_vis$plot_1d)


#------------
# Gating Main 
#------------

# Because we're working with cytosets instead of flowsets, workflow 
# changes from flowset -> GatingSet -> transformation to
# cytoset -> transformation -> GatinSet
# which can make flowJo integration harder but produces simpler code for us

rm(list = setdiff(ls(), lsf.str()))
# Define a function for loading, gating, and plotting
process_controls <- function(pattern, path = "data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA") {
  tmp <- read_files(path = path, pattern = pattern)
  cs <- load_cytoset(tmp)
  
  cs_gating <- realize_view(cs)
  cs_gating <- apply_transform(cs_gating)
  
  if (pattern == "DG") {
    res <- gating_dg(cs_gating)
    gs <- res$gs
    global_cutoff <- res$global_cutoff
    gate_cutoffs <- res$gate_cutoffs
    p <- plot_1d_density(cs_gating, channel = "FL1.H", gs = gs, gate_alias = "nonDebris", per_gate = gate_cutoffs$cutoff, global_gate = global_cutoff)
    return(list(pattern = pattern, plot = p, gs = gs, gate_ranges = gate_cutoffs, global_range = global_cutoff))
  } else {
    res <- gating_mmb(cs_gating)
    gs <- res$gs
    global_range <- res$global_range
    gate_ranges <- res$gate_ranges
    
    # For plotting: pass gs only, as before
    p <- plot_1d_density(cs_gating, channel = "FSC.H", gs = gs, gate_alias = "nonDebris", per_gate = c(gate_ranges$min_val, gate_ranges$max_val), global_gate = global_range)
    
    # Return all useful info
    return(list(pattern = pattern, plot = p, gs = gs, 
                gate_ranges = gate_ranges, global_range = global_range))
  }
}


patterns <- c("DG", "3MMB", "6MMB", "8MMB")

results <- lapply(patterns, process_controls)
names(results) <- patterns
gate_summary <- lapply(results, function(x) x$global_range)
saveRDS(gate_summary, "gate_summary_by_type.rds")


for (ptype in names(results)) {
  cat("\nShowing plot for", ptype, "...\n")
  print(results[[ptype]]$plot)
  readline("Press [enter] to continue")
}


