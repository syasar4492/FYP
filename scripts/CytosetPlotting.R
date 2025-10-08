library(flowWorkspace)
library(flowCore)
library(flowViz)

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
plot_1d_density <- function(cs, channel = 'FL1.H', gs = NULL, gate_alias = "nonDebris"){
  fs <- cytoset_to_flowSet(cs)
  f <- as.formula(paste("~", channel))
  
  if (!is.null(gs)){
    gate_obj <- gs_pop_get_gate(gs, gate_alias)[[1]]
    if (gate_obj@max == Inf) { 
      threshold <- gate_obj@min  
      cat("Min threshold by mindensity:", threshold, "\n")
      p <- densityplot(f, fs, overlap = 0, refline = threshold)
    }
    else {
      cat("Min:", gate_obj@min, "\n")
      cat("Max:", gate_obj@max, "\n")
      min_val = gate_obj@min
      max_val = gate_obj@max
      p <- densityplot(f, fs, overlap = 0) +
        latticeExtra::layer(
        panel.abline(v = c(min_val, max_val), col = "red", lty = 2),
        data = list(min_val = min_val, max_val = max_val))
    }
  }
   else {
    p <- densityplot(f, fs, overlap = 0)
  }
  p
}

#gating dg beads
gating_dg <- function(cs_dg){
  gs_dg <- GatingSet(cs_dg)
  
  gs_add_gating_method(
    gs_dg,
    alias = "nonDebris",
    pop = "+",
    parent = "root",
    dims = "FL1.H",
    gating_method = "mindensity"
  )
  recompute(gs_dg)
  gate_obj <- gs_pop_get_gate(gs_dg, "nonDebris")[[1]]
  cat("Type of gate is:", class(gate_obj), "\n")
  cat("Slot names are:", slotNames(gate_obj), "\n")
  if (class(gate_obj) == "rectangleGate"){
    min_val <- gate_obj@min
    max_val <- gate_obj@max
    filterId <- gate_obj@filterId
    cat("Min: ", min_val, "\n")
    cat("Max: ", max_val, "\n")
    cat("FilterID:", filterId, "\n")
  }
  return(gs_dg)
}

gating_mmb <- function(cs_mmb){
  
  gs_mmb <- GatingSet(cs_mmb)
  
  
  # for register_plugins to work later:
  #pp_res is preprocessing results, unused here, channels is multiple even though we use one
  peak_density_gate <- function(fr, pp_res, yChannel = "FSC.H", filterId = "peak_density_gate", smooth_ratio = 1, drop_frac = 0.05, ...) { 
    vars = exprs(fr)[, yChannel]
    dens = density(vars, adjust = smooth_ratio)
    
    peak_idx <- which.max(dens$y)
    
    idx_left <- peak_idx
    while (idx_left > 1 && dens$y[idx_left] > drop_frac * max(dens$y)){
      idx_left <- idx_left - 1
    }
    
    idx_right <- peak_idx
    while (idx_right < length(dens$y) && dens$y[idx_right] > drop_frac * max(dens&y)){
      idx_right <- idx_right + 1
    }
    
    min_val = dens$x[idx_left]
    max_val = dens$x[idx_right]
    
    gate <- rectangleGate(
      filterId = "peak_density_gate",
      .gate = structure(list(c(min_val, max_val)), names = yChannel)
    )
    
    return(gate)
  }
  
  register_plugins(fun = peak_density_gate, methodName = "peak_density_gate")
  
  gs_add_gating_method(
    gs_mmb,
    alias = "nonDebris",
    pop = "+",
    parent = "root",
    dims = "FSC.H",
    gating_method = "peak_density_gate",
    gating_args = list(smooth_ratio = 1, drop_frac = 0.05)
  )
  recompute(gs_mmb)
  
  gate_obj <- gs_pop_get_gate(gs_mmb, "nonDebris")[[1]]
  cat("Type of gate is:", class(gate_obj), "\n")
  cat("Slot names are:", slotNames(gate_obj), "\n")
  if (class(gate_obj) == "rectangleGate"){
    min_val <- gate_obj@min
    max_val <- gate_obj@max
    cat("Peak density gate limits:", min_val, "to", max_val, "\n")
    filterId <- gate_obj@filterId
    cat("FilterID:", filterId)
  }
  return(gs_mmb)
}


#ungated main
rm(list = setdiff(ls(), lsf.str()))
tmp_dg <- read_files(path = 'data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA', pattern = 'DG')
cs_dg <- load_cytoset(tmp_dg)
cs_dg_t <- realize_view(cs_dg)
cs_dg_t <- apply_transform(cs_dg_t)
plot_cytoset(cs_dg_t)
plot_1d_density(cs_dg_t)

tmp_mmb <- read_files(path = 'data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA', pattern = 'MMB')
cs_mmb <- load_cytoset(tmp_mmb)
cs_mmb_t <- realize_view(cs_mmb)
cs_mmb_t <- apply_transform(cs_mmb_t)
plot_cytoset(cs_mmb_t)
plot_1d_density(cs_mmb_t, "FSC.H")


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
  # Load files matching the pattern
  tmp <- read_files(path = path, pattern = pattern)
  
  # Convert to cytoset
  cs <- load_cytoset(tmp)
  
  # Apply transformations
  cs_gating <- realize_view(cs)
  cs_gating <- apply_transform(cs_gating)
  
  # Apply gating
  if (pattern == "DG"){
    gs <- gating_dg(cs_gating)
    p <- plot_1d_density(cs_gating, channel = "FL1.H", gs = gs, gate_alias = "nonDebris")
  }
  else {
    gs <- gating_mmb(cs_gating)
    p <- plot_1d_density(cs_gating, channel = "FSC.H", gs = gs, gate_alias = "nonDebris")
  }
  
  # Plot
  return(p)
}

# List of particle types
patterns <- c("DG", "3MMB", "6MMB", "8MMB")

# Loop over particle types and generate plots
plots <- lapply(patterns, process_controls)

# Show each plot one at a time
for (p in plots) {
  print(p)
  # optional pause if you want to step through
  readline("Press [enter] to continue")
}





