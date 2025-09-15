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

apply_transform <- function(cs, channels = c("FSC-H", "FL1-H")) {
  lgcl <- estimateLogicle(cs[[1]], channels)  # estimate from first sample
  cs_t <- transform(cs, lgcl)                 # apply to whole cytoset
  return(cs_t)
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


plot_1d_histograms <- function(cs, channel){
  sn <- sampleNames(cs)
  for (s in sn){
    cf <- cs[[s]]
    cat("Cytoframe", s, "dimensions (events parameters):", dim(cf), "\n")
    plot(cf, "FL1-H", smooth = TRUE)
  }
}


tmp_dg <- read_files(path = 'data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA', pattern = 'DG')
cs_dg <- load_cytoset(tmp_dg)
cs_dg_t <- realize_view(cs_dg)
cs_dg_t <- apply_transform(cs_dg_t)
plot_cytoset(cs_dg_t)

tmp_mmb <- read_files(path = 'data/raw/09.09.2025_separated controls_0.25M NaCl TE Tween BSA', pattern = 'MMB')
cs_mmb <- load_cytoset(tmp_mmb)
cs_mmb_t <- realize_view(cs_mmb)
cs_mmb_t <- apply_transform(cs_mmb_t)
plot_cytoset(cs_mmb_t)

