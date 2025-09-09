library(flowCore)
library(flowViz)
library(ggplot2)
library(ggcyto)
library(flowClust)
library(mclust)

read_fcs <- function(filename) {
  if (isFCSfile(filename)) {
    ff <- read.FCS(filename, emptyValue = TRUE, dataset = 1, alter.names = TRUE, truncate_max_range = FALSE)
  }
  else {
    stop("File is not of .fcs type")
  }
  return(ff)
  
}

process_ff <- function(ff){
  tf <- arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=1/150, c=0)
  ff <- transform(ff, transformList(c("FSC.H", "FSC.A", "FL1.H", "FL1.A"), tf))
  return(ff)
}

downsampling <- function(ff, size){
  n <- nrow(ff)
  idx <- sample(n, size)
  ff_sub <- ff[idx, ]
  return(ff_sub)
}

clustering <- function(ff, k = 7){
  
  transformed_ff <- process_ff(ff)
  clusters <- flowClust(transformed_ff, varNames = c("FL1.H", "FSC.H"), K = k, lambda = 0.3)
  labels <- clusters@label
  
  df <- as.data.frame(exprs(ff))
  df$cluster <- factor(labels)
  return(list(df = df, labels = labels))
}

plot <- function(ff){
  plot_ff <- ggplot(ff, aes(x = FL1.H, y=FSC.H, color = cluster))+
    geom_point(alpha=0.2, size=0.5)+
    geom_hex(bins = 500)+
    scale_x_logicle()+
    scale_y_logicle()
  return(plot_ff)
}

split_by_cluster <- function(df, labels){
  
  split_df <- split(df, df$cluster)
  return(split_df)
  
}

gmm <- function(cluster_df){
  gmm_model <- Mclust(cluster_df[, c("FL1.H", "FSC.H")], G = 1:7)
  return(gmm_model)
}

main <- function(){
  filename <- "06.11.2025_TE varying NaCl concs_EFG/1M NaCl_Beads.fcs"
  ff <- read_fcs(filename)
  
  
  res <- clustering(ff)
  clustered_df <- res$df
  labels <- res$labels
  
  p <- plot(clustered_df)
  split_df <- split_by_cluster(clustered_df, labels)
  cluster2 <- split_df[["2"]]
  gmm_model <- gmm(cluster2)
  summary(gmm_model)
  
}

main()








