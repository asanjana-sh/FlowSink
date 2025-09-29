# BiocManager::install("ComplexHeatmap")
# devtools::install_github("saeyslab/PeacoQC")

# Loading required libraries
load_lib <- c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", 
              "poweRlaw", "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", 
              "dismo", "lctools", "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", 
              "RColorBrewer", "this.path", "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "proxy", 
              "sampling", "collections", "umap", "gsignal", "e1071", "caTools", "caret", "DescTools", "PerformanceAnalytics", 
              "randomForest", "tsne", "rdist", "Rcpp", "Matrix", "mvtnorm", "gridExtra", "matrixcalc", "msos", "reticulate", 
              "data.table", "ggalt", "intergraph", "kohonen", "rdflib", "cluster", "autoimage", "ggforce", "ggpubr", "stringr",
              "grid", "scales")

# Install missing libraries
install_lib <- load_lib[!load_lib %in% installed.packages()]
if(length(install_lib) > 0) { install.packages(install_lib, dependencies = TRUE)}
# Load libraries
sapply(load_lib, require, character=TRUE)

# Additional libraries
library(flowCore)
library(FlowSOM)
library(flowWorkspace)
library(kohonen)
library(PeacoQC)


folder_path = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition"
folder_name = "TP3"
file_list = list.files(path = paste(folder_path, "/", folder_name, sep = ""), full.names = TRUE)

for(file in file_list){
  print(file)
  # General pipeline for preprocessing and quality control with PeacoQC
  # Read in raw fcs file
  ff <- flowCore::read.FCS(file)
  
  data = as.data.frame(ff@exprs)
  twoDFCSPlot(data, "Time", "FSC-A", FALSE, FALSE, "No Gate")
  
  # Define channels where the margin events should be removed
  # and on which the quality control should be done
  channels <- c(1: ncol(ff@exprs))
  ff <- RemoveMargins(ff=ff, channels=channels, output="frame")
  
  data = as.data.frame(ff@exprs)
  twoDFCSPlot(data, "Time", "FSC-A", FALSE, FALSE, "No Gate")
  
  spill_matrix = read.csv(paste(folder_path, "/", "Co-inhibition.csv", sep=""))
  spill_matrix = spill_matrix[, 2:ncol(spill_matrix)]
  colnames(spill_matrix) = colnames(ff@description$SPILL)
  
  # Compensate and transform the data
  marker_channels = colnames(flowCore::keyword(ff)$SPILL)
  ff <- flowCore::compensate(ff, spill_matrix) # ff@description$SPILL
  
  asinhTrans <- flowCore::arcsinhTransform(transformationId="asinhTrans", a=0, b=(1/150), c=0)
  ff <- flowCore::transform(ff, transformList(marker_channels, asinhTrans) )
  
  # FSC-A rescaling
  # Calculate the combined range for all marker channels
  marker_ranges <- sapply(marker_channels, function(channel) range(exprs(ff)[, channel]))
  combined_range <- range(marker_ranges)
  exprs(ff)[, "FSC-A"] <- scales::rescale(exprs(ff)[, "FSC-A"], to = combined_range)
  
  #Run PeacoQC
  PeacoQC_res <- PeacoQC(ff, channels,
                         determine_good_cells="all",
                         save_fcs=TRUE,
                         output_directory = folder_path,
                         name_directory= paste("PeacoQC_", folder_name, sep=""))
}



