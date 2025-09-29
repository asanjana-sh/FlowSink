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


#### This code is only for the Co-stimulation panel, so folder path will always be fixed
#### Folder name will change for the 3 time point data: TP1, TP2, or TP3
folder_path = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation"
folder_name = "TP3"

#### Day 1 data file names start with 'Specimen' and Day 2 data file names start with 'co-stimb'
#### As they have separate compensation matrix, we process them separately
#### note that, the raw data folder also provides separate .wsp (workspace) files that we don't use in R
file_list_1 = list.files(path = paste(folder_path, "/", folder_name, sep = ""), pattern = "Specimen.*.fcs", full.names = TRUE)
file_list_2 = list.files(path = paste(folder_path, "/", folder_name, sep = ""), pattern = "co-stim.*.fcs", full.names = TRUE)

#### Day 1, file_list_1, 'Co-stimulation.csv' compensation matrix file
for(file in file_list_1){
  print(file)
  # General pipeline for preprocessing and quality control with PeacoQC
  # Read in raw fcs file
  ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)
  
  # Define channels where the margin events should be removed and on which the quality control should be done
  marker_channels = colnames(flowCore::keyword(ff)$SPILL)
  channels = c("FSC-A", "SSC-A", marker_channels)
  ff <- RemoveMargins(ff=ff, channels=channels, output="frame")
  
  spill_matrix = read.csv(paste(folder_path, "/", "Co-stimulation.csv", sep=""))
  spill_matrix = spill_matrix[, 2:ncol(spill_matrix)]
  colnames(spill_matrix) = marker_channels
  
  # Compensate and transform the data
  ff <- flowCore::compensate(ff, spill_matrix) # ff@description$SPILL
  asinhTrans <- flowCore::arcsinhTransform(transformationId="asinhTrans", a=0, b=(1/150), c=0)
  ff <- flowCore::transform(ff, transformList(marker_channels, asinhTrans) )
  
  # FSC-A rescaling by calculate the combined range for all marker channels
  marker_ranges <- sapply(marker_channels, function(channel) range(exprs(ff)[, channel]))
  combined_range <- range(marker_ranges)
  exprs(ff)[, "FSC-A"] <- scales::rescale(exprs(ff)[, "FSC-A"], to = combined_range)
  
  #Run PeacoQC
  PeacoQC_res <- PeacoQC(ff, channels,
                         determine_good_cells="all",
                         plot = TRUE,
                         save_fcs=TRUE,
                         output_directory = folder_path,
                         name_directory= paste("PeacoQC_", folder_name, sep=""))
}

#### Day 2, file_list_2, 'Co-stimulation_2.csv' compensation matrix file
for(file in file_list_2){
  print(file)
  # General pipeline for preprocessing and quality control with PeacoQC
  # Read in raw fcs file
  ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)
  
  # Define channels where the margin events should be removed and on which the quality control should be done
  marker_channels = colnames(flowCore::keyword(ff)$SPILL)
  channels = c("FSC-A", "SSC-A", marker_channels)
  ff <- RemoveMargins(ff=ff, channels=channels, output="frame")
  
  spill_matrix = read.csv(paste(folder_path, "/", "Co-stimulation_2.csv", sep=""))
  spill_matrix = spill_matrix[, 2:ncol(spill_matrix)]
  colnames(spill_matrix) = marker_channels
  
  # Compensate and transform the data
  ff <- flowCore::compensate(ff, spill_matrix) # ff@description$SPILL
  asinhTrans <- flowCore::arcsinhTransform(transformationId="asinhTrans", a=0, b=(1/150), c=0)
  ff <- flowCore::transform(ff, transformList(marker_channels, asinhTrans) )
  
  # FSC-A rescaling by calculate the combined range for all marker channels
  marker_ranges <- sapply(marker_channels, function(channel) range(exprs(ff)[, channel]))
  combined_range <- range(marker_ranges)
  exprs(ff)[, "FSC-A"] <- scales::rescale(exprs(ff)[, "FSC-A"], to = combined_range)
  
  #Run PeacoQC
  PeacoQC_res <- PeacoQC(ff, channels,
                         determine_good_cells="all",
                         plot = TRUE,
                         save_fcs=TRUE,
                         output_directory = folder_path,
                         name_directory= paste("PeacoQC_", folder_name, sep=""))
}

#### After saving the quality controlled fcs files, we merge all the plots in a single pdf file
png_dir_path = paste(folder_path, "/PeacoQC_", folder_name, "/PeacoQC_plots", sep="")
peacoQCPDF <- function(png_dir){
  # Load necessary libraries
  library(grid)
  library(gridExtra)
  
  # Get a list of all PNG files in the directory
  png_files <- list.files(png_dir, pattern = "\\.png$", full.names = TRUE)
  
  # Create a PDF file to save the output with landscape orientation
  pdf(paste(png_dir, "/output.pdf", sep=""), width = 11, height = 8.5)
  
  # Loop through each PNG file
  for (png_file in png_files) {
    # Extract the file name without the extension
    file_name <- tools::file_path_sans_ext(basename(png_file))
    
    # Read the PNG file
    img <- rasterGrob(png::readPNG(png_file), interpolate = TRUE)
    
    # Create a page with the image and the title
    grid.newpage()
    grid.draw(img)
    grid.text(file_name, x = 0.5, y = 0.95, gp = gpar(fontsize = 20, fontface = "bold"))
  }
  
  # Close the PDF file
  dev.off()
}
peacoQCPDF(png_dir_path)


