#### The file naming convention of the original FCS files was quite confusing. This code renames the FCS files
#### for better understanding what the files contain. Then separates into 3 subfolders for the vaccination timepoints

# Loading required libraries
load_lib <- c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", 
              "poweRlaw", "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", 
              "dismo", "lctools", "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", 
              "RColorBrewer", "this.path", "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "proxy", 
              "sampling", "collections", "umap", "gsignal", "e1071", "caTools", "caret", "DescTools", "PerformanceAnalytics", 
              "randomForest", "tsne", "rdist", "Rcpp", "Matrix", "mvtnorm", "gridExtra", "matrixcalc", "msos", "reticulate", 
              "data.table", "ggalt", "intergraph", "kohonen", "rdflib", "cluster")

# Install missing libraries
install_lib <- load_lib[!load_lib %in% installed.packages()]
if(length(install_lib) > 0) { install.packages(install_lib, dependencies = TRUE)}
# Load libraries
sapply(load_lib, require, character=TRUE)


#### the name of the sheet here needs to be changed for each panel (Sheet1-Sheet6)
metadata = read.xlsx("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/220707_metadata.xlsx",
                     sheet = "Sheet6")[-c(22:24, 31:33), c(2:6)]

# Define the folder path (needs to changed for each panel as well)
folder_path <- "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/DCMyeloid"

# List all .fcs files in the subfolders
fcs_files <- list.files(path = folder_path, pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE)
fcs_files_2 <- list.files(path = folder_path, pattern = "\\.fcs$", recursive = TRUE, full.names = FALSE)


# Loop over each .fcs file
for (f in c(1:length(fcs_files_2))) {
  if(fcs_files_2[f] %in% metadata$Filename){
    print(fcs_files_2[f])
    matching_index <- which(metadata$Filename == fcs_files_2[f])
    print(matching_index)
    print(metadata$New_Filename[matching_index])
    
    new_file_name = paste(strsplit(fcs_files[f], "DCMyeloid/")[[1]][1], "DCMyeloid/", metadata$New_Filename[matching_index], 
                          sep = "")
    file.rename(fcs_files[f], new_file_name)
  }else{
    print("meta data not found!\n")
  }
}

fcs_files <- list.files(path = folder_path, pattern = "\\.fcs$", recursive = TRUE, full.names = TRUE)
# Define the folder paths
tp1_folder <- "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/DCMyeloid/TP1/"
tp2_folder <- "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/DCMyeloid/TP2/"
tp3_folder <- "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/DCMyeloid/TP3/"

# Create the folders if they don't exist
dir.create(tp1_folder, showWarnings = FALSE)
dir.create(tp2_folder, showWarnings = FALSE)
dir.create(tp3_folder, showWarnings = FALSE)

# Loop over each file in fcs_files
for (file in fcs_files) {
  if (grepl("VAC_1", file)) {
    # Move the file to TP1 folder
    file.rename(file, file.path(tp1_folder, basename(file)))
  } else if (grepl("VAC_2", file)) {
    # Move the file to TP2 folder
    file.rename(file, file.path(tp2_folder, basename(file)))
  } else if (grepl("VAC_3", file)) {
    # Move the file to TP3 folder
    file.rename(file, file.path(tp3_folder, basename(file)))
  }
}

