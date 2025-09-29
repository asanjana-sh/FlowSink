# Loading required libraries
load_lib <- c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", "ggnetwork", "ggplot2", 
              "poweRlaw", "imager", "viridis", "plotrix", "openxlsx", "tidyr", "spdep", "maptools", "tmap", "OpenImageR", 
              "dismo", "lctools", "officer", "rvg", "truncnorm", "emdist", "ks", "rlist", "readxl", "OneR", "MASS", 
              "RColorBrewer", "this.path", "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "proxy", 
              "sampling", "collections", "umap", "gsignal", "e1071", "caTools", "caret", "DescTools", "PerformanceAnalytics", 
              "randomForest", "tsne", "rdist", "Rcpp", "Matrix", "mvtnorm", "gridExtra", "matrixcalc", "msos", "reticulate", 
              "data.table", "ggalt", "intergraph", "kohonen", "rdflib", "cluster", "MetaCyto")
# Install missing libraries
install_lib <- load_lib[!load_lib %in% installed.packages()]
if(length(install_lib) > 0) { install.packages(install_lib, dependencies = TRUE)}
# Load libraries
sapply(load_lib, require, character=TRUE)
#BiocManager::install("MetaCyto")
#BiocManager::install("Biobase")


library(devtools)
library(flowCore)
library(FlowSOM)
library(CytoNorm)
library(ggpubr)
library(MetaCyto)
library(Biobase)


set2FrameSelf <- function(flowSet){
  expr = flowCore::fsApply(flowSet, function(x){v=flowCore::exprs(x); return(v)})
  
  eventN = flowCore::fsApply(flowSet,function(x){ v=flowCore::exprs(x); n=nrow(v); return(n)})
  
  Label = flowCore::fsApply(flowSet,function(x){v=flowCore::pData(flowCore::parameters(x)); return(v)})
  
  #annotate expr colnames
  channels = Label[[1]]$name
  antibodies = Label[[1]]$desc
  
  antibodies = toupper(antibodies)
  channels = toupper(channels)
  
  antibodies = sapply(1:length(antibodies), function(i){
    if(is.na(antibodies[i])|antibodies[i]=="NA"){return(channels[i])}else{antibodies[i]}
  })
  
  # colnames(expr) = antibodies
  
  #add sample_id
  sample_id = lapply(1:length(eventN),function(i){rep(i,eventN[i])})
  sample_id = unlist(sample_id)
  expr = cbind(expr,"sample_id"=sample_id)
  
  parameters_df <- data.frame(name = colnames(expr), 
                              desc = c(unname(antibodies), "sample_id"), 
                              range = apply(expr, 2, max) - apply(expr, 2, min),
                              minRange = apply(expr, 2, min),
                              maxRange = apply(expr, 2, max),
                              stringsAsFactors = FALSE)
  parameters_annotated <- AnnotatedDataFrame(parameters_df)
  
  fFrame = flowCore::flowFrame(expr, parameters = parameters_annotated)
  rownames(fFrame@parameters@data) = paste0("$P", 1:nrow(fFrame@parameters@data))
  
  return(fFrame)
}


# Function to merge FCS files by sample ID
merge_fcs_files <- function(directory) {
  # List all FCS files in the directory
  files <- list.files(directory, pattern = "\\.fcs$", full.names = TRUE)
  
  # Extract sample IDs from file names
  sample_ids <- unique(sub("export_(\\d+)_.*\\.fcs", "\\1", basename(files)))
  
  # Loop through each sample ID
  for (sample_id in sample_ids) {
    # Get files with the same sample ID
    sample_files <- files[grep(paste0("export_", sample_id, "_"), files)]
    
    # Read and merge the FCS files
    flow_set = read.flowSet(sample_files)
    flow_frame = set2FrameSelf(flow_set)
    
    # Write the merged FCS file
    output_file <- file.path(directory, paste0("export_", sample_id, "_Lymphocytes.fcs"))
    write.FCS(flow_frame, output_file)
  }
  write.csv(list.files(directory, pattern = "\\.fcs$", full.names = F), paste(directory, "/temp.csv", sep = ""))
  
}

# Specify the directory containing the FCS files
directory <- "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Lymphocytes TP1-2-3/Cytokine_tubes_uncompensated"
# Call the function to merge FCS files
merge_fcs_files(directory)

naming_file = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Lymphocytes TP1-2-3/Cytokine_transformed_what_is_what.xls"
# Function to rename files based on an Excel mapping file
rename_files <- function(directory, excel_file) {
  # Read the Excel file
  naming_data <- read_excel(excel_file, sheet = "Sheet2")
  
  # Loop through each row in the naming data
  for (i in 1:nrow(naming_data)) {
    current_name <- naming_data$`Name after exporting`[i]
    new_name <- naming_data$Rename[i]
    
    # Construct full file paths
    current_path <- file.path(directory, current_name)
    new_path <- file.path(directory, new_name)
    
    # Rename the file
    if (file.exists(current_path)) {
      file.rename(current_path, new_path)
    } else {
      warning(paste("File", current_name, "does not exist in the directory."))
    }
  }
}

# Call the function to rename the files
rename_files(directory, naming_file)

merge_fcs_files(directory)
