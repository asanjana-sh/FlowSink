#' Extract and Compile Channel Information from Co-inhibition/ Co-stimulation/ Cytokine FCS Files
#'
#' This script reads raw FCS files from a panel for a specified time point 
#' (e.g., TP1) and compiles information about each channel, for fluor::marker pair verification. Specifically, it:
#' 
#' 1. Initializes an empty data frame to store channel names and descriptions.
#' 2. Iterates over all FCS files in the designated folder.
#' 3. Reads each FCS file using `flowCore::read.FCS`.
#' 4. Extracts the channel name and description for each parameter and stores them 
#'    as "name::description" strings.
#' 5. Builds a row-wise data frame where each row corresponds to a sample/FCS file, 
#'    and columns correspond to channels.
#' 6. Assigns informative row names based on the time point and MCV identifier, and 
#'    generic column names for the channels.
#'
#' @note The script relies on `MPMUtilities.R` for auxiliary functions (if any) 
#'       and requires the `flowCore` package to read FCS files.

setwd(dirname(this.path::this.path()))
source("MPMUtilities.R")
setwd(dirname(this.path::this.path()))


channel_info = data.frame(matrix(nrow = 0, ncol=21))
row_names = c()

folder_path = "../Data/Raw Data/Cytokine"
folder_name = "TP1"
file_list = list.files(path = paste(folder_path, "/", folder_name, sep = ""), full.names = TRUE)

for(file in file_list){
  cn = paste0(folder_name, "_MCV", strsplit(strsplit(file, "MCV")[[1]][2], "_")[[1]][1])
  print(cn)
  row_names = c(row_names, cn)

    # Read in raw fcs file
  ff <- flowCore::read.FCS(file)
  
  channel_info = rbind(channel_info, 
                       paste0(ff@parameters@data$name, "::", ff@parameters@data$desc) )
}

rownames(channel_info) = row_names
colnames(channel_info) = paste0("col", c(1:ncol(channel_info)))

