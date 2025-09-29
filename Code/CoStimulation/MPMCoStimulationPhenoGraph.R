#' DISCLAIMER: The main goal of this script is to compute the Sinkhorn distances 
#' between every pair of cell populations in a sample and save them for future usage. 
#' Since these files are already saved, it is recommended not to run this script 
#' unless required.
#'
#' Compute and Visualize Phenotype Distances Using Sinkhorn Metrics
#'
#' This script analyzes cell-type composition and relationships in a Co-stimulation panel 
#' by computing phenotype distances and visualizing filtered graphs. Specifically, it:
#' 
#' 1. Loads a pre-defined edge filter and cell-type phenotype data from Excel files.
#' 2. Calls `phenotype_distance()` to compute:
#'    - Fully-connected and filtered phenotype graphs.
#'    - Node and edge data structures.
#'    - Filtered graph layout and adjacency matrix.
#' 3. Iterates through all patient sample Excel files for a specified time point.
#'    - Filters out unwanted cell types and selects relevant marker columns.
#'    - Computes cell-type composition and percentages.
#'    - Maps cell types to numeric labels corresponding to the phenotype layout.
#' 4. Computes Sinkhorn distance matrices between sample cell distributions and phenotype layout using `computeSinkhornDistances()`.
#' 5. Optionally saves Sinkhorn distance matrices to Excel files for future usage.
#' 6. Adds size information to the phenotype layout and visualizes the filtered phenotype graph using `visualizeFilteredGraph()`.
#'
#' @note The script relies on `MPMUtilities.R` and assumes that cell-type labels in the data 
#'       match those in the phenotype layout. The filtered graph highlights relationships 
#'       between cell types weighted by Sinkhorn distances.

setwd(dirname(this.path::this.path()))
source("../MPMUtilities.R")
setwd(dirname(this.path::this.path()))


pheno_file = "../../Data/Raw Data/Co-stimulation/Normalized Data/Co-stimulation_cell_type.xlsx"
edge_filter = read.xlsx(pheno_file, sheet = "Test_edge_filter", rowNames = FALSE)[, -1]

# svglite(paste("../../Output/", "pheno1.svg", sep=""), width = 3.5, height = 3.5)
# pheno = phenotype_distance(pheno_file, graph_kind=1, edge_filter_matrix=NULL)
# dev.off()
# 
# svglite(paste("../../Output/", "pheno2.svg", sep=""), width = 3.5, height = 3.5)
# pheno = phenotype_distance(pheno_file, graph_kind=2, edge_filter_matrix=NULL)
# dev.off()
# 
# svglite(paste("../../Output/", "pheno3.svg", sep=""), width = 3.5, height = 3.5)
# pheno = phenotype_distance(pheno_file, graph_kind=3, edge_filter_matrix=NULL)
# dev.off()

# svglite(paste("../../Output/", "pheno4.svg", sep=""), width = 3.5, height = 3.5)
pheno = phenotype_distance(pheno_file, graph_kind=4, edge_filter_matrix=edge_filter)
# dev.off()

pheno_graph_obj = pheno[[1]] # the graph object for the fully-connected phenotype graph
filtered_graph_obj = pheno[[2]] # the graph object for the filtered phenotype graph
vertices_pheno =pheno[[3]] # the data structures for the nodes and edges
edges_pheno = pheno[[4]]
pheno_layout = pheno[[5]] # the layout of that filtered phenotype graph
pheno_filter_matrix = pheno[[6]] # the 0-1 adjacency matrix of the filtered phenotype graph


#### folder_path is the path to the panel's final cell type assigned data; this is for Co-stimulation panel
#### folder_name is the time point specific folder name; it can be: TP1, TP2, TP3
#### file_list has all the fcs files in the folder (both group 1 and group 2)
folder_path = "../../Data/Raw Data/Co-stimulation/Final Cell Assigned Data"
folder_name = "TP3"
file_list = list.files(path = paste(folder_path, "/", folder_name, sep = ""), pattern = ".xlsx", full.names = TRUE)

for(file_path in file_list){
  # file_path = file_list[1]
  print(file_path)
  
  file_name = strsplit(file_path, "/")[[1]][length(strsplit(file_path, "/")[[1]])]
  sample = read.xlsx(file_path, sheet="Sheet 1", rowNames = FALSE)
  
  #### filter for cells and marker channels we want
  sample = subset(sample, !(cell_type %in% c("Other cell")))
  sample = sample[, c(8:ncol(sample))]
  
  #### check assigned cell type percentage 
  t1 = data.frame(table(sample$cell_type))
  colnames(t1) = c("Cell Type", "Count")
  t1$Percent = round(t1$Count / sum(t1$Count) * 100, 3)
  
  cell_type_map = data.frame(cell_type=pheno_layout$v_label, numeric_label=c(1:length(pheno_layout$v_label)))
  
  #### compute and visualize sinkhorn space for true cell types
  res = computeSinkhornDistances(sample, cell_type_map)
  sinkhorn_distance_matrix = res[[1]]
  ds = res[[2]]
  
  #### saving the sinkhorn distance matrix
  sink_file = paste0(folder_path, "/", folder_name, "/Sink_", strsplit(file_name, "\\.")[[1]][1], ".xlsx")
  write.xlsx(sinkhorn_distance_matrix, sink_file, rowNames=FALSE)
}

pheno_layout = add_size_column(sample, pheno_layout)
par(mar=c(0, 0, 0, 0))
visualizeFilteredGraph(sinkhorn_distance_matrix, 
                       ds, 
                       graph_lo=pheno_layout, 
                       pheno_filter_matrix, 
                       file=file_name)
  
