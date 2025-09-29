#' Visualize and Compare Cell Population Graphs Across Time Points
#'
#' This script generates graph-based visualizations of cell populations for a single panel 
#' (e.g., Co-stimulation) across multiple time points. Graphs are constructed using 
#' precomputed Sinkhorn distance matrices for each sample, which quantify differences 
#' in cell-type distributions. 
#'
#' The script visualizes each time point individually, overlays graphs from different 
#' time points to highlight changes in population sizes and connectivity, and computes 
#' the Graph Edit Distance (GED) as a quantitative measure of these changes. Outputs 
#' include filtered graph visualizations (SVG files), overlay comparisons, and GED values 
#' for all analyzed samples.
#'
#' Inputs:
#' - Phenotype layout and edge filter matrices
#' - Precomputed Sinkhorn distance matrices for each sample
#' - Final cell-assigned data for each time point
#'
#' Outputs:
#' - Graph visualizations of individual time points
#' - Overlays comparing graphs across time points
#' - GED values summarizing structural differences between graphs

setwd(dirname(this.path::this.path()))
source("../MPMUtilities.R")
setwd(dirname(this.path::this.path()))


pheno_file = "../../Data/Raw Data/Co-stimulation/Normalized Data/Co-stimulation_cell_type.xlsx"
edge_filter = read.xlsx(pheno_file, sheet = "Test_edge_filter", rowNames = FALSE)[, -1]

svglite(paste("../../Output/", "pheno4.svg", sep=""), width = 3.5, height = 3.5)
pheno = phenotype_distance(pheno_file, graph_kind=4, edge_filter_matrix=edge_filter)
dev.off()

pheno_graph_obj = pheno[[1]] # the graph object for the fully-connected phenotype graph
filtered_graph_obj = pheno[[2]] # the graph object for the filtered phenotype graph
vertices_pheno =pheno[[3]] # the data structures for the nodes and edges
edges_pheno = pheno[[4]]
pheno_layout = pheno[[5]] # the layout of that filtered phenotype graph
pheno_filter_matrix = pheno[[6]] # the 0-1 adjacency matrix of the filtered phenotype graph

ged = data.frame(matrix(nrow=0, ncol=3))
colnames(ged) = c("Baseline_After1Vac",	"After1Vac_After3Vac",	"Baseline_After3Vac")

#### folder_path is the path to the panel's final cell type assigned data; this is for Co-stimulation panel
#### folder_name is the time point specific folder name; it can be: TP1, TP2, TP3
#### file_list has all the .xlsx files in the folder (both group 1 and group 2)
folder_path = "../../Data/Raw Data/Co-stimulation/Final Cell Assigned Data"

folder_name_1 = "TP1"
file_list_1 = list.files(path = paste(folder_path, "/", folder_name_1, sep = ""), pattern = ".xlsx", full.names = TRUE)
folder_name_2 = "TP2"
file_list_2 = list.files(path = paste(folder_path, "/", folder_name_2, sep = ""), pattern = ".xlsx", full.names = TRUE)
folder_name_3 = "TP3"
file_list_3 = list.files(path = paste(folder_path, "/", folder_name_3, sep = ""), pattern = ".xlsx", full.names = TRUE)

file_path_1 = file_list_1[1]
file_path_2 = file_list_2[1]
file_path_3 = file_list_3[1]

file_name_1 = strsplit(file_path_1, "/")[[1]][length(strsplit(file_path_1, "/")[[1]])]
file_name_2 = strsplit(file_path_2, "/")[[1]][length(strsplit(file_path_2, "/")[[1]])]
file_name_3 = strsplit(file_path_3, "/")[[1]][length(strsplit(file_path_3, "/")[[1]])]
cat(file_name_1, "\n", file_name_2, "\n ", file_name_3)

sample_1 = read.xlsx(file_path_1, sheet="Sheet 1", rowNames = FALSE)
sample_2 = read.xlsx(file_path_2, sheet="Sheet 1", rowNames = FALSE)
sample_3 = read.xlsx(file_path_3, sheet="Sheet 1", rowNames = FALSE)

#### filter for cells and marker channels we want
sample_1 = subset(sample_1, !(cell_type %in% c("Other cell")))
sample_1 = sample_1[, c(8:ncol(sample_1))]
sample_2 = subset(sample_2, !(cell_type %in% c("Other cell")))
sample_2 = sample_2[, c(8:ncol(sample_2))]
sample_3 = subset(sample_3, !(cell_type %in% c("Other cell")))
sample_3 = sample_3[, c(8:ncol(sample_3))]

#### extracting the saved sinkhorn distance matrix
sink_file_1 = paste0(folder_path, "/", folder_name_1, "/Sink_", strsplit(file_name_1, "\\.")[[1]][1], ".xlsx")
sinkhorn_distance_matrix_1 = as.matrix(read.xlsx(sink_file_1, sheet = "Sheet 1", rowNames=FALSE))

sink_file_2 = paste0(folder_path, "/", folder_name_2, "/Sink_", strsplit(file_name_2, "\\.")[[1]][1], ".xlsx")
sinkhorn_distance_matrix_2 = as.matrix(read.xlsx(sink_file_2, sheet = "Sheet 1", rowNames=FALSE))

sink_file_3 = paste0(folder_path, "/", folder_name_3, "/Sink_", strsplit(file_name_3, "\\.")[[1]][1], ".xlsx")
sinkhorn_distance_matrix_3 = as.matrix(read.xlsx(sink_file_3, sheet = "Sheet 1", rowNames=FALSE))


pheno_layout_1 = add_size_column(sample_1, pheno_layout)
graph_1 = visualizeFilteredGraph(sinkhorn_distance_matrix_1, 
                       sample_1, 
                       graph_lo=pheno_layout_1, 
                       pheno_filter_matrix, 
                       file=file_name_1)
graph_obj_1 = graph_1[[1]]
vertices_1 = graph_1[[2]] 
edges_1 = graph_1[[3]]

pheno_layout_2 = add_size_column(sample_2, pheno_layout)
graph_2 = visualizeFilteredGraph(sinkhorn_distance_matrix_2, 
                                 sample_2, 
                                 graph_lo=pheno_layout_2, 
                                 pheno_filter_matrix, 
                                 file=file_name_2)
graph_obj_2 = graph_2[[1]]
vertices_2 = graph_2[[2]] 
edges_2 = graph_2[[3]]

pheno_layout_3 = add_size_column(sample_3, pheno_layout)
graph_3 = visualizeFilteredGraph(sinkhorn_distance_matrix_3, 
                                 sample_3, 
                                 graph_lo=pheno_layout_3, 
                                 pheno_filter_matrix, 
                                 file=file_name_3)
graph_obj_3 = graph_3[[1]]
vertices_3 = graph_3[[2]] 
edges_3 = graph_3[[3]]

w_range = range(c(edges_1$weight, edges_2$weight, edges_3$weight))

svglite(paste("../../Output/", "Graph_", strsplit(file_name_1, "\\.")[[1]][1], ".svg", sep=""), width = 5, height = 5)
visualizeFilteredGraph_2(vertices_1, edges_1, w_range)
dev.off()
svglite(paste("../../Output/", "Graph_", strsplit(file_name_2, "\\.")[[1]][1], ".svg", sep=""), width = 5, height = 5)
visualizeFilteredGraph_2(vertices_2, edges_2, w_range)
dev.off()
svglite(paste("../../Output/", "Graph_", strsplit(file_name_3, "\\.")[[1]][1], ".svg", sep=""), width = 5, height = 5)
visualizeFilteredGraph_2(vertices_3, edges_3, w_range)
dev.off()

svglite(paste("../../Output/", "GraphCompare_1_2.svg", sep=""), width = 5.5, height = 5.5)
overlayFilteredGraph(vertices_1, edges_1, vertices_2, edges_2,
                     graph_lo=pheno_layout,
                     file_name_1, file_name_2,
                     w_range)
dev.off()

svglite(paste("../../Output/", "GraphCompare_2_3.svg", sep=""), width = 5.5, height = 5.5)
overlayFilteredGraph(vertices_2, edges_2, vertices_3, edges_3,
                     graph_lo=pheno_layout,
                     file_name_2, file_name_3,
                     w_range)
dev.off()

svglite(paste("../../Output/", "GraphCompare_1_3.svg", sep=""), width = 5.5, height = 5.5)
overlayFilteredGraph(vertices_1, edges_1, vertices_3, edges_3,
                     graph_lo=pheno_layout,
                     file_name_1, file_name_3,
                     w_range)
dev.off()


el1 = data.frame(V1=as_edgelist(graph_obj_1)[,1], 
                 V2=as_edgelist(graph_obj_1)[, 2], 
                 weight=E(graph_obj_1)$weight)

el2 = data.frame(V1=as_edgelist(graph_obj_2)[,1], 
                 V2=as_edgelist(graph_obj_2)[, 2], 
                 weight=E(graph_obj_2)$weight)

el3 = data.frame(V1=as_edgelist(graph_obj_3)[,1], 
                 V2=as_edgelist(graph_obj_3)[, 2], 
                 weight=E(graph_obj_3)$weight)

labels1 = data.frame(id=vertices_1$id, label=vertices_1$label)
labels2 = data.frame(id=vertices_2$id, label=vertices_2$label)
labels3 = data.frame(id=vertices_3$id, label=vertices_3$label)

vsizes1 = data.frame(id=vertices_1$id, size=vertices_1$size)
vsizes2 = data.frame(id=vertices_2$id, size=vertices_2$size)
vsizes3 = data.frame(id=vertices_3$id, size=vertices_3$size)

ged12 = sum(abs(vsizes1$size - vsizes2$size)) + sum(abs(el1$weight - el2$weight))
ged23 = sum(abs(vsizes2$size - vsizes3$size)) + sum(abs(el2$weight - el3$weight))
ged13 = sum(abs(vsizes1$size - vsizes3$size)) + sum(abs(el1$weight - el3$weight))
ged = rbind(ged, c(ged12, ged23, ged13))


