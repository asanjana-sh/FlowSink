#' Panel GED Analysis and MDS Visualization
#'
#' Computes pairwise Graph Edit Distances (GED) between single-cell samples
#' across three time points, performs non-metric MDS for visualization, and
#' conducts adonis2 tests to evaluate temporal and response-based differences.
#'
#' Workflow:
#' - Load phenotype definitions and filtered adjacency matrix.
#' - Read normalized cell type assignments and precomputed Sinkhorn distance matrices.
#' - Filter unwanted cells and ensure matrices match phenotype graph size.
#' - Compute pairwise GED and symmetrize the matrix.
#' - Perform MDS and annotate plots with patient, time, and response data.
#' - Conduct statistical tests on GED distances.
#'
#' Outputs:
#' - GED matrix Excel file.
#' - MDS plots with response/time annotations.
#' - Statistical summaries (adonis2).
#'
#' Note: Precomputed Sinkhorn distance matrices are required; recalculation is optional.


setwd(dirname(this.path::this.path()))
source("../MPMUtilities.R")
setwd(dirname(this.path::this.path()))


pheno_file = "../../Data/Raw Data/Co-stimulation/Normalized Data/Co-stimulation_cell_type.xlsx"
edge_filter = read.xlsx(pheno_file, sheet = "Test_edge_filter", rowNames = FALSE)[, -1]

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

folder_name_1 = "TP1"
file_list_1 = list.files(path = paste(folder_path, "/", folder_name_1, sep = ""), pattern = "^Norm.*\\.xlsx$", full.names = TRUE)
sink_list_1 = list.files(path = paste(folder_path, "/", folder_name_1, sep = ""), pattern = "^Sink.*\\.xlsx$", full.names = TRUE)

folder_name_2 = "TP2"
file_list_2 = list.files(path = paste(folder_path, "/", folder_name_2, sep = ""), pattern = "^Norm.*\\.xlsx$", full.names = TRUE)
sink_list_2 = list.files(path = paste(folder_path, "/", folder_name_2, sep = ""), pattern = "^Sink.*\\.xlsx$", full.names = TRUE)

folder_name_3 = "TP3"
file_list_3 = list.files(path = paste(folder_path, "/", folder_name_3, sep = ""), pattern = "^Norm.*\\.xlsx$", full.names = TRUE)
sink_list_3 = list.files(path = paste(folder_path, "/", folder_name_3, sep = ""), pattern = "^Sink.*\\.xlsx$", full.names = TRUE)

file_list = c(file_list_1, file_list_2, file_list_3)
sink_list = c(sink_list_1, sink_list_2, sink_list_3)
file_name_list = list()
sample_list = list()
sd_mat_list = list()

for(i in c(1:length(file_list))){
  file_path = file_list[i]
  file_name = strsplit(file_path, "/")[[1]][length(strsplit(file_path, "/")[[1]])]
  file_name_list[[i]] = file_name
  cat(file_name, "\n")
  
  sample = read.xlsx(file_path, sheet="Sheet 1", rowNames = FALSE)
  sample = subset(sample, !(cell_type %in% c("Other cell")))
  sample = sample[, c(8:ncol(sample))]
  sample_list[[i]] = sample
  
  sink_file = sink_list[i]
  sinkhorn_distance_matrix = as.matrix(read.xlsx(sink_file, sheet = "Sheet 1", rowNames=FALSE))
  sd_mat_list[[i]] = sinkhorn_distance_matrix
}

# Find indices of matrices that are not 13x13
invalid_indices <- which(sapply(sd_mat_list, function(mat) { dim(mat)[1] != 13 || dim(mat)[2] != 13}))
if(length(invalid_indices) != 0){
  # Remove invalid matrices from the list
  sd_mat_list <- sd_mat_list[-invalid_indices]
  sample_list = sample_list[-invalid_indices]
  file_name_list = file_name_list[-invalid_indices]
}

# sample_count = length(sample_list)
# ged_mat = matrix(0, sample_count, sample_count, byrow = TRUE)
# 
# for(i in c(1:sample_count)){
#   for(j in c(i:sample_count)){
#     if(i != j){
#       sample_1 = sample_list[[i]]
#       sinkhorn_distance_matrix_1 = sd_mat_list[[i]]
#       file_name_1 = file_name_list[[i]]
#       pheno_layout_1 = add_size_column(sample_1, pheno_layout)
#       graph_1 = visualizeFilteredGraph(sinkhorn_distance_matrix_1, sample_1, graph_lo=pheno_layout_1, 
#                                        pheno_filter_matrix, file=file_name_1)
#       graph_obj_1 = graph_1[[1]]
#       vertices_1 = graph_1[[2]] 
#       edges_1 = graph_1[[3]]
#       
#       sample_2 = sample_list[[j]]
#       sinkhorn_distance_matrix_2 = sd_mat_list[[j]]
#       file_name_2 = file_name_list[[j]]
#       pheno_layout_2 = add_size_column(sample_2, pheno_layout)
#       graph_2 = visualizeFilteredGraph(sinkhorn_distance_matrix_2, sample_2, graph_lo=pheno_layout_2, 
#                                        pheno_filter_matrix, file=file_name_2)
#       graph_obj_2 = graph_2[[1]]
#       vertices_2 = graph_2[[2]] 
#       edges_2 = graph_2[[3]]
#       
#       el1 = data.frame(V1=as_edgelist(graph_obj_1)[,1], V2=as_edgelist(graph_obj_1)[, 2], 
#                        weight=E(graph_obj_1)$weight)
#       
#       el2 = data.frame(V1=as_edgelist(graph_obj_2)[,1], V2=as_edgelist(graph_obj_2)[, 2], 
#                        weight=E(graph_obj_2)$weight)
#       
#       
#       labels1 = data.frame(id=vertices_1$id, label=vertices_1$label)
#       labels2 = data.frame(id=vertices_2$id, label=vertices_2$label)
#       
#       vsizes1 = data.frame(id=vertices_1$id, size=vertices_1$size)
#       vsizes2 = data.frame(id=vertices_2$id, size=vertices_2$size)
#       
#       ged12 = sum(abs(vsizes1$size - vsizes2$size)) + sum(abs(el1$weight - el2$weight))
#       ged_mat[i, j] = ged12
#     }
#   }
# }
# 
# write.xlsx(ged_mat, "../../Output/ged/CoStim_GED_matrix.xlsx")
ged_mat = read.xlsx("../../Output/ged/CoStim_GED_matrix.xlsx", sheet = "Sheet 1")

r = nrow(ged_mat)
c = ncol(ged_mat)

#### creating the symmetric matrix
d1 = ged_mat
d2 = t(ged_mat)
d = matrix(rep(0, r*c), nrow = r)
d[upper.tri(d)] = d1[upper.tri(d1)]
d[lower.tri(d)] = d2[lower.tri(d2)]
isSymmetric(d)
diag(d)

#### non-metric MDS
# Convert the matrix to a distance object
dist_object <- as.dist(d)
# Perform non-metric multidimensional scaling
nmds_result <- isoMDS(dist_object, k = 2)

#### creating data frame with all data
df = data.frame(x=nmds_result$points[,1], y=nmds_result$points[,2])
df = as.data.frame(lapply(df, minMaxNormalize))

df$file_name = unlist(file_name_list)
df$Patient = unlist(lapply(file_name_list, function(x) { strsplit(x, "_")[[1]][7]}))
df$Time = unlist(lapply(file_name_list, function(x) { strsplit(x, "_")[[1]][9]}))
df$Time <- factor(df$Time, levels = c(1, 2, 3), labels = c("Baseline", "After 1 vac", "After 3 vac"))

os_pfs = read.xlsx("../../Data/Raw Data/OS_PFS.xlsx", sheet = "Sheet1")
df = merge(df, os_pfs, by = "Patient", all.x = TRUE)

mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                    brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )

sqr_range_2 = range(min(range(df$x)[1], range(df$y)[1]), max(range(df$x)[2], range(df$y)[2]))

g = ggplot(df)+
  geom_point(aes(x=x, y=y, colour=Patient, shape=Time), cex = 3) +
  #geom_text_repel(aes(x=x, y=y, colour=Patient, label = PFS_cens)) +
  ggforce::geom_mark_ellipse(aes(x = x, y = y, colour = Patient), linetype="dotted", linewidth=0.25) +
  theme(legend.position="right", legend.text=element_text(size=10), legend.title = element_text(size=12),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
  
  scale_colour_manual(values = mycolor) +
  xlim(sqr_range_2) + ylim(sqr_range_2) +
  xlab("MDS dimension 1") + ylab("MDS dimension 2")
plot(g)

svglite("../../Output/ged/CoStim_GED_MDS.svg", width = 7, height = 6)
plot(g)
dev.off()


ged_mat = ged_mat[-c(3, 11, 16, 24, 36, 37), -c(3, 11, 16, 24, 36, 37)]
file_name_list = file_name_list[-c(3, 11, 16, 24, 36, 37)]

r = nrow(ged_mat)
c = ncol(ged_mat)

#### creating the symmetric matrix
d1 = ged_mat
d2 = t(ged_mat)
d = matrix(rep(0, r*c), nrow = r)
d[upper.tri(d)] = d1[upper.tri(d1)]
d[lower.tri(d)] = d2[lower.tri(d2)]
isSymmetric(d)
diag(d)
d = as.dist(d)

group = factor(c(rep("Baseline", 11), rep("After 1 vac", 11), rep("After 3 vac", 11)))
print(adonis2(d ~ group))
