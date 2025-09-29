#' Dimensionality Reduction and Visualization of Flow Cytometry Data
#'
#' This script performs dimensionality reduction and visualization for flow cytometry 
#' data from a Co-inhibition panel. Specifically, it:
#' 
#' 1. Loads pre-processed, cell-assigned data from an Excel file.
#' 2. Extracts marker expression columns while excluding FSC and SSC channels.
#' 3. Performs PCA on the marker data and constructs a data frame with the first two PCs.
#' 4. Removes duplicate marker rows and runs t-SNE on the deduplicated data to capture 
#'    non-linear structures in two dimensions.
#' 5. Runs UMAP on the full marker data to generate a 2D embedding preserving local 
#'    and global structure.
#' 6. Generates scatter plots for PCA, t-SNE, and UMAP embeddings, coloring points 
#'    by cell type and customizing plot aesthetics.
#' 7. Saves the resulting plots as high-resolution PNG files in the output directory.
#'
#' @note The script assumes the last column in the input data contains cell type labels 
#'       and that the marker columns are located in positions 8 to 21. The color palette 
#'       for plotting should be defined in the variable `mycolor`.


setwd(dirname(this.path::this.path()))
source("MPMUtilities.R")
setwd(dirname(this.path::this.path()))


# Load your data
data_file = "../Data/Raw Data/Co-inhibition/Final Cell Assigned Data/TP1/Norm_Norm_coinhib_g1_Tube_013_MCV005_VAC_1_QC_Lymph_CellType.xlsx"
data = read.xlsx(data_file) # Co-inhibition patient MCV005 vac_1, replace with other file or path if needed

# Assume the last column is the cell type label, and filtering out the FSC and SSC columns
markers <- data[, 8:21]
cell_types <- data[, ncol(data)]

# Remove duplicate rows
unique_markers <- unique(markers)
unique_cell_types <- cell_types[!duplicated(markers)]


# PCA
pca_result <- prcomp(unique_markers, scale. = TRUE)
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], CellType = unique_cell_types)

# Run t-SNE on the de-duplicated data
set.seed(42)
tsne_result <- Rtsne(unique_markers, dims = 2, perplexity = 30)
tsne_df <- data.frame(tSNE1 = tsne_result$Y[,1], tSNE2 = tsne_result$Y[,2], CellType = unique_cell_types)

# UMAP
umap_result <- umap(unique_markers)
umap_df <- data.frame(UMAP1 = umap_result$layout[,1], UMAP2 = umap_result$layout[,2], CellType = unique_cell_types)


# Plot PCA
g1 = ggplot(pca_df, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(alpha = 0.7, cex = 0.3) +
  theme(legend.position="none", legend.text=element_text(size=10), legend.title = element_text(size=11),
        legend.key=element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "white", linewidth=0.25, linetype=2)) +
  scale_colour_manual(values = mycolor) +
  guides(color = guide_legend(nrow=2, override.aes = list(size = 4))) +
  ggtitle("PCA of FC data")
ggsave("../Output/pca_example.png", plot = g1, width = 2.5, height = 2.5, dpi = 1200)


# Plot t-SNE
g2 = ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = CellType)) +
  geom_point(alpha = 0.7, cex = 0.3) +
  theme(legend.position="none", legend.text=element_text(size=21), legend.title = element_text(size=28),
        legend.key=element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "white", linewidth=0.25, linetype=2)) +  
  scale_colour_manual(values = mycolor) +
  ggtitle("t-SNE of FC data")
ggsave("../Output/tsne_example.png", plot = g2, width = 2.5, height = 2.5, dpi = 1200)


# Plot UMAP
g3 = ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(alpha = 0.7, cex = 0.3) +
  theme(legend.position="none", legend.text=element_text(size=21), legend.title = element_text(size=28),
        legend.key=element_blank(),
        #legend.box.margin=margin(-10,-10,-10,-10),
        plot.title = element_text(hjust = 0.5, size=10),
        plot.subtitle = element_text(hjust = 0.5, size=10),
        axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
        panel.background = element_rect(fill='white', colour='black'),
        panel.grid.major = element_line(color = "white", linewidth=0.25, linetype=2)) + 
  scale_colour_manual(values = mycolor) +
  ggtitle("UMAP of FC data")
ggsave("../Output/umap_example.png", plot = g3, width = 2.5, height = 2.5, dpi = 1200)
