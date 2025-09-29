#' Project Initialization
#'
#' Sets up the project environment by installing/loading required R packages, 
#' setting the working directory, defining a custom color palette, and enabling 
#' Python integration via `reticulate`.

# Loading required libraries
load_lib <- c("deldir", "spatstat", "magrittr", "dplyr", "igraph", "scales", "httr", "tidyverse", 
              "ggnetwork", "ggplot2", "poweRlaw", "imager", "viridis", "plotrix", "openxlsx", "tidyr", 
              "spdep", "maptools", "tmap", "OpenImageR", "dismo", "lctools", "officer", "rvg", "truncnorm", 
              "emdist", "ks", "rlist", "readxl", "OneR", "MASS", "RColorBrewer", "this.path", 
              "causaloptim", "RBGL", "svglite", "ggrepel", "devtools", "geosphere", "proxy", "sampling", 
              "collections", "umap", "gsignal", "e1071", "caTools", "caret", "DescTools", 
              "PerformanceAnalytics", "randomForest", "tsne", "rdist", "Rcpp", "Matrix", "mvtnorm", 
              "gridExtra", "matrixcalc", "msos", "reticulate", "data.table", "ggalt", "intergraph", 
              "kohonen", "rdflib", "cluster", "autoimage", "ggforce", "ggpubr", "stringr","grid", 
              "cowplot", "svglite", "umap", "see", "multcomp", "vegan", "PERMANOVA", "Rtsne")
# Install missing libraries
install_lib <- load_lib[!load_lib %in% installed.packages()]
if(length(install_lib) > 0) { install.packages(install_lib, dependencies = TRUE)}
# Load libraries
sapply(load_lib, require, character=TRUE)

# how to install the introdataviz package to get split and half violin plots
# devtools::install_github("psyteachr/introdataviz")
# devtools::install_github("saeyslab/PeacoQC")

# Additional libraries
library(flowCore)
library(FlowSOM)
library(flowWorkspace)
library(kohonen)
library(flowDensity)

setwd(dirname(this.path::this.path()))

# Python integration
path_to_python <- "C:/Users/AbidaSanjanaShemonti/PycharmProjects/OT-GED/.venv"
reticulate::use_virtualenv(path_to_python, required = TRUE)
reticulate::source_python('C:/Users/AbidaSanjanaShemonti/PycharmProjects/OT-GED/TestOTSinkhorn.py')
reticulate::source_python('C:/Users/AbidaSanjanaShemonti/PycharmProjects/OT-GED/TestGED.py')




#' Project Color Palette
#'
#' Defines a unique custom color palette (`mycolor`) by combining shades from 
#' several RColorBrewer palettes: "Dark2", "Set1", "Accent", and "Set2".
mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                    brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )




#' Min-Max Normalization of a Numeric Vector
#'
#' Scales the input numeric vector to a [0, 1] range by applying min-max normalization.
#' The minimum value of the vector is mapped to 0, and the maximum value is mapped to 1.
#'
#' @param x A numeric vector to be normalized.
#' @param na.rm Logical, indicating whether `NA` values should be removed before 
#'   computing the minimum and maximum. Defaults to `TRUE`.
#'
#' @return A numeric vector of the same length as `x`, with values scaled to [0, 1].
minMaxNormalize <- function(x, na.rm = TRUE) {
  return((x - min(x, na.rm = na.rm)) / (max(x, na.rm = na.rm) - min(x, na.rm = na.rm)))
}




#' Min-Max Normalization with Specified Range
#'
#' Normalizes a numeric vector to the [0, 1] range using specified minimum and maximum values.
#'
#' @param x A numeric vector to be normalized.
#' @param size_min Numeric value representing the minimum for normalization.
#' @param size_max Numeric value representing the maximum for normalization.
#'
#' @return A numeric vector with values scaled to [0, 1] based on the provided min and max.
minmax2 <-function(x, size_min, size_max){
  return( (x- size_min) /(size_max-size_min) )
}




#' Add Cell Population Sizes to Phenotype Layout
#'
#' Adds a `size` column to the phenotype layout data frame, representing the 
#' number of occurrences of each cell type in the sample data.
#'
#' @param sample A data frame containing at least a `cell_type` column.
#' @param pheno_layout A data frame representing the phenotype layout, with a `v_label` column.
#'
#' @return The `pheno_layout` data frame with an additional `size` column indicating 
#'         the count of each cell type in the sample.
add_size_column <- function(sample, pheno_layout) {
  row_ord = pheno_layout$v_label
  
  # Count the occurrences of each cell_type in sample
  cell_type_counts <- as.data.frame(table(sample$cell_type))
  colnames(cell_type_counts) <- c("cell_type", "size")
  
  # Merge the count information with pheno_layout
  pheno_layout <- merge(pheno_layout, cell_type_counts, by.x = "v_label", by.y = "cell_type", all.x = TRUE)
  
  # Replace NA values in the size column with 0
  pheno_layout$size[is.na(pheno_layout$size)] <- 0
  
  # Reorder rows to match the original pheno_layout order
  pheno_layout <- pheno_layout[match(row_ord, pheno_layout$v_label), ]
  
  return(pheno_layout)
}




#' Generate Shades of a Color
#'
#' Generates a set of lighter shades from a given base color.
#'
#' @param color A color name or hex code representing the base color.
#' @param n Integer specifying the number of shades to generate. Defaults to 6.
#'
#' @return A character vector of selected shades of the input color.
generate_shades <- function(color, n = 6) {
  shades <- colorRampPalette(c(color, "white"))(n)
  return(shades[c(5, 3, 1)])
}




#' Extract Legend from a ggplot Object
#'
#' Retrieves the legend grob from a ggplot object.
#'
#' @param myggplot A ggplot object.
#'
#' @return A grob representing the legend of the input ggplot.
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




#' Add Alpha Transparency to Colors
#'
#' Applies an alpha (transparency) value to one or more colors.
#'
#' @param color A color name, hex code, or vector of colors.
#' @param alpha Numeric value between 0 and 1 specifying the transparency. Defaults to 0.5.
#'
#' @return A character vector of colors with the specified alpha transparency applied.
add_alpha <- function(color, alpha = 0.5){
  apply(sapply(color, col2rgb)/255, 2, function(x) 
    rgb(x[1], x[2], x[3], alpha = alpha))
}




#' Compute Phenotype Distance Graph from Cell Population Data
#'
#' Computes a phenotype distance graph based on Hamming distances between cell 
#' populations in a panel. The function supports multiple graph types, including 
#' fully-connected, threshold-based, minimum spanning tree (default), and custom-edge graphs.
#'
#' @param pheno_file Path to an Excel file containing phenotype information. Must have a sheet named "Phenotype".
#' @param graph_kind Integer specifying the type of graph to generate:
#'        1 = fully-connected, 2 = threshold-based, 3 = minimum spanning tree (default), 4 = edge-customized.
#' @param edge_filter_matrix A matrix used to filter edges when `graph_kind = 4`.
#'
#' @return A list containing:
#' \enumerate{
#'   \item The fully-connected phenotype graph object.
#'   \item The filtered phenotype graph object (MST or custom).
#'   \item A data frame of vertices (nodes) with coordinates, labels, colors, and sizes.
#'   \item A data frame of edges with source and target node information.
#'   \item Layout information (x, y coordinates and labels) for the filtered graph.
#'   \item The 0-1 adjacency matrix of the filtered phenotype graph.
#' }
phenotype_distance <- function(pheno_file, graph_kind=3, edge_filter_matrix){
  # Read the Excel file
  data <- read_excel(pheno_file, sheet = "Phenotype")
  print(data)
  
  v_label = unname(unlist(c(data[, 1])))
  
  # Exclude the first column
  data <- data[, -1]
  
  # Convert each row to a comma-separated string
  strings <- apply(data, 1, function(row) paste(row, collapse = ","))
  
  hamming_distance <- function(str1, str2) {
    # Define the custom distances
    token_distances <- list(
      "++" = list("+" = 1, "-" = 2, "dc" = 1, "low" = 1),
      "+" = list("++" = 1, "-" = 1, "dc" = 1, "low" = 1),
      "-" = list("++" = 2, "+" = 1, "dc" = 1, "low" = 1),
      "dc" = list("++" = 1, "+" = 1, "-" = 1, "low" = 1),
      "low" = list("++" = 1, "+" = 1, "-" = 1, "dc" = 1)
    )
    
    # Split the strings into tokens
    substrings1 <- unlist(strsplit(str1, ","))
    substrings2 <- unlist(strsplit(str2, ","))
    
    # Check if the lengths are equal
    if (length(substrings1) != length(substrings2)) {
      stop("Strings must have the same number of tokens")
    }
    
    # Compute the custom Hamming distance
    distance <- 0
    for (i in seq_along(substrings1)) {
      if (substrings1[i] != substrings2[i]) {
        distance <- distance + token_distances[[substrings1[i]]][[substrings2[i]]]
      }
    }
    
    return(distance)
  }
  
  # Compute pairwise distances (upper triangle only)
  n <- length(strings)
  pheno_distance_matrix <- matrix(0, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      pheno_distance_matrix[i, j] <- hamming_distance(strings[i], strings[j])
    }
  }
  
  # Fill the lower triangle of the distance matrix
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      pheno_distance_matrix[i, j] <- pheno_distance_matrix[j, i]
    }
  }
  
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  #### create graph from the phenotype distances 
  pheno_graph_obj = graph_from_adjacency_matrix(pheno_distance_matrix, weighted=TRUE, mode="undirected") %>% 
    set_vertex_attr("label", value = v_label)
  pheno_graph_obj = set_vertex_attr(pheno_graph_obj, "size", value = 6.5)
  V(pheno_graph_obj)$cluster.col = mycolor[1: length(v_label)]
  V(pheno_graph_obj)$cluster.col.alpha = add_alpha(V(pheno_graph_obj)$cluster.col, alpha = 0.5)
  V(pheno_graph_obj)$name = v_label
  
  #### filter the fully connected graph into a tree with the biologically meaningful edges
  #### the default is constructing an MST
  if(graph_kind == 1){
    filtered_graph = pheno_graph_obj
    
  }else if(graph_kind == 2){
    # Calculate the average edge weight
    edge_weights <- E(pheno_graph_obj)$weight
    average_weight <- mean(edge_weights)
    
    # Filter edges with weights below the average
    filtered_edges <- E(pheno_graph_obj)[weight < average_weight]
    
    # Create a new graph with the filtered edges
    filtered_graph <- subgraph.edges(pheno_graph_obj, filtered_edges)
    
  }else if(graph_kind==3){
    filtered_graph = mst(pheno_graph_obj)
    
  }else if(graph_kind == 4){
    # Get the adjacency matrix of the original graph
    adj_matrix <- as_adjacency_matrix(pheno_graph_obj, attr = "weight", sparse = FALSE)
    
    # Apply the edge filter matrix to the adjacency matrix
    filtered_adj_matrix <- adj_matrix * as.matrix(edge_filter_matrix)
    
    # Create a new graph from the filtered adjacency matrix
    filtered_graph <- graph_from_adjacency_matrix(filtered_adj_matrix, mode = "undirected", weighted = TRUE)
    # Copy vertex attributes from the original graph to the new graph
    vertex_attr_names <- vertex_attr_names(pheno_graph_obj)
    for (attr in vertex_attr_names) {
      vertex_attr(filtered_graph, attr) <- vertex_attr(pheno_graph_obj, attr)
    }
    
  }else{
    filtered_graph = mst(pheno_graph_obj)
    
  }
  
  mst_adj_matrix <- as_adjacency_matrix(filtered_graph, sparse = FALSE)
  # Convert the adjacency matrix to binary (0 and 1)
  mst_adj_matrix[mst_adj_matrix > 0] <- 1
  # Print the adjacency matrix
  print(mst_adj_matrix)
  
  # Apply custom layout
  set.seed(42)
  # graph_lo = layout_with_dh(filtered_graph, weight.node.dist = 0.2, weight.edge.lengths = 0.8)
  graph_lo = layout_with_kk(filtered_graph, weights = E(filtered_graph)$weight) # with layout_with_kk larger edge weights will result longer edges
  rownames(graph_lo) <- seq_len(nrow(graph_lo))
  
  # Convert to data frames
  vertices <- data.frame(
    id = V(filtered_graph)$name,
    x = graph_lo[,1],
    y = graph_lo[,2],
    label = V(filtered_graph)$label,
    color = V(filtered_graph)$cluster.col,
    size = V(filtered_graph)$size
  )
  
  edges <- igraph::as_data_frame(filtered_graph, what = "edges")
  edges <- merge(edges, vertices, by.x = "from", by.y = "id")
  edges <- merge(edges, vertices, by.x = "to", by.y = "id", suffixes = c(".from", ".to"))
  
  g_plot = ggplot() +
    geom_curve(data = edges, aes(x = x.from, y = y.from, xend = x.to, yend = y.to), 
               linewidth = 0.5, 
               color = "grey", curvature = 0.1, ncp = 100, show.legend = F) +
    
    geom_point(data = vertices, aes(x = x, y = y, size=size), fill = "wheat", color = "grey", pch=21,
               show.legend = FALSE) +
    
    geom_text_repel(data = vertices, aes(x = x, y = y, label = label), color = "black",
                    size = 3.5, max.overlaps = Inf) +
    
    scale_size_continuous(range = c(min(vertices$size), max(vertices$size))) +
    scale_color_identity() +
    theme_void() +
    theme(legend.position.inside = "right")+
    theme(legend.position.inside = c(1, 0.1), legend.direction = "vertical",
          plot.title = element_text(hjust = 0.5, size=16),
          plot.subtitle = element_text(hjust = 0.5, size=10))
  #labs(title="Phenotype distance graph")
  print(g_plot)
  
  layout_info = data.frame(x=graph_lo[, 1], y=graph_lo[, 2], v_label=v_label)
  
  #### returns respectively: the graph object for the fully-connected phenotype graph,
  #### the graph object for the filtered phenotype graph (mst or custom)
  #### the data structures for the nodes and edges
  #### the layout of that filtered phenotype graph (x and y coordinates of the nodes) and label together, for convenience
  #### the 0-1 adjacency matrix of the filtered phenotype graph
  return(list(pheno_graph_obj, filtered_graph, vertices, edges, layout_info, mst_adj_matrix))
}




#' Compute Sinkhorn Distances Between Clusters
#'
#' Calculates pairwise Sinkhorn distances between clusters in a dataset based 
#' on selected features, excluding time, FSC, and SSC measurements. Cluster 
#' labels are mapped to numeric values for computation.
#'
#' @param data_sample A data frame containing sample data, with the last column 
#'                    representing cluster labels.
#' @param label_map A data frame mapping cluster names to numeric labels.
#'
#' @return A list containing:
#' \enumerate{
#'   \item `sink_dist_mat`: A matrix of pairwise Sinkhorn distances between clusters.
#'   \item `data_sample`: The input data frame with cluster labels mapped to numeric values (excluding the numeric label column).
#' }
computeSinkhornDistances <- function(data_sample, label_map){
  #### extract how many features are there; note that time, fsc, ssc features should not be in here
  feature_count = ncol(data_sample) - 1
  clust_count = length(unique(data_sample[, ncol(data_sample)]))
  
  #### convert the cluster names to numericals
  col_ord = c(colnames(data_sample), colnames(label_map)[2])
  data_sample[, ncol(data_sample)] = as.factor(data_sample[, ncol(data_sample)])
  data_sample <- merge(data_sample, label_map, by = "cell_type", all.x = TRUE)
  data_sample = data_sample[, col_ord]
  
  # #### save the sample data
  # cat("Saving data file...\n")
  # write.csv(data_sample[, c(1:feature_count, ncol(data_sample))], paste("./PycharmProjects/SOT/", "data.csv", sep=""), row.names=FALSE)
  # 
  setDT(data_sample)  # Convert the dataframe to a data.table
  s_sample = data_sample # just another variable name, was something different before
  
  # Convert the result back to a data.table
  s_sample$numeric_label = as.numeric(s_sample$numeric_label)
  setDT(s_sample)
  table(s_sample$numeric_label)
  
  
  #### Sinkhorn test
  setkey(s_sample, "numeric_label")
  
  sink_dist_mat = matrix(0, clust_count, clust_count, byrow = TRUE)
  for(i in c(1:clust_count)){
    for(j in c(i:clust_count)){
      if(i != j){
        cat(i, j , "\n")
        source = s_sample[J(i)] # s_sample[s_sample$cluster.labels==i, ]
        target = s_sample[J(j)] # s_sample[s_sample$cluster.labels==j, ]
        
        # Set seed for reproducibility
        # set.seed(123)
        if(nrow(source) > 10000){source = source[sample(.N, 10000)]}
        if(nrow(target) > 10000){target = target[sample(.N, 10000)]}
        
        n = nrow(source)
        m = nrow(target)
        a = matrix(rep(1/n, n)) # initial uniform source mass
        b = matrix(rep(1/m, m)) # initial uniform target mass
        computed_cost_mat = proxy::dist(source[, 1:feature_count], target[, 1:feature_count], method = "Euclidean")
        # cdist(source[, 1:feature_count], target[, 1:feature_count])
        
        lambd = 0.1
        res = testSinkhorn(a, b, computed_cost_mat, lambd)
        sinkhorn_dist = sum(res*computed_cost_mat)
        
        sink_dist_mat[i, j] = sinkhorn_dist
      }
    }
  }
  
  data_sample = as.data.frame(data_sample)
  return(list(sink_dist_mat, data_sample[, -ncol(data_sample)]))
}




#' Visualize Filtered Graph Based on Sinkhorn Distances
#'
#' Generates a filtered network graph from a Sinkhorn distance matrix, applying 
#' a phenotype-based filter to edges. Computes community structure and overlays 
#' vertices and edges on a layout, producing a ggplot visualization.
#'
#' @param sink_dist_mat A matrix of pairwise Sinkhorn distances between clusters.
#' @param sample A data frame containing sample data with cluster labels.
#' @param graph_lo A data frame containing x and y coordinates and node labels for the graph layout.
#' @param pheno_filter_matrix A 0-1 matrix used to filter edges in the graph.
#' @param file A string used for labeling or referencing the output (optional).
#'
#' @return A list containing:
#' \enumerate{
#'   \item `graph_obj`: The filtered igraph object.
#'   \item `vertices`: A data frame of vertex information, including coordinates, size, color, and community.
#'   \item `edges`: A data frame of edge information with source and target node coordinates and weights.
#' }
visualizeFilteredGraph <- function(sink_dist_mat, sample, graph_lo, pheno_filter_matrix, file){
  sink_dist_mat[is.na(sink_dist_mat)] <- 1e-10  # Replace NA with very small value
  
  r = nrow(sink_dist_mat)
  c = ncol(sink_dist_mat)
  
  #### creating the symmetric matrix
  d1 = sink_dist_mat
  d2 = t(sink_dist_mat)
  d = matrix(rep(0, r*c), nrow = r)
  d[upper.tri(d)] = d1[upper.tri(d1)]
  d[lower.tri(d)] = d2[lower.tri(d2)]
  isSymmetric(d)
  diag(d)
  
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  #### create graph from the Sinkhorn measures 
  graph_obj = graph_from_adjacency_matrix(d, weighted=TRUE, mode="undirected") %>% 
    set_vertex_attr("label", value = graph_lo$v_label)
  graph_obj = set_vertex_attr(graph_obj, "size", value = (graph_lo$size/sum(graph_lo$size)*100))
  V(graph_obj)$cluster.col = mycolor[1: length(graph_lo$v_label)]
  V(graph_obj)$cluster.col.alpha = add_alpha(V(graph_obj)$cluster.col, alpha = 0.5)
  V(graph_obj)$name = graph_lo$v_label
  
  #### graph filtering step
  # Get the adjacency matrix of the original graph
  original_adj_matrix <- as_adjacency_matrix(graph_obj, attr = "weight", sparse = FALSE)
  
  # Multiply the original adjacency matrix with the reordered MST adjacency matrix
  new_adj_matrix <- original_adj_matrix * pheno_filter_matrix
  
  # Create a new graph from the updated adjacency matrix
  new_graph <- graph_from_adjacency_matrix(new_adj_matrix, weighted = TRUE, mode = "undirected")
  
  # Copy vertex attributes from the original graph to the new graph
  for (attr in vertex_attr_names(graph_obj)) {
    new_graph <- set_vertex_attr(new_graph, attr, value = vertex_attr(graph_obj, attr))
  }
  graph_obj = new_graph
  
  # Detect communities
  community <- cluster_fast_greedy(graph_obj, weights = exp(-E(graph_obj)$weight))
  vertex_community_map = data.frame(vertex_label = V(graph_obj)$label, community_id=membership(community))
  vertex_community_map <- vertex_community_map %>%
    group_by(community_id) %>%
    mutate(community_name = first(vertex_label)) %>%
    ungroup()%>%
    mutate(community_name = word(community_name, 1))
  V(graph_obj)$community <- vertex_community_map$community_id
  print(membership(community))
  
  # Convert to data frames
  vertices <- data.frame(
    id = V(graph_obj)$name,
    x = graph_lo$x,
    y = graph_lo$y,
    label = V(graph_obj)$label,
    color = V(graph_obj)$cluster.col,
    size = V(graph_obj)$size,
    community = V(graph_obj)$community
  )
  print(vertices)
  
  edges <- igraph::as_data_frame(graph_obj, what = "edges")
  edges <- merge(edges, vertices, by.x = "from", by.y = "id")
  edges <- merge(edges, vertices, by.x = "to", by.y = "id", suffixes = c(".from", ".to"))
  
  #### ploting the filtered graph
  plot(graph_obj, layout=as.matrix(graph_lo[, 2:3]))
  
  g_plot = ggplot() +
    geom_curve(data = edges, aes(x = x.from, y = y.from, xend = x.to, yend = y.to, color=weight),
               linewidth = 1, curvature = 0.2, ncp = 100) +
    scale_color_gradient(low = "turquoise", high = "salmon", name = "Edge Weight",
                         guide = guide_colorbar(direction = "vertical",
                                                title.position = "top",
                                                title.theme = element_text(size = 8),
                                                title.vjust = 0.75)) +
    
    geom_circle(data = vertices, aes(x0 = x, y0 = y, r=minmax2(size, min(vertices$size), max(vertices$size))), # minmax2(size, min(vertices$size), max(vertices$size))
                fill = "wheat", color = "wheat3", show.legend = FALSE) +
    
    geom_text_repel(data = vertices, aes(x = x, y = y, label = paste(str_wrap(label, width = 15), "\n", round(size,2), "%", sep="")),
                    color = "black", size = 2, max.overlaps = Inf) +
    
    coord_fixed() +
    theme_void() +
    theme(legend.direction = "vertical", legend.position = "right",
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          plot.margin = unit(c(0, 0, 0, 0), "mm"))
  #labs(title = paste("Graph-based visualization", file, "\n", sep = "\n"))
  print(g_plot)
  
  return(list(graph_obj, vertices, edges))
}




#' Visualize Filtered Graph (Alternative Layout)
#'
#' Plots a network graph with vertices and edges, where edge widths are 
#' inversely proportional to weights and edge colors represent weights. Node 
#' sizes and labels are visualized with circles and text annotations.
#'
#' @param vertices A data frame containing vertex information, including coordinates, size, and labels.
#' @param edges A data frame containing edge information with source and target coordinates and weights.
#' @param range A numeric vector of length 2 specifying the range of edge weights for the color gradient.
#'
#' @return None. The function prints the ggplot visualization of the graph.
visualizeFilteredGraph_2 <- function(vertices, edges, range){
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  g_plot = ggplot() +
    geom_curve(data = edges, aes(x = x.from, y = y.from, xend = x.to, yend = y.to, color=weight, size = 1/(weight)), # linewidth = 1
               curvature = 0.2, ncp = 100, show.legend = c(size = FALSE)) +
    #scale_size_continuous(range = range(10/edges$weight)) +
    scale_color_gradient(low = "blue", high = "#FFC20A", name = "Edge Weight",
                         limits = c(range[1], range[2]),  
                         guide = guide_colorbar(direction = "vertical",
                                                title.position = "top",
                                                title.theme = element_text(size = 8),
                                                title.vjust = 0.75)) +
    
    geom_circle(data = vertices, aes(x0 = x, y0 = y, r=size/200), # minmax2(size, min(vertices$size), max(vertices$size))
                fill = "wheat", color = "black", show.legend = FALSE) +
    
    geom_text_repel(data = vertices, aes(x = x, y = y, label = paste(str_wrap(label, width = 15), "\n", round(size,2), "%", sep="")),
                    color = "black", size = 10, max.overlaps = Inf) +
    
    coord_fixed() +
    theme_void() +
    theme(legend.direction = "vertical", legend.position = "right",
          legend.text = element_text(size = 12), legend.title = element_text(size = 14),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          plot.margin = unit(c(0, 0, 0, 0), "mm"))
  #labs(title = paste("Graph-based visualization", file, "\n", sep = "\n"))
  print(g_plot)
}




#' Overlay Filtered Graphs for Comparative Visualization
#'
#' Compares two network graphs (baseline vs. modified) by overlaying vertices and edges. 
#' Highlights added, removed, and changed nodes and edges using colors, linetypes, and labels.
#'
#' @param vertices_b Data frame of baseline graph vertices with coordinates, size, and labels.
#' @param edges_b Data frame of baseline graph edges with source/target coordinates and weights.
#' @param vertices_t Data frame of target/modified graph vertices.
#' @param edges_t Data frame of target/modified graph edges.
#' @param graph_lo Optional layout data frame (x/y coordinates).
#' @param file_b Filename string for baseline graph (used to infer timepoint).
#' @param file_t Filename string for target graph (used to infer timepoint).
#' @param range Numeric vector specifying edge weight range for visualization scaling.
#'
#' @return None. Prints a ggplot visualization highlighting changes between the two graphs.
#'         Added nodes/edges, removed nodes/edges, and modified nodes/edges are visually distinguished.
overlayFilteredGraph <- function(vertices_b, edges_b, vertices_t, edges_t, graph_lo=NULL, file_b, file_t, range){
  if(grepl("VAC_1", file_b)){bottom = "Baseline"}
  else if(grepl("VAC_2", file_b)){bottom = "After 1 vac"}
  else if(grepl("VAC_3", file_b)){bottom = "After 3 vac"}
  else{bottom = "Unknown"}
  
  if(grepl("VAC_1", file_t)){top = "Baseline"}
  else if(grepl("VAC_2", file_t)){top = "After 1 vac"}
  else if(grepl("VAC_3", file_t)){top = "After 3 vac"}
  else{top = "Unknown"}
  
  linetypes <- setNames(c("dotted", "solid"), c(bottom, top))
  linetype_levels <- c("Baseline", "After 1 vac", "After 3 vac")
  
  panel = strsplit(strsplit(file_b, "_Tube")[[1]][1], "Norm_Norm_")[[1]][2]
  patient = strsplit(strsplit(file_b, "_VAC")[[1]][1], "_")[[1]][7]
  file = paste0("Panel: ", panel, ", Patient ID: ", patient)
  
  vertices_b$size = round(vertices_b$size, 2)
  vertices_t$size = round(vertices_t$size, 2)
  #edges_b$weight=round(edges_b$weight, 2)
  #edges_t$weight=round(edges_t$weight, 2)
  
  size_min = min(min(vertices_b$size), min(vertices_t$size))
  size_max = max(max(vertices_b$size), max(vertices_t$size))
  
  # Identify added, removed, and changed vertices
  added_vertices <- setdiff(vertices_t$id, vertices_b$id)
  removed_vertices <- setdiff(vertices_b$id, vertices_t$id)
  common_vertices <- intersect(vertices_b$id, vertices_t$id)
  
  # Create unique identifiers for edges with sorted nodes
  edges_b$edge_id <- apply(edges_b[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
  edges_t$edge_id <- apply(edges_t[, c("from", "to")], 1, function(x) paste(sort(x), collapse = "_"))
  
  # Identify added, removed, and changed edges
  added_edges <- edges_t[!edges_t$edge_id %in% edges_b$edge_id, ]
  removed_edges <- edges_b[!edges_b$edge_id %in% edges_t$edge_id, ]
  common_edges <- merge(edges_b, edges_t, by = "edge_id", suffixes = c(".base", ".mod"))
  
  # Identify edges with changed weights
  changed_edges <- common_edges[common_edges$weight.base != common_edges$weight.mod, ]
  
  # Create data frames for added and removed vertices
  vertices_added <- vertices_t[vertices_t$id %in% added_vertices, ]
  vertices_removed <- vertices_b[vertices_b$id %in% removed_vertices, ]
  
  # Create data frames for added and removed edges
  edges_added <- added_edges
  edges_removed <- removed_edges
  
  # Identify changed vertices (location or size)
  changed_vertices <- vertices_b[vertices_b$id %in% common_vertices, ]
  changed_vertices <- merge(changed_vertices, vertices_t, by = "id", suffixes = c(".base", ".mod"))
  changed_vertices <- changed_vertices[changed_vertices$x.base != changed_vertices$x.mod | 
                                         changed_vertices$y.base != changed_vertices$y.mod | 
                                         changed_vertices$size.base != changed_vertices$size.mod, ]
  changed_vertices$label.new <- paste(changed_vertices$label.mod,
                                      paste(ifelse((changed_vertices$size.base - changed_vertices$size.mod) >= 0, "(- by ", "(+ by "), 
                                            round(abs(changed_vertices$size.base - changed_vertices$size.mod),2), "%)", sep=""),
                                      sep = "\n")
  
  # Combine the IDs of added, removed, and changed vertices
  excluded_ids <- unique(c(vertices_added$id, vertices_removed$id, changed_vertices$id))
  # Filter vertices_b to exclude these IDs
  vertices_b_filtered <- vertices_b[!vertices_b$id %in% excluded_ids, ]
  
  g_plot = ggplot() +
    # base vertices
    geom_circle(data=vertices_b, aes(x0=x, y0=y, r=size/100, linetype=factor(bottom, levels = linetype_levels)), fill="grey", color="grey", alpha=0.3, size=1) +
    
    geom_curve(data = edges_added, aes(x = x.from, y = y.from, xend = x.to, yend = y.to, color = weight),
               linewidth = 5,curvature = 0.2, ncp = 100) +
    geom_curve(data = edges_removed, aes(x = x.from, y = y.from, xend = x.to, yend = y.to),
               linewidth = 5, linetype="dotted", color = "red3", curvature = 0.2, ncp = 100) +
    geom_curve(data = changed_edges, aes(x = x.from.base, y = y.from.base, xend = x.to.base, yend = y.to.base, color=(weight.base-weight.mod), size = 1/exp(weight.base-weight.mod)), #linewidth = 5, 
               curvature = 0.2, ncp = 100, show.legend = c(size=F)) +
    scale_size_continuous(range = 2*range(1/exp(changed_edges$weight.base - changed_edges$weight.mod))) +
    
    # added vertices if any
    geom_circle(data = vertices_added, aes(x0 = x, y0 = y, r = size/100, linetype=factor(top, levels = linetype_levels)), 
                fill = NA, color = "salmon") +
    geom_text_repel(data = vertices_added, aes(x = x, y = y, label = paste(str_wrap(15), "\n", size, "%", sep="")),
                    color="black", size =2.5, max.overlaps = Inf) +
    # removed vertices
    geom_circle(data = vertices_removed, aes(x0 = x, y0 = y, r = size/100, linetype=factor(top, levels = linetype_levels)), 
                fill = "snow", color = "red3") +
    geom_text_repel(data = vertices_removed, aes(x = x, y = y, label = str_wrap(label, width = 15)), color="black",
                    size = 2.5, max.overlaps = Inf) +
    # changed vertices
    geom_circle(data = changed_vertices, aes(x0 = x.mod, y0 = y.mod, r = size.mod/100, linetype=factor(top, levels = linetype_levels)), 
                fill = "wheat", color = "black") +
    geom_text_repel(data = changed_vertices, aes(x = x.mod, y = y.mod, label = label.new), # paste(str_wrap(label.mod, width = 15), "\n", size.mod, "%", sep="")),
                    color="black", size = 10, max.overlaps = Inf) +
    # base vertices outline
    geom_circle(data=vertices_b, aes(x0=x, y0=y, r=size/100, linetype=factor(bottom, levels = linetype_levels)), color="black", size=1) +
    
    scale_color_gradient(low = "orange", high = "black", name = "Edge Weight (difference)",
                         #limits = c(-range[2], range[2]),  
                         guide = guide_colorbar(direction = "vertical",
                                                title.position = "top",
                                                title.theme = element_text(size = 8),
                                                title.vjust = 0.75)) +
    scale_linetype_manual(values = linetypes, 
                          name = "Time", guide = guide_legend(direction = "vertical",
                                                              title.position = "top",
                                                              title.theme = element_text(size = 8))) +
    coord_fixed() +
    theme_void() +
    theme(legend.direction = "vertical", legend.position = "right",
          #legend.box = "vertical",  legend.box.just = "left", 
          legend.text = element_text(size = 9),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10),
          plot.margin = unit(c(0, 0, 0, 0), "mm")) 
    #labs(title = "Cell Population Comparison", subtitle = file)
  
  # table <- ggtexttable(changed_vertices$label.new, rows = NULL, 
  #                      theme = ttheme(base_size = 10))
  # ggarrange(g_plot, table, ncol = 2, nrow = 1, widths = c(3,2))
  print(g_plot)
}




#' Scale Edge Weights Across Multiple Graphs
#'
#' Normalizes the edge weights of three edge data frames to a common 0b
