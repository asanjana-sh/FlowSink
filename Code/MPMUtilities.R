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
#' Normalizes the edge weights of three edge data frames to a common 0–1 range
#' based on the combined minimum and maximum weights across all three datasets.
#'
#' @param edges_1 A data frame of edges with a `weight` column (first graph).
#' @param edges_2 A data frame of edges with a `weight` column (second graph).
#' @param edges_3 A data frame of edges with a `weight` column (third graph).
#'
#' @details 
#' The function computes the global minimum and maximum weight across all three 
#' edge data frames, then linearly scales the weights in each data frame to the 
#' 0–1 range using the `minmax2` function.
#'
#' @return A list of three data frames (`edges_1`, `edges_2`, `edges_3`) with 
#'         normalized `weight` columns.
scaleEdgeWeights <- function(edges_1, edges_2, edges_3){
  w_range = range(c(edges_1$weight, edges_2$weight, edges_3$weight))
  edges_1$weight = minmax2(edges_1$weight, w_range[1], w_range[2])
  edges_2$weight = minmax2(edges_2$weight, w_range[1], w_range[2])
  edges_3$weight = minmax2(edges_3$weight, w_range[1], w_range[2])
  
  return( list(edges_1, edges_2, edges_3) )
}




#' Visualize UMAP Embedding of Flow Cytometry Data
#'
#' Computes a 2D UMAP embedding of the input dataset and generates a scatter plot
#' colored by cell type.
#'
#' @param data_sample A numeric data frame or matrix where rows are cells and
#'                    columns are marker measurements. The last column is ignored
#'                    during UMAP computation.
#' @param cell_type A vector of cell type labels corresponding to each row of `data_sample`.
#' @param file A string representing the file or sample name to include in the plot title.
#'
#' @details
#' Uses the `umap` function with 30 nearest neighbors and a minimum distance of 0.3
#' to generate the 2D embedding. The points are plotted using `ggplot2` with distinct
#' colors per cell type. The plot limits are set to a square range to preserve aspect ratio.
#'
#' @return Generates a ggplot scatter plot of the UMAP embedding; does not return a value.
visualizeUMAP <- function(data_sample, cell_type, file){
  config <- umap.defaults
  config$n_neighbors <- 30
  config$min_dist <- 0.3
  
  umap_result <- umap(data_sample[, -ncol(data_sample)], config = config)
  umap_data <- as.data.frame(umap_result$layout)
  umap_data$class_labels <- cell_type
  
  mycolor <- unique(c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")))
  
  sqr_range <- range(c(umap_data$V1, umap_data$V2))
  
  g <- ggplot(umap_data, aes(x = V1, y = V2, color = class_labels)) +
    geom_point(pch = 19, cex = 0.8) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      panel.background = element_rect(fill = 'white', colour = 'black'),
      panel.grid.major = element_line(color = "grey", linewidth = 0.25, linetype = 2)
    ) +
    scale_colour_manual(values = mycolor) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    xlim(sqr_range) + ylim(sqr_range) +
    xlab("UMAP dimension 1") + ylab("UMAP dimension 2") +
    labs(title = paste("UMAP of FC Data", file, sep = "\n"))
  
  plot(g)
}




#' Compare GED Distribution Across Time Points
#'
#' Takes a data frame of GED (Graph Edit Distance) values for different time point pairs,
#' reshapes it into long format, and visualizes the distribution using half-violin plots,
#' boxplots, and individual points.
#'
#' @param ged A data frame where each column corresponds to a specific time point pair 
#'            and contains GED values.
#'
#' @return A ggplot object showing the distribution of GED values across time point pairs.
compareGEDDistribution <- function(ged){
  # First, you need to reshape the dataframe to a long format
  ged_long <- tidyr::pivot_longer(ged, cols = everything(), names_to = "TP_pair", values_to = "GED")
  
  # Specify the order of the time variable
  ged_long$TP_pair <- factor(ged_long$TP_pair, 
                             levels = c("Baseline_After1Vac", "After1Vac_After3Vac", "Baseline_After3Vac"))  # Adjust the levels as per your data
  
  
  # Create the boxplot
  print(ggplot(ged_long, aes(x = TP_pair, y = GED, fill=TP_pair)) +
    geom_violinhalf(trim=FALSE, show.legend = FALSE, alpha=0.2) +
    geom_boxplot(color = "black", width = 0.2, outlier.size = 0.08, size = 0.4,
                 staplewidth = 0, outliers = FALSE, show.legend = FALSE) +
    geom_point(color="black", size=0.6, alpha=0.9, show.legend = FALSE) +
    theme(legend.position = "right", legend.text = element_text(size = 7), legend.title = element_blank(),
          legend.key=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          axis.text.x = element_text(size = 7, angle=0), axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = 2)) +
    labs(title = "", x = "Time Point Pair", y = "GED"))
  
}




#' Overlay Filtered Graphs for AML Comparison
#'
#' Compares two graphs representing Normal and AML conditions by overlaying their
#' vertices and edges. The function highlights added, removed, and changed vertices 
#' and edges, with visual cues such as color, size, and line type to indicate differences.
#'
#' @param vertices_b Data frame of vertices for the baseline (Normal) graph.
#' @param edges_b Data frame of edges for the baseline (Normal) graph.
#' @param vertices_t Data frame of vertices for the target (AML) graph.
#' @param edges_t Data frame of edges for the target (AML) graph.
#' @param graph_lo Optional graph object for additional plotting reference (default NULL).
#' @param file_b Filename or identifier for the baseline graph.
#' @param file_t Filename or identifier for the target graph.
#' @param range Numeric range used for scaling edge weight differences.
#'
#' @return A ggplot object visualizing the overlay of baseline and AML graphs, 
#'         showing added, removed, and changed vertices and edges.
overlayFilteredGraph_AML <- function(vertices_b, edges_b, vertices_t, edges_t, graph_lo=NULL, file_b, file_t, range){
  bottom = "Normal"
  top = "AML"
  
  linetypes <- setNames(c("dotted", "solid"), c(bottom, top))
  linetype_levels <- c("Normal", "AML")
  
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
    geom_circle(data=vertices_b, aes(x0=x, y0=y, r=size/200, linetype=factor(bottom, levels = linetype_levels)), fill="grey", color="grey", alpha=0.3, size=1) +
    
    geom_curve(data = edges_added, aes(x = x.from, y = y.from, xend = x.to, yend = y.to, color = weight),
               linewidth = 5,curvature = 0.2, ncp = 100) +
    geom_curve(data = edges_removed, aes(x = x.from, y = y.from, xend = x.to, yend = y.to),
               linewidth = 5, linetype="dotted", color = "red3", curvature = 0.2, ncp = 100) +
    geom_curve(data = changed_edges, aes(x = x.from.base, y = y.from.base, xend = x.to.base, yend = y.to.base, color=(weight.base-weight.mod), size = 1/(weight.base-weight.mod)), #linewidth = 5, 
               curvature = 0.2, ncp = 100, show.legend = c(size=F)) +
    #scale_size_continuous(range = 2*range(1/exp(changed_edges$weight.base - changed_edges$weight.mod))) +
    
    # added vertices if any
    geom_circle(data = vertices_added, aes(x0 = x, y0 = y, r = size/200, linetype=factor(top, levels = linetype_levels)), 
                fill = NA, color = "salmon") +
    geom_text_repel(data = vertices_added, aes(x = x, y = y, label = paste(str_wrap(15), "\n", size, "%", sep="")),
                    color="black", size =2.5, max.overlaps = Inf) +
    # removed vertices
    geom_circle(data = vertices_removed, aes(x0 = x, y0 = y, r = size/200, linetype=factor(top, levels = linetype_levels)), 
                fill = "snow", color = "red3") +
    geom_text_repel(data = vertices_removed, aes(x = x, y = y, label = str_wrap(label, width = 15)), color="black",
                    size = 2.5, max.overlaps = Inf) +
    # changed vertices
    geom_circle(data = changed_vertices, aes(x0 = x.mod, y0 = y.mod, r = size.mod/200, linetype=factor(top, levels = linetype_levels)), 
                fill = "wheat", color = "black") +
    geom_text_repel(data = changed_vertices, aes(x = x.mod, y = y.mod, label = label.new), # paste(str_wrap(label.mod, width = 15), "\n", size.mod, "%", sep="")),
                    color="black", size = 4, max.overlaps = Inf) +
    # base vertices outline
    geom_circle(data=vertices_b, aes(x0=x, y0=y, r=size/200, linetype=factor(bottom, levels = linetype_levels)), color="black", size=1) +
    
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




#' Create Pie Chart Grob Without Legend
#'
#' Generates a pie chart grob from class data for use in composite graphs, 
#' without displaying a legend.
#'
#' @param class_data A data frame containing the features and values for the pie chart. 
#'                   The first four columns are excluded from the chart.
#' @param mycolor A vector of colors to use for the pie chart slices.
#'
#' @return A grob object representing the pie chart, suitable for embedding in other plots.
create_pie_grob <- function(class_data, mycolor) {
  pie_data <- data.frame(
    feature = names(class_data)[-c(1:4)],
    percentage = as.numeric(class_data[-c(1:4)])
  )
  
  pie_grob <- ggplot(pie_data, aes(x = "", y = percentage, fill = feature)) +
    geom_bar(stat = "identity", width = 1, color = "black", size = 0.15) +
    coord_polar("y") +
    scale_fill_manual(values = mycolor) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 8)) +
    ggtitle(class_data$v_label)  # Add cell_type as the title
  
  ggplotGrob(pie_grob)
}




#' Create Individual Pie Chart Grob with Legend
#'
#' Generates a pie chart grob for a single class with a legend, and also prints 
#' a corresponding bar chart of marker percentages.
#'
#' @param class_data A data frame containing feature values and metadata. The first four columns are excluded from the chart.
#' @param mycolor A vector of colors to use for the pie chart slices.
#'
#' @return A grob object representing the pie chart with a legend. The bar chart is printed to the current graphics device.
create_individual_pie_grob <- function(class_data, mycolor) {
  pie_data <- data.frame(Marker = names(class_data)[-c(1:4)],
                         percentage = round(as.numeric(class_data[-c(1:4)]) / class_data$size, 2)*100 )
  
  pie_grob <- ggplot(pie_data, aes(x = "", y = percentage, fill = Marker)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 0.15, alpha=0.8) +
    coord_polar("y") +
    scale_fill_manual(values = mycolor) +
    theme_void() +
    theme(legend.direction = "vertical", legend.position = "right",
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size=10),
          plot.subtitle = element_text(hjust = 0.5, size=10)) + 
    ggtitle(class_data$v_label)  # Add cell_type as the title
  
  bar_chart = ggplot(pie_data, aes(x = Marker, y = percentage, fill = Marker)) +
    geom_bar(stat = "identity", width = 1, color = "black", size = 0.15, alpha=0.8) +
    scale_fill_manual(values = mycolor) +
    theme(legend.position="none", legend.text=element_text(size=8), legend.title = element_blank(),
          #legend.box.margin=margin(-10,-10,-10,-10),
          plot.title = element_text(hjust = 0.5, size=8),
          plot.subtitle = element_text(hjust = 0.5, size=8),
          axis.text.x = element_text(size = 7, angle=90), axis.text.y = element_text(size = 7),
          axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
          panel.background = element_rect(fill='white', colour='black'),
          panel.grid.major = element_line(color = "grey", linewidth=0.1, linetype=2)) +
    labs(title = class_data$v_label, x = "Marker", y = "% of + marker expression")
  print(bar_chart)
}




#' Visualize Pie Charts on a Graph Layout
#'
#' Generates pie charts for each cell type and overlays them on a network graph. 
#' Each pie chart represents the frequency of features above their median values 
#' in the sample data. Individual pie charts are saved as SVG files.
#'
#' @param data_sample A data frame of sample data, with the last column representing cell type labels.
#' @param pie_table A data frame containing node positions and sizes for the pie charts.
#' @param v A data frame of vertex information for the network graph.
#' @param e A data frame of edges for the network graph, including coordinates and weights.
#'
#' @return None. The function prints the combined plot and saves individual pie charts as SVG files.
visualizePieChart <- function(data_sample, pie_table, v, e) {
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  # Compute the threshold for each feature (median value)
  thresholds <- apply(data_sample[, -ncol(data_sample)], 2, median)
  
  # Convert data_sample to data.table for easier manipulation
  data_sample <- as.data.table(data_sample)
  
  # Compute the frequency of each feature above the threshold for each cell_type
  for (feature in names(thresholds)) {
    threshold <- thresholds[[feature]]
    freq_table <- data_sample[, .(freq = sum(get(feature) >= threshold)), by = cell_type]
    setnames(freq_table, "freq", feature)
    pie_table <- merge(pie_table, freq_table, by.x = "v_label", by.y = "cell_type", all.x = TRUE)
  }
  
  # Replace NA values with 0 (if any)
  pie_table[is.na(pie_table)] <- 0
  
  # Normalize the size values to a suitable range
  min_size <- min(pie_table$size)
  max_size <- max(pie_table$size)
  normalized_size <- (pie_table$size - min_size) / (max_size - min_size) + 0.5  # Adding 0.5 to ensure a minimum size
  
  # Create the base plot with edges and vertices
  g_plot <- ggplot() +
    geom_curve(data = e, aes(x = x.from, y = y.from, xend = x.to, yend = y.to, color = weight),
               linewidth = 1, curvature = 0.2, ncp = 100) +
    scale_color_gradient(low = "turquoise", high = "salmon", name = "Edge Weight",
                         guide = guide_colorbar(direction = "horizontal",
                                                title.position = "left",
                                                title.theme = element_text(size = 8),
                                                title.vjust = 0.75)) +
    coord_fixed() +
    theme_void() +
    theme(legend.direction = "horizontal", legend.position = "top",
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 10),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          plot.margin = unit(c(0, 0, 0, 0), "mm"))
  
  # Add pie charts to the plot
  for (i in 1:nrow(pie_table)) {
    class_data <- pie_table[i, ]
    pie_grob <- create_pie_grob(class_data, mycolor)
    size_factor <- normalized_size[i] 
    
    g_plot <- g_plot + annotation_custom(pie_grob, 
                                         xmin = class_data$x - 0.5*size_factor, xmax = class_data$x + 0.5*size_factor, 
                                         ymin = class_data$y - 0.5*size_factor + 0.2, ymax = class_data$y + 0.5*size_factor + 0.2)
    
    svglite(paste("../Output/", "PieGraph_Indiv_", class_data$v_label, ".svg", sep=""), width = 2.5, height = 2.5)
    plot(create_individual_pie_grob(class_data, mycolor))
    dev.off()
  }
  
  # Create a dummy plot to extract the legend
  dummy_plot <- ggplot(data.frame(Marker = names(thresholds), percentage = rep(1, length(thresholds))), 
                       aes(x = "", y = percentage, fill = Marker)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = mycolor) +
    theme_void() +
    theme(legend.position = "right", legend.text = element_text(size = 7), legend.title = element_text(size = 8))
  legend <- get_legend(dummy_plot)
  
  # Combine the main plot and the legend
  combined_plot <- plot_grid(g_plot, legend, ncol = 2, rel_widths = c(9, 1))
  
  # Display the combined plot
  print(combined_plot)
  
}




#' Visualize Marker Statistics by Cell Type
#'
#' Generates boxplots of marker values for each cell type in the dataset and 
#' saves them as SVG files.
#'
#' @param data A data frame containing marker values (columns 1 to n-1) and a cell type column (last column).
#' @param data_name A string used to name the output SVG files.
#'
#' @return None. The function prints plots and saves them as SVG files in the output folder.
visualizeMarkerStat <- function(data, data_name){
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  features <- colnames(data)[1:(ncol(data) - 1)]
  cell_type <- colnames(data)[ncol(data)]
  
  # Loop through each cell type and create a plot
  unique_cell_types <- unique(data[[cell_type]])
  
  for (ct in unique_cell_types) {
    subset_data <- data[data$cell_type==ct, ]
    meltData <- reshape::melt(subset_data)
    
    svglite(paste("../Output/marker stat/", data_name, "_", ct, ".svg", sep=""), width = 2.5, height = 2.5)
    p1 <- ggplot(meltData, aes(variable, value, fill=variable, alpha=0.3)) + 
      geom_boxplot(width=0.5, outlier.size = 0.08, size = 0.1,
                   staplewidth = 0.5, outliers = FALSE) + 
      theme(legend.position="none", legend.text=element_text(size=8), legend.title = element_blank(),
            #legend.box.margin=margin(-10,-10,-10,-10),
            plot.title = element_text(hjust = 0.5, size=8),
            plot.subtitle = element_text(hjust = 0.5, size=8),
            axis.text.x = element_text(size = 7, angle=90), axis.text.y = element_text(size = 7),
            axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
            panel.background = element_rect(fill='white', colour='black'),
            panel.grid.major = element_line(color = "grey", linewidth=0.1, linetype=2)) +
      scale_fill_manual(values = mycolor) +
      labs(title = ct, x = "Marker", y = "Marker value")
    print(p1)
    dev.off()
    
    # avg_subset_data = apply(subset_data[, 1:length(features)], 2, mean)
    # avg_subset_data = data.frame(variable=names(avg_subset_data), value=avg_subset_data)
    # 
    # p2 <- ggplot(avg_subset_data, aes(x=variable, y=value, fill=variable, alpha=0.3)) + 
    #   geom_bar(stat = "identity") +
    #   theme(legend.position="none", legend.text=element_text(size=8), legend.title = element_blank(),
    #         #legend.box.margin=margin(-10,-10,-10,-10),
    #         plot.title = element_text(hjust = 0.5, size=12),
    #         plot.subtitle = element_text(hjust = 0.5, size=10),
    #         axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
    #         axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
    #         panel.background = element_rect(fill='white', colour='black'),
    #         panel.grid.major = element_line(color = "grey", linewidth=0.25, linetype=2)) +
    #   scale_fill_manual(values = mycolor) +
    #   labs(title = paste("Marker Average for Cell Type:", ct), x = "Markers", y = "Marker Values")
    # print(p2)
  }
}




#' Compare Marker Statistics Across Multiple Time Points
#'
#' Generates boxplots of marker values for each cell type across three datasets 
#' representing different time points or conditions, and saves them as SVG files.
#'
#' @param data1 A data frame of marker values for the first condition (e.g., Baseline).
#' @param data2 A data frame of marker values for the second condition (e.g., After 1 vac).
#' @param data3 A data frame of marker values for the third condition (e.g., After 3 vac).
#'
#' @return None. The function prints plots and saves them as SVG files in the output folder.
compareMarkerStat <- function(data1, data2, data3){
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  features <- colnames(data1)[1:(ncol(data1) - 1)]
  cell_type <- colnames(data1)[ncol(data1)]
  
  # Loop through each cell type and create a plot
  unique_cell_types <- unique(data1[[cell_type]])
  
  for (ct in unique_cell_types) {
    subset_data1 <- data1[data1$cell_type==ct, ]
    avg_subset_data1 = apply(subset_data1[, 1:length(features)], 2, mean)
    avg_subset_data1 = data.frame(variable=names(avg_subset_data1), value=avg_subset_data1, time="Baseline")
    subset_data1$time = "Baseline"
    
    subset_data2 <- data2[data2$cell_type==ct, ]
    avg_subset_data2 = apply(subset_data2[, 1:length(features)], 2, mean)
    avg_subset_data2 = data.frame(variable=names(avg_subset_data2), value=avg_subset_data2, time="After 1 vac")
    subset_data2$time = "After 1 vac"
    
    subset_data3 <- data3[data3$cell_type==ct, ]
    avg_subset_data3 = apply(subset_data3[, 1:length(features)], 2, mean)
    avg_subset_data3 = data.frame(variable=names(avg_subset_data3), value=avg_subset_data3, time="After 3 vac")
    subset_data3$time = "After 3 vac"
    
    avg_subset_data = rbind(avg_subset_data1, avg_subset_data2, avg_subset_data3)
    subset_data = rbind(subset_data1, subset_data2, subset_data3)
    
    # Convert the data from wide to long format
    long_data <- subset_data %>%
      pivot_longer(cols = colnames(subset_data)[1:length(features)], names_to = "Marker", values_to = "Value")
    # Specify the order of the time variable
    long_data$time <- factor(long_data$time, 
                             levels = c("Baseline", "After 1 vac", "After 3 vac"))  # Adjust the levels as per your data
    
    markers <- colnames(subset_data)[1:length(features)]
    # Generate shades for each marker
    marker_colors <- lapply(mycolor[1:length(markers)], generate_shades)
    
    # Loop through each marker and create a plot
    for (i in seq_along(markers)) {
      marker <- markers[i]
      
      svglite(paste("../Output/compare marker stat/", "Compare_", ct, "_", marker, ".svg", sep=""), 
              width = 2.5, height = 2.5)
      
      plot <- ggplot(long_data[long_data$Marker == marker, ]) +
        geom_boxplot(aes(x = time, y = Value, fill = time),
                     color = "black", width = 0.5, outlier.size = 0.08, size = 0.1,
                     staplewidth = 0.5, outliers = FALSE) +
        facet_wrap(~ Marker) +
        theme(legend.position = "none", legend.text = element_text(size = 7), legend.title = element_blank(),
              legend.key=element_blank(),
              plot.title = element_text(hjust = 0.5, size = 8),
              plot.subtitle = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(size = 7, angle=0), axis.text.y = element_text(size = 7),
              axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = 2)) +
        scale_fill_manual(values = marker_colors[[i]]) +
        labs(title = ct, x = "Time", y = "Value")       
      print(plot)
      dev.off()
    }
  }
}




#' Export Ki67 Expression Percentages to Excel
#'
#' Extracts Ki67 expression percentages from precomputed summary tables for three datasets
#' and writes them to separate Excel files.
#'
#' @param tables A list of three summary data frames (table1, table2, table3) containing marker counts per cell type.
#' @param file_name_1 File name associated with the first dataset.
#' @param file_name_2 File name associated with the second dataset.
#' @param file_name_3 File name associated with the third dataset.
#'
#' @details 
#' The function computes Ki67-positive cell percentages by dividing the Ki67 counts by the total
#' cell count (`size`) for each cell type. It then writes the resulting percentages to Excel files
#' in the "../Output/compare marker expression/" directory. The file names are prefixed with "Ki67_".
#'
#' @return None (Excel files are saved to disk).
expressionData <- function(tables, file_name_1, file_name_2, file_name_3){
  table1 = tables[[1]]
  table2 = tables[[2]]
  table3 = tables[[3]]
  
  column_names = table1$v_label
  ki67_ind = which(tolower(colnames(table1))=="ki67")
  
  tab1 = data.frame(matrix(nrow=0, ncol=length(column_names)))
  tab1 = rbind(tab1, table1[, ki67_ind]/table1$size*100)
  colnames(tab1) = column_names
  
  tab2 = data.frame(matrix(nrow=0, ncol=length(column_names)))
  tab2 = rbind(tab2, table2[, ki67_ind]/table2$size*100)
  colnames(tab2) = column_names
  
  tab3 = data.frame(matrix(nrow=0, ncol=length(column_names)))
  tab3 = rbind(tab3, table3[, ki67_ind]/table3$size*100)
  colnames(tab3) = column_names
  
  write.xlsx(tab1, paste("../Output/compare marker expression/", "Ki67_", strsplit(file_name_1, "\\.")[[1]][1], ".xlsx", sep=""))
  write.xlsx(tab2, paste("../Output/compare marker expression/", "Ki67_", strsplit(file_name_2, "\\.")[[1]][1], ".xlsx", sep=""))
  write.xlsx(tab3, paste("../Output/compare marker expression/", "Ki67_", strsplit(file_name_3, "\\.")[[1]][1], ".xlsx", sep=""))
  
}




#' Compare Marker Expression Across Three Conditions
#'
#' Computes the frequency of cells expressing each marker above feature-specific thresholds 
#' for three datasets (e.g., baseline and post-vaccination samples) and visualizes the results 
#' using grouped bar plots.
#'
#' @param data1 Data frame for first condition (columns: markers + cell_type).
#' @param data2 Data frame for second condition (columns: markers + cell_type).
#' @param data3 Data frame for third condition (columns: markers + cell_type).
#' @param table1 Predefined summary table for data1 (must include "v_label" column for cell types).
#' @param table2 Predefined summary table for data2.
#' @param table3 Predefined summary table for data3.
#' @param file1 Filename string corresponding to data1 (used to extract threshold info for Ki67).
#' @param file2 Filename string corresponding to data2.
#' @param file3 Filename string corresponding to data3.
#'
#' @return A list containing updated summary tables: `table1`, `table2`, `table3`.
#'         Each table includes the computed frequencies of marker-positive cells per cell type.
#'         Additionally, grouped bar plots are generated for each cell type comparing the three conditions.
compareMarkerEx <- function(data1, data2, data3, table1, table2, table3, file1, file2, file3){
  mycolor = unique( c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), 
                      brewer.pal(8, "Accent"), brewer.pal(8, "Set2")) )
  
  degated_thresholds = read.xlsx("../Data/Raw Data/BB515_A_Threshold.xlsx", sheet = "Sheet 1")
  ind1 = paste0(strsplit(file1, "_CellType")[[1]][1], ".fcs")
  ind2 = paste0(strsplit(file2, "_CellType")[[1]][1], ".fcs")
  ind3 = paste0(strsplit(file3, "_CellType")[[1]][1], ".fcs")
  
  th1 = degated_thresholds[degated_thresholds$file_name==ind1, 2]
  th2 = degated_thresholds[degated_thresholds$file_name==ind2, 2]
  th3 = degated_thresholds[degated_thresholds$file_name==ind3, 2]
  cat(th1, th2, th3, "\n")
  
  th = median(th1, th2, th3)
  
  # Compute the threshold for each feature (median value)
  thresholds_1 <- apply(data1[, -ncol(data1)], 2, median)
  #### to match the results of the paper we only look into ki67 expression
  #### the expert opinion on ki67 threshold is 1.7
  #### from deGate function and averaging over 3 panels it is 1.267983
  #### so we use that, for the rest of the markers we use median value
  ki67_ind = which(tolower(colnames(data1))=="ki67")
  thresholds_1[ki67_ind] = th
  
  # Convert data_sample to data.table for easier manipulation
  data1 <- as.data.table(data1)
  # Compute the frequency of each feature above the threshold for each cell_type
  for (feature_1 in names(thresholds_1)) {
    threshold_1 <- thresholds_1[[feature_1]]
    
    freq_table_1 <- data1[, .(freq = sum(get(feature_1) >= threshold_1)), by = cell_type]
    
    setnames(freq_table_1, "freq", feature_1)
    table1 <- merge(table1, freq_table_1, by.x = "v_label", by.y = "cell_type", all.x = TRUE)
  }
  
  # Compute the threshold for each feature (median value)
  thresholds_2 <- apply(data2[, -ncol(data2)], 2, median)
  thresholds_2[ki67_ind] = th
  # Convert data_sample to data.table for easier manipulation
  data2 <- as.data.table(data2)
  # Compute the frequency of each feature above the threshold for each cell_type
  for (feature_2 in names(thresholds_2)) {
    threshold_2 <- thresholds_2[[feature_2]]
    freq_table_2 <- data2[, .(freq = sum(get(feature_2) >= threshold_2)), by = cell_type]
    setnames(freq_table_2, "freq", feature_2)
    table2 <- merge(table2, freq_table_2, by.x = "v_label", by.y = "cell_type", all.x = TRUE)
  }
  
  # Compute the threshold for each feature (median value)
  thresholds_3 <- apply(data3[, -ncol(data3)], 2, median)
  thresholds_3[ki67_ind] = th
  # Convert data_sample to data.table for easier manipulation
  data3 <- as.data.table(data3)
  # Compute the frequency of each feature above the threshold for each cell_type
  for (feature_3 in names(thresholds_3)) {
    threshold_3 <- thresholds_3[[feature_3]]
    freq_table_3 <- data3[, .(freq = sum(get(feature_3) >= threshold_3)), by = cell_type]
    setnames(freq_table_3, "freq", feature_3)
    table3 <- merge(table3, freq_table_3, by.x = "v_label", by.y = "cell_type", all.x = TRUE)
  }
  
  for (i in 1:nrow(table1)) {
    print(i)
    class_data_1 <- table1[i, ]
    class_data_2 <- table2[i, ]
    class_data_3 <- table3[i, ]
    
    create_grouped_barplot(class_data_1, class_data_2, class_data_3, mycolor)
  }
  
  return(list(table1, table2, table3))
}




#' Compare Ki67 Expression Across Timepoints
#'
#' Creates boxplots comparing Ki67 expression percentages for each cell type
#' across three timepoints: Baseline, After 1 vaccination, and After 3 vaccinations.
#'
#' @param baseline A data frame containing Ki67 expression percentages at Baseline.
#' @param vac1 A data frame containing Ki67 expression percentages After 1 vaccination.
#' @param vac3 A data frame containing Ki67 expression percentages After 3 vaccinations.
#' @param mycolor A vector of colors used to generate shades for each cell type.
#'
#' @details 
#' The function combines the three datasets into a single long-format data frame, assigns
#' a "Time" factor variable, and then loops through each cell type to create individual
#' boxplots showing the distribution of Ki67 expression across the three timepoints.
#' Shades for each cell type are generated using the `generate_shades` function.
#'
#' @return None (plots are displayed in the current R graphics device).
compareMarkerEx_2 <- function(baseline, vac1, vac3, mycolor) {
  all_exp = rbind(baseline, vac1, vac3)
  all_exp$Time = c(rep("Baseline", times = nrow(baseline)) , 
                   rep("After 1 vac", times = nrow(vac1)), 
                   rep("After 3 vac", times = nrow(vac3)))
  
  # Specify the order of the time variable
  all_exp$Time <- factor(all_exp$Time, 
                         levels = c("Baseline", "After 1 vac", "After 3 vac"))  # Adjust the levels as per your data
  
  cell_types = colnames(baseline)
  # Generate shades for each cell type
  cell_colors <- lapply(mycolor[1:length(cell_types)], generate_shades)
  
  for(i in c(1:ncol(baseline))){
    data = all_exp[, c(i, ncol(all_exp))]
    cell = colnames(data)[1]
    colnames(data) = c("Ki67_Expression", "Time")
    
    g <- ggplot(data) +
      geom_boxplot(aes(x = Time, y = Ki67_Expression, fill = Time),
                   color = "black", width = 0.5, outlier.size = 0.08, size = 0.4,
                   staplewidth = 0.5, outliers = F) +
      theme(legend.position = "right", legend.text = element_text(size = 7), legend.title = element_blank(),
            legend.key=element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9),
            plot.subtitle = element_text(hjust = 0.5, size = 8),
            axis.text.x = element_text(size = 7, angle=0), axis.text.y = element_text(size = 7),
            axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = 2)) +
      scale_fill_manual(values = cell_colors[[i]]) +
      labs(title = cell, x = "Time", y = "Percentage of Ki67+ expression")
    print(g)
  }
}




#' Create Grouped Barplots for Marker Expression Comparison
#'
#' Generates grouped barplots comparing the percentage of marker-positive cells across
#' three conditions (e.g., baseline, after 1 vaccination, after 3 vaccinations) 
#' for each marker in a given cell type.
#'
#' @param class_data_1 Summary data frame for the first condition (columns: v_label, size, and marker counts).
#' @param class_data_2 Summary data frame for the second condition (same format as class_data_1).
#' @param class_data_3 Summary data frame for the third condition (same format as class_data_1).
#' @param mycolor A vector of base colors to use for generating marker-specific shades.
#'
#' @details 
#' The function calculates the percentage of positive cells for each marker by dividing 
#' the counts by the total cell count (`size`). It then creates separate barplots per marker,
#' facetted by marker name, with the three time points as grouped bars. Plots are saved as SVG files.
#'
#' @return None (plots are saved to "../Output/compare marker expression/").
create_grouped_barplot <- function(class_data_1, class_data_2, class_data_3, mycolor){
  # Combine the data from the three datasets
  combined_data <- data.frame(
    Marker = rep(names(class_data_1)[-c(1:4)], 3),
    percentage = c(
      round(as.numeric(class_data_1[-c(1:4)]) / class_data_1$size, 2) * 100,
      round(as.numeric(class_data_2[-c(1:4)]) / class_data_2$size, 2) * 100,
      round(as.numeric(class_data_3[-c(1:4)]) / class_data_3$size, 2) * 100
    ),
    Time = rep(c("Baseline", "After 1 vac", "After 3 vac"), each = length(names(class_data_1)[-c(1:4)]))
  )
  
  # Specify the order of the time variable
  combined_data$Time <- factor(combined_data$Time, 
                               levels = c("Baseline", "After 1 vac", "After 3 vac"))  # Adjust the levels as per your data
  
  markers <- names(class_data_1)[-c(1:4)]
  # Generate shades for each marker
  marker_colors <- lapply(mycolor[1:length(markers)], generate_shades)
  
  for (i in seq_along(markers)) {
    marker <- markers[i]
    d = combined_data[combined_data$Marker == marker, ]
    
    if(all(d$percentage==0)){
      print("All zero values, no barplot to show.")
    }else{
      svglite(paste("../Output/compare marker expression/", "Compare_", class_data_1$v_label, "_", marker, ".svg", sep=""), 
              width = 2.5, height = 2.5)
      
      # Create the grouped barplot
      bar_chart <- ggplot(d) +
        geom_bar(aes(x = Time, y = percentage, fill = Time),
                 stat = "identity", position = "dodge", width = 0.7, color = "black", size = 0.15, alpha = 0.8) +
        facet_wrap(~ Marker) +
        theme(legend.position = "none", legend.text = element_text(size = 7), legend.title = element_blank(),
              legend.key=element_blank(),
              plot.title = element_text(hjust = 0.5, size = 8),
              plot.subtitle = element_text(hjust = 0.5, size = 8),
              axis.text.x = element_text(size = 7, angle=0), axis.text.y = element_text(size = 7),
              axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
              panel.background = element_rect(fill = 'white', colour = 'black'),
              panel.grid.major = element_line(color = "grey", linewidth = 0.1, linetype = 2)) +
        scale_fill_manual(values = marker_colors[[i]]) +
        labs(title = class_data_1$v_label, x = "Time", y = "% of + marker expression")
      plot(bar_chart)
      dev.off()
    }
  }
}




#' Compute Marker Thresholds for Flow Cytometry Data
#'
#' Processes FCS files from Co-stimulation, Co-inhibition, and Cytokine panels to 
#' compute thresholds for each marker using the deGate method. The function reads 
#' pre-processed, normalized flow cytometry data, calculates per-channel thresholds, 
#' saves the threshold data to Excel files, and generates density plots for each marker.
#'
#' @return A data frame containing the median thresholds for selected markers across 
#'         Co-stimulation, Co-inhibition, and Cytokine panels.
computeMarkerThreshold <- function(){
  ####co-stimulation
  costim_threshold_data = data.frame(matrix(nrow=0, ncol=24))
  
  folder_path = "../Data/Raw Data/Co-stimulation/Normalized Data"
  folder_name = "TP3"
  file_list = list.files(path = paste(folder_path, "/", folder_name, sep = ""), pattern = ".fcs", full.names = TRUE)
  
  for(file_path in file_list){
    file_name = strsplit(file_path, "/")[[1]][length(strsplit(file_path, "/")[[1]])]
    print(file_name)
    
    #### load the data; the data is already pre-processed (PeacoQC, compensated, transformed, normalized)
    set.seed(42)
    ff = suppressWarnings(flowCore::read.FCS(file_path, truncate_max_range=FALSE))  # applies a linearize transformation
    
    chan_th_list = c()
    chan_th_list = c(chan_th_list, file_name)
    for(chan in unname(ff@parameters@data$name)){
      chan_th = deGate(ff,channel=chan, verbose = FALSE, upper = TRUE, 
                       twin.factor = 0.98, after.peak = 1, alpha = 0.15, slope.w = 4, tinypeak.removal = 0.2)
      chan_th_list = c(chan_th_list, chan_th)
    }
    costim_threshold_data = rbind(costim_threshold_data, chan_th_list)
  }
  colnames(costim_threshold_data) = c("file_name", "FSC-A::FSC-A", "FSC-H::FSC-H", "FSC-W::FSC-W", 
                                      "SSC-A::SSC-A", "SSC-H::SSC-H", "SSC-W::SSC-W", 
                                      "APC-A::PD-1", "APC-Cy7-A::CD3", "Alexa 700-A::CD8", 
                                      "BB515-A::KI67", "BB700-A::CD137", "BV421-A::CCR7", "BV480-A::LD", 
                                      "BV605-A::CD56", "BV650-A::ICOS", "BV711-A::HLA-DR", "BV786-A::CD4", 
                                      "Original_ID::ORIGINAL_ID", "PE-A::FOXP3", "PE-CF594-A::CD45RA", 
                                      "PE-Cy7-A::CD28", "Time::TIME", "sample_id::sample_id")
  
  costim_threshold_data[-1] <- lapply(costim_threshold_data[-1], as.numeric)
  write.xlsx(costim_threshold_data,
             "../Data/Raw Data/costim_threshold_data.xlsx")
  
  
  ####co-inhibition
  coinhib_threshold_data = data.frame(matrix(nrow=0, ncol=24))
  
  folder_path = "../Data/Raw Data/Co-inhibition/Normalized Data"
  folder_name = "TP3"
  file_list = list.files(path = paste(folder_path, "/", folder_name, sep = ""), pattern = ".fcs", full.names = TRUE)
  
  for(file_path in file_list){
    file_name = strsplit(file_path, "/")[[1]][length(strsplit(file_path, "/")[[1]])]
    print(file_name)
    
    #### load the data; the data is already pre-processed (PeacoQC, compensated, transformed, normalized)
    set.seed(42)
    ff = suppressWarnings(flowCore::read.FCS(file_path, truncate_max_range=FALSE))  # applies a linearize transformation
    
    chan_th_list = c()
    chan_th_list = c(chan_th_list, file_name)
    for(chan in unname(ff@parameters@data$name)){
      chan_th = deGate(ff,channel=chan, verbose = FALSE, upper = TRUE, 
                       twin.factor = 0.98, after.peak = 1, alpha = 0.15, slope.w = 4, tinypeak.removal = 0.2)
      chan_th_list = c(chan_th_list, chan_th)
    }
    coinhib_threshold_data = rbind(coinhib_threshold_data, chan_th_list)
  }
  colnames(coinhib_threshold_data) = c("file_name", "FSC-A::FSC-A", "FSC-H::FSC-H", "FSC-W::FSC-W", 
                                       "SSC-A::SSC-A", "SSC-H::SSC-H", "SSC-W::SSC-W", "APC-A::PD1", 
                                       "APC-Cy7-A::CD3", "Alexa 700-A::CD8", "BB515-A::KI67", 
                                       "BB700-A::CTLA-4", "BV421-A::CCR7", "BV480-A::LIVE DEAD", 
                                       "BV605-A::CD56", "BV650-A::TIM3", "BV711-A::CD39", "BV786-A::CD4", 
                                       "Original_ID::ORIGINAL_ID", "PE-A::FOXP3", "PE-CF594-A::CD45RA", 
                                       "PE-Cy7-A::LAG3", "Time::TIME", "sample_id::sample_id")
  
  
  coinhib_threshold_data[-1] <- lapply(coinhib_threshold_data[-1], as.numeric)
  write.xlsx(coinhib_threshold_data,
             "../Data/Raw Data/coinhib_threshold_data.xlsx")
  
  
  #### cytokine 
  cytokine_threshold_data = data.frame(matrix(nrow=0, ncol=24))
  
  folder_path = "../Data/Raw Data/Cytokine/Normalized Data"
  folder_name = "TP3"
  file_list = list.files(path = paste(folder_path, "/", folder_name, sep = ""), pattern = ".fcs", full.names = TRUE)
  
  for(file_path in file_list){
    file_name = strsplit(file_path, "/")[[1]][length(strsplit(file_path, "/")[[1]])]
    print(file_name)
    
    #### load the data; the data is already pre-processed (PeacoQC, compensated, transformed, normalized)
    set.seed(42)
    ff = suppressWarnings(flowCore::read.FCS(file_path, truncate_max_range=FALSE))  # applies a linearize transformation
    
    chan_th_list = c()
    chan_th_list = c(chan_th_list, file_name)
    for(chan in unname(ff@parameters@data$name)){
      chan_th = deGate(ff,channel=chan, verbose = FALSE, upper = TRUE, 
                       twin.factor = 0.98, after.peak = 1, alpha = 0.15, slope.w = 4, tinypeak.removal = 0.2)
      chan_th_list = c(chan_th_list, chan_th)
    }
    cytokine_threshold_data = rbind(cytokine_threshold_data, chan_th_list)
  }
  colnames(cytokine_threshold_data) = c("file_name", "FSC-A::FSC-A", "FSC-H::FSC-H", "FSC-W::FSC-W", 
                                        "SSC-A::SSC-A", "SSC-H::SSC-H", "SSC-W::SSC-W", "APC-A::PD1", 
                                        "APC-Cy7-A::CD3", "Alexa 700-A::CD8", "BB515-A::GRANZYME B", 
                                        "BB700-A::TNFA", "BV421-A::CCR7", "BV480-A::LIVE DEAD", 
                                        "BV605-A::CD56", "BV650-A::IL2", "BV711-A::IFNY", "BV786-A::CD4", 
                                        "Original_ID::ORIGINAL_ID", "PE-A::FOXP3", "PE-CF594-A::CD45RA", 
                                        "PE-Cy7-A::IL-10", "Time::TIME", "sample_id::sample_id")
  
  cytokine_threshold_data[-1] <- lapply(cytokine_threshold_data[-1], as.numeric)
  write.xlsx(cytokine_threshold_data,
             "../Data/Raw Data/cytokine_threshold_data.xlsx")
  
  #### channel density and stat
  for (i in c(2: ncol(costim_threshold_data))) {
    plot(density(costim_threshold_data[, i]), main=names(costim_threshold_data)[i])
  }
  for (i in c(2: ncol(coinhib_threshold_data))) {
    plot(density(coinhib_threshold_data[, i]), main=names(coinhib_threshold_data)[i])
  }
  for (i in c(2: ncol(cytokine_threshold_data))) {
    plot(density(cytokine_threshold_data[, i]), main=names(cytokine_threshold_data)[i])
  }
  
  costim_threshold_mean = as.data.frame(colMeans(costim_threshold_data[, c(-1)]))
  coinhib_threshold_mean = as.data.frame(colMeans(coinhib_threshold_data[, c(-1)]))
  cytokine_threshold_mean = as.data.frame(colMeans(cytokine_threshold_data[, c(-1)]))
  
  costim_threshold_median = as.data.frame(apply(costim_threshold_data[-1], 2, median, na.rm = TRUE))
  coinhib_threshold_median = as.data.frame(apply(coinhib_threshold_data[-1], 2, median, na.rm = TRUE))
  cytokine_threshold_median = as.data.frame(apply(cytokine_threshold_data[-1], 2, median, na.rm = TRUE))
  
  all_median = data.frame(coinhib=coinhib_threshold_median[c(8, 9 , 12, 14, 17, 19, 20),], 
                          costim=costim_threshold_median[c(8, 9 , 12, 14, 17, 19, 20),], 
                          cytokine=cytokine_threshold_median[c(8, 9 , 12, 14, 17, 19, 20),])
  rownames(all_median)= rownames(coinhib_threshold_median)[c(8, 9 , 12, 14, 17, 19, 20)]
  
  all_median$med_med = apply(all_median, 1, median, na.rm = TRUE)
}
