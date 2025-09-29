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

library(devtools)
library(flowCore)
library(FlowSOM)
library(CytoNorm)
library(ggpubr)

#### For T-cells we have 3 panels and 2 groups of patients
#### Panel 1: co-inhibition, Panel 2: co-stimulation, Panel 3: cytokine
#### Group 1: day 1, Group 2: day 2, and this is consistent for all 3 time points
#### now we organize the preprocessed data in batches

files_P1_G1 = c(list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition/Normalized_Cell_Type_Marker/TP1", pattern="Norm_coinhib_g1.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition/Normalized_Cell_Type_Marker/TP2", pattern="Norm_coinhib_g1.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition/Normalized_Cell_Type_Marker/TP3", pattern="Norm_coinhib_g1.*.fcs", full.names = TRUE))

files_P1_G2 = c(list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition/Normalized_Cell_Type_Marker/TP1", pattern="Norm_coinhib_g2.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition/Normalized_Cell_Type_Marker/TP2", pattern="Norm_coinhib_g2.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-inhibition/Normalized_Cell_Type_Marker/TP3", pattern="Norm_coinhib_g2.*.fcs", full.names = TRUE))


files_P2_G1 = c(list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation/Normalized_Cell_Type_Marker/TP1", pattern="Norm_costim_g1.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation/Normalized_Cell_Type_Marker/TP2", pattern="Norm_costim_g1.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation/Normalized_Cell_Type_Marker/TP3", pattern="Norm_costim_g1.*.fcs", full.names = TRUE))

files_P2_G2 = c(list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation/Normalized_Cell_Type_Marker/TP1", pattern="Norm_costim_g2.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation/Normalized_Cell_Type_Marker/TP2", pattern="Norm_costim_g2.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Co-stimulation/Normalized_Cell_Type_Marker/TP3", pattern="Norm_costim_g2.*.fcs", full.names = TRUE))


files_P3_G1 = c(list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Normalized_Cell_Type_Marker/TP1", pattern="Norm_cytokine_g1.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Normalized_Cell_Type_Marker/TP2", pattern="Norm_cytokine_g1.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Normalized_Cell_Type_Marker/TP3", pattern="Norm_cytokine_g1.*.fcs", full.names = TRUE))

files_P3_G2 = c(list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Normalized_Cell_Type_Marker/TP1", pattern="Norm_cytokine_g2.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Normalized_Cell_Type_Marker/TP2", pattern="Norm_cytokine_g2.*.fcs", full.names = TRUE),
                list.files("C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Cytokine/Normalized_Cell_Type_Marker/TP3", pattern="Norm_cytokine_g2.*.fcs", full.names = TRUE))


#### Read in 1 file per panel to get channel information
ff_P1 <- read.FCS(files_P1_G1[1])
ff_P2 <- read.FCS(files_P2_G2[1])
ff_P3 <- read.FCS(files_P3_G2[1])

##### Retrieve channel information
#### the cell type markers are common in the 3 panels, so extracting from the first panel is ok
cellTypeMarkers <- c("CD56", "CD3", "CD4", "FOXP3", "CD8", "CD45RA", "CCR7")
cellTypeChannels <- GetChannels(object = ff_P1, markers = cellTypeMarkers)

#### the cell state markers are different in the 3 panels, extracting all separately
P1_cellStateMarkers <- c("LAG3", "PD1", "TIM3", "CD39", "KI67", "CTLA-4")
P1_cellStateChannels <- GetChannels(object = ff_P1, markers = P1_cellStateMarkers)

P2_cellStateMarkers <- c("CD28", "CD137", "PD-1", "HLA-DR", "ICOS", "KI67")
P2_cellStateChannels <- GetChannels(object = ff_P2, markers = P2_cellStateMarkers)

P3_cellStateMarkers <- c("PD1", "GRANZYME B", "IL-10", "TNFA", "IL2", "IFNY")
P3_cellStateChannels <- GetChannels(object = ff_P3, markers = P3_cellStateMarkers)

##### Create aggregate per batch with 10,000 cells from each file
set.seed(42)
agg_P1_G1 <- AggregateFlowFrames(fileNames = files_P1_G1, cTotal = length(files_P1_G1)*10000)
agg_P1_G2 <- AggregateFlowFrames(fileNames = files_P1_G2, cTotal = length(files_P1_G2)*10000)

agg_P2_G1 <- AggregateFlowFrames(fileNames = files_P2_G1, cTotal = length(files_P2_G1)*10000)
agg_P2_G2 <- AggregateFlowFrames(fileNames = files_P2_G2, cTotal = length(files_P2_G2)*10000)

agg_P3_G1 <- AggregateFlowFrames(fileNames = files_P3_G1, cTotal = length(files_P3_G1)*10000)
agg_P3_G2 <- AggregateFlowFrames(fileNames = files_P3_G2, cTotal = length(files_P3_G2)*10000)

#### train and normalize the cell state markers ###############################################################
model_cellstate_P1_G1 <- CytoNorm.train(files = flowSet(agg_P1_G1),
                                        labels = c("P1"),
                                        channels = P1_cellStateChannels,
                                        transformList = NULL,
                                        seed = 1, verbose = TRUE, plot = FALSE,
                                        FlowSOM.params = list(nCells = 1e+06, xdim = 10, ydim = 10, nClus = 6, scale = FALSE))
# Normalize files
CytoNorm.normalize(model = model_cellstate_P1_G1,
                   files = c(files_P1_G1),
                   labels = c(rep("P1", length(files_P1_G1)) ),
                   transformList = NULL,
                   verbose = TRUE,
                   refix = "",
                   transformList.reverse = NULL,
                   outputDir = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Normalized_Tcell_cellState_P1_G1")

###############################################################################################################
model_cellstate_P1_G2 <- CytoNorm.train(files = flowSet(agg_P1_G2),
                                        labels = c("P1"),
                                        channels = P1_cellStateChannels,
                                        transformList = NULL,
                                        seed = 1, verbose = TRUE, plot = FALSE,
                                        FlowSOM.params = list(nCells = 1e+06, xdim = 10, ydim = 10, nClus = 6, scale = FALSE))
# Normalize files
CytoNorm.normalize(model = model_cellstate_P1_G2,
                   files = c(files_P1_G2),
                   labels = c(rep("P1", length(files_P1_G2)) ),
                   transformList = NULL,
                   verbose = TRUE,
                   refix = "",
                   transformList.reverse = NULL,
                   outputDir = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Normalized_Tcell_cellState_P1_G2")

###############################################################################################################
model_cellstate_P2_G1 <- CytoNorm.train(files = flowSet(agg_P2_G1),
                                        labels = c("P2"),
                                        channels = P2_cellStateChannels,
                                        transformList = NULL,
                                        seed = 1, verbose = TRUE, plot = FALSE,
                                        FlowSOM.params = list(nCells = 1e+06, xdim = 10, ydim = 10, nClus = 6, scale = FALSE))
# Normalize files
CytoNorm.normalize(model = model_cellstate_P2_G1,
                   files = c(files_P2_G1),
                   labels = c(rep("P2", length(files_P2_G1)) ),
                   transformList = NULL,
                   verbose = TRUE,
                   refix = "",
                   transformList.reverse = NULL,
                   outputDir = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Normalized_Tcell_cellState_P2_G1")

###############################################################################################################
model_cellstate_P2_G2 <- CytoNorm.train(files = flowSet(agg_P2_G2),
                                        labels = c("P2"),
                                        channels = P2_cellStateChannels,
                                        transformList = NULL,
                                        seed = 1, verbose = TRUE, plot = FALSE,
                                        FlowSOM.params = list(nCells = 1e+06, xdim = 10, ydim = 10, nClus = 6, scale = FALSE))
# Normalize files
CytoNorm.normalize(model = model_cellstate_P2_G2,
                   files = c(files_P2_G2),
                   labels = c(rep("P2", length(files_P2_G2)) ),
                   transformList = NULL,
                   verbose = TRUE,
                   refix = "",
                   transformList.reverse = NULL,
                   outputDir = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Normalized_Tcell_cellState_P2_G2")

###############################################################################################################
model_cellstate_P3_G1 <- CytoNorm.train(files = flowSet(agg_P3_G1),
                                        labels = c("P3"),
                                        channels = P3_cellStateChannels,
                                        transformList = NULL,
                                        seed = 1, verbose = TRUE, plot = FALSE,
                                        FlowSOM.params = list(nCells = 1e+06, xdim = 10, ydim = 10, nClus = 6, scale = FALSE))
# Normalize files
CytoNorm.normalize(model = model_cellstate_P3_G1,
                   files = c(files_P3_G1),
                   labels = c(rep("P3", length(files_P3_G1)) ),
                   transformList = NULL,
                   verbose = TRUE,
                   refix = "",
                   transformList.reverse = NULL,
                   outputDir = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Normalized_Tcell_cellState_P3_G1")

###############################################################################################################
model_cellstate_P3_G2 <- CytoNorm.train(files = flowSet(agg_P3_G2),
                                        labels = c("P3"),
                                        channels = P3_cellStateChannels,
                                        transformList = NULL,
                                        seed = 1, verbose = TRUE, plot = FALSE,
                                        FlowSOM.params = list(nCells = 1e+06, xdim = 10, ydim = 10, nClus = 6, scale = FALSE))
# Normalize files
CytoNorm.normalize(model = model_cellstate_P3_G2,
                   files = c(files_P3_G2),
                   labels = c(rep("P3", length(files_P3_G2)) ),
                   transformList = NULL,
                   verbose = TRUE,
                   refix = "",
                   transformList.reverse = NULL,
                   outputDir = "C:/Users/AbidaSanjanaShemonti/Documents/R Projects/FlowViz/Data/Raw Data/Normalized_Tcell_cellState_P3_G2")

  
