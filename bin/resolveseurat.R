rm(list = ls())
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)

# https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system

### Data files to read
CELL_COORD_PATH <- "Projects/Resolve_Seurat/data//Resolve_Mm_BHA1_DAPI.tiff_coordinates.csv"
CELL_COUNT_PATH <- "Projects/Resolve_Seurat/data/Resolve_Mm_BHA1_DAPI.tiff_measurements.csv"
RAW_TRANSCRIPT_PATH <- "Projects/Resolve_Seurat/data/Resolve_Mm_BHA1_raw_locations.v2.txt"
CELL_ROI_PATH <- "Projects/Resolve_Seurat/data/Resolve_Mm_BHA1_ROI_via_DAPI.zip"

#### Process the cell coordinates for the @meta.data and the @images slots  
cell_coordinates <- fread(CELL_COORD_PATH) 

cell_coordinates <- cell_coordinates[, .(cell_id = V1, x = Y, y = X,
  x_tile = Y, y_tile = X, size = Width * Height)]

# Rotate the coordinates by 90° so they match the images
# Might be better to do the swap directly in the plots as a ggplot2 annotation
#
DAPI_HEIGHT <- 32160 # 32160 = Height of the DAPI image
cell_coordinates[, x := DAPI_HEIGHT - x] 
cell_coordinates[, x_tile := DAPI_HEIGHT - x_tile]

cell_coordinates[, cell_id := paste0("Cell_", cell_id - 1)]

### Process the measurements to build the count matrix
measurements <- fread(CELL_COUNT_PATH)
setnames(measurements, "V33823", "Cell_33821") # To Do: better parsing
measurements <- measurements[, c(cell_coordinates$cell_id, "Gene-ID"), with = F]

markers <- measurements[, `Gene-ID`]
count_matrix <- as.matrix(measurements[, -"Gene-ID"])
row.names(count_matrix) <- markers
colnames(count_matrix) <- colnames(measurements)[colnames(measurements) != "Gene-ID"]  

### Process the coordinates of the raw transcripts for the @image slot
# Only for visualizing the raw transcripts (can be skipped)
transcripts <- fread(RAW_TRANSCRIPT_PATH)
transcripts <- transcripts[, .(x = DAPI_HEIGHT - V2, y = V1, gene = V4)]
transcripts <- as.data.frame(transcripts)

### Process the coordinates of the segmentation ROI polygons for the @image slot
# Can read directly from the zip file in imgeJ input
cell_roi <- read.ijzip(CELL_ROI_PATH, names = TRUE, list.files = FALSE, verbose = FALSE)
cell_roi <- lapply(cell_roi, `[[`, "coords")
cell_roi <- lapply(cell_roi, as.data.table)
cell_roi <- rbindlist(cell_roi, idcol = "cell")

cell_roi[, x2 := y] # reorienting the coordinates
cell_roi[, y := x]
cell_roi[, x := DAPI_HEIGHT - x2]
cell_roi[, x2 := NULL]

cell_roi <- as.data.frame(cell_roi)

### Create the Seurat Object
resolve.obj <- CreateSeuratObject(counts = count_matrix, assay="Spatial")

resolve.obj@meta.data["cell_id"] <- rownames(resolve.obj@meta.data)
resolve.obj@meta.data <- as.data.table(resolve.obj@meta.data)
# Assumes all the cells come from one image !!!!
resolve.obj@meta.data[, region := 1]
resolve.obj@meta.data[, z := 0]
resolve.obj@meta.data[, tile_num := 0]
resolve.obj@meta.data <- merge(resolve.obj@meta.data, cell_coordinates,
 by = "cell_id", sort = F)

resolve.obj@meta.data <- as.data.frame(resolve.obj@meta.data)
rownames(resolve.obj@meta.data) <- resolve.obj@meta.data$cell_id 

fov_coordinates <- as.data.frame(cell_coordinates[, .(x, y, cell = cell_id)])
rownames(fov_coordinates) <- fov_coordinates$cell

# Adding the spatial information for both centroids and segmentation ROIs.
# If we load the segmentation ROIs, the centroids are not really necessary. 
resolve.obj@images <- list()
resolve.obj@images["cen"] <- CreateFOV(
  coords = fov_coordinates,
  type = "centroids",
  nsides = 0L,
  radius = 1L,
  theta = 0L,
  molecules = transcripts, # Only for visualising the raw transcripts, can be skipped
  assay = "Spatial",
  key = NULL,
  name = NULL)

resolve.obj@images["seg"] <- CreateFOV(
  coords = cell_roi,
  type = "segmentation",
  nsides = 0L,
  radius = 1L,
  theta = 0L,
  molecules = transcripts, # Only for visualising the raw transcripts, can be skipped
  assay = "Spatial",
  key = NULL,
  name = NULL)

#### Quick test single cell analysis, not meaningful just to have something to plot
resolve.obj <- subset(resolve.obj, subset = nCount_Spatial > 10) # Skip cells with no transcripts 
resolve.obj <- SCTransform(resolve.obj, assay = "Spatial", verbose = FALSE)
resolve.obj <- FindVariableFeatures(resolve.obj, assay = "Spatial")
resolve.obj <- RunPCA(object = resolve.obj, npcs = 20, verbose = FALSE)
resolve.obj <- RunUMAP(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindNeighbors(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindClusters(object = resolve.obj, verbose = FALSE,
  resolution = 0.1)

resolve.obj@active.assay <- "Spatial" # To show the real counts in the plots


### Plotting

# No spatial data
nos_dp <- DimPlot(resolve.obj, label = TRUE, label.box = TRUE) + NoLegend() 
nos_fp <- FeaturePlot(resolve.obj, features = c("Slc17a7", "Bcl6", "Foxp2",
  "Sox10"), min.cutoff = "q10", max.cutoff = "q90") 

# Plus spatial data (cells as circles on their centroids)
cen_dp <- ImageDimPlot(resolve.obj, cols = "parade", fov = "cen")
cen_fp <- ImageFeaturePlot(resolve.obj, fov = "cen", features = c("Slc17a7",
  "Bcl6", "Foxp2", "Sox10"))

# Plus spatial data (cells as segmentation outlines)
seg_dp <- ImageDimPlot(resolve.obj, cols = "parade", fov = "seg")
seg_fp <- ImageFeaturePlot(resolve.obj, features = c("Slc17a7", "Bcl6",
  "Foxp2", "Sox10"), min.cutoff = "q10", max.cutoff = "q90", fov = "seg")

# Raw transcript visualization
# - size = 0 avoids displaying the cells
# - nmols = nrow(transcripts) displays all the transcripts, instead of a sample 
cen_mp <- ImageFeaturePlot(resolve.obj, fov = "cen", features = "Sox10", size = 0,
  molecules = c("Slc17a7", "Bcl6", "Foxp2", "Sox10"), mols.cols = c("#FF0000",
 "#00FF00", "#0000FF", "#FFFF00"), nmols = nrow(transcripts))

seg_mp <- ImageFeaturePlot(resolve.obj, fov = "seg", features = "Bcl6", 
  molecules = c("Foxp2"), mols.cols = c("#0000FF"), nmols = nrow(transcripts),
  min.cutoff = "q10", max.cutoff = "q90")


####################
# Find Spatially variable features
# The coordinates for the spatial analysis can come from either the centroids or
# the polygons from the cell ROIs. The  GetTissueCoordinates.Centroids or
# GetTissueCoordinates.Segmentation need to be modified. See below.
#
# Both the "moransi" and the "markvariogram" work but the "moransi" method
# requires the installation of the "ape" or the "Rfast2" package

###############
# Change the definitions of these 2 functions from the SeuratObject package
# This is not good practice, but it works as a quick workaround for the
# FindSpatiallyVariableFeatures function which require the returned object
# to have rownames. 
#
# see the fortify.Centroids and the FindSpatiallyVariableFeatures functions
GetTissueCoordinates.Centroids <- function(object, full = TRUE, ...){
  a <- as.data.frame(object@coords)
  a$cell <- object@cells
  rownames(a) <- object@cells
  return(a)
}

assignInNamespace("GetTissueCoordinates.Centroids",
                  GetTissueCoordinates.Centroids, ns = "SeuratObject")

# see the fortify.Segmentation function and the FindSpatiallyVariableFeatures functions
GetTissueCoordinates.Segmentation <- function(object, full = TRUE, ...){
  a <- rbindlist(lapply(object@polygons , function(o){
    as.data.frame(t(o@labpt))
  }), id = "cell")
  b <- data.frame(x = a$V1, y = a$V2, cell = a$cell) # Columns need to be in this order
  row.names(b) <- b$cell
  return(b)
}

assignInNamespace("GetTissueCoordinates.Segmentation",
                  GetTissueCoordinates.Segmentation, ns = "SeuratObject")
################

##### Remove this line
# Just to reduce the size of the dataset for testing locally
#resolve.obj <- subset(resolve.obj, subset = nCount_Spatial > 1000)
########
resolve.obj <- ScaleData(resolve.obj) # 

sfs_method <- "markvariogram" # moransi also runs 
resolve.obj <- FindSpatiallyVariableFeatures(resolve.obj, assay = "Spatial", image = "seg",
  selection.method = sfs_method)
top_features <- SpatiallyVariableFeatures(resolve.obj, decreasing = T,
  selection.method = sfs_method)

top_fp <- ImageFeaturePlot(resolve.obj, fov = "cen", features = "Vip", size = 0,
  molecules = top_features[1:3], nmols = nrow(transcripts)) 
