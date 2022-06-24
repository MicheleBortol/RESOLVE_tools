rm(list = ls())
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)
library(Giotto)
library(magick)

### Data files to read
CELL_COORD_PATH <- "Resolve_Mm_BHA1.v2/Resolve_Mm_BHA1_DAPI.tiff_coordinates.csv"
CELL_COUNT_PATH <- "Resolve_Mm_BHA1.v2/Resolve_Mm_BHA1_DAPI.tiff_measurements.csv"
RAW_TRANSCRIPT_PATH <- "Resolve_Mm_BHA1.v2/Resolve_Mm_BHA1_raw_locations.v2.txt"
CELL_ROI_PATH <- "Resolve_Mm_BHA1.v2/Resolve_Mm_BHA1_ROI_via_DAPI.zip"
DAPI_IMAGE_PATH <- "Resolve_Mm_BHA1.v2/Resolve_Mm_BHA1_DAPI.tiff"

### Load and process the DAPI image
# Giotto images are very memory heavy, so we need to downscale
DAPI_scale_factor <- 0.33 # Between 0 and 1
DAPI_img <- image_read(DAPI_IMAGE_PATH)

# Just for visualization, equalize the image
DAPI_img <- image_equalize(DAPI_img)

# Dimensions before rescaling, used to fix cell coordinates
DAPI_width <- image_info(DAPI_img)$width
DAPI_height <- image_info(DAPI_img)$height

DAPI_img <- image_scale(DAPI_img, paste0("%", DAPI_scale_factor * 100))
DAPI_boundaries <- c("xmax_adj" = 0, "xmin_adj" = 0,
	"ymax_adj" = 0,	"ymin_adj" = 0)

#### Process the cell coordinates for the @meta.data and the @images slots
cell_coordinates <- fread(CELL_COORD_PATH)

cell_coordinates <- cell_coordinates[, .(cell_id = V1, x = X,
	y = DAPI_height - Y, x_tile = X, y_tile = DAPI_height - Y,
	size = Width * Height)]

cell_coordinates[, cell_id := paste0("Cell_", cell_id - 1)]

### Process the coordinates of the segmentation ROI polygons for the @image slot
# Can read directly from the zip file in imgeJ input
cell_roi <- read.ijzip(CELL_ROI_PATH, names = TRUE, list.files = FALSE, verbose = FALSE)
cell_roi <- lapply(cell_roi, `[[`, "coords")
cell_roi <- lapply(cell_roi, as.data.table)
cell_roi <- rbindlist(cell_roi, idcol = "cell")
cell_roi[, y := DAPI_height - y]

cell_roi <- as.data.frame(cell_roi)

### Process the measurements to build the count matrix
measurements <- fread(CELL_COUNT_PATH)

setnames(measurements, "V33823", "Cell_33821") # To Do: better parsing
measurements <- measurements[, c(cell_coordinates$cell_id, "Gene-ID"), with = F]

measurements[`Gene-ID` == "_FP", `Gene-ID` := "FP"]

markers <- measurements[, `Gene-ID`]
count_matrix <- as.matrix(measurements[, -"Gene-ID"])
row.names(count_matrix) <- markers
colnames(count_matrix) <- colnames(measurements)[colnames(measurements) != "Gene-ID"]

# Process the coordinates of the raw transcripts for the @image slot
# Only for visualizing the raw transcripts (can be skipped)
transcripts <- fread(RAW_TRANSCRIPT_PATH)
transcripts <- transcripts[, .(x = V1, y = DAPI_height - V2, gene = V4)]
transcripts[gene == "_FP", gene := "FP"]
transcripts <- as.data.frame(transcripts)

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
row.names(fov_coordinates) <- fov_coordinates$cell

# Adding the spatial information
resolve.obj@images <- list()
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

##### Quick test single cell analysis, not meaningful just to have something to plot
resolve.obj <- subset(resolve.obj, subset = nCount_Spatial > 10) # Skip cells with no transcripts 
resolve.obj <- SCTransform(resolve.obj, assay = "Spatial", verbose = FALSE)
resolve.obj <- FindVariableFeatures(resolve.obj, assay = "Spatial")
resolve.obj <- RunPCA(object = resolve.obj, npcs = 20, verbose = FALSE)
resolve.obj <- RunUMAP(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindNeighbors(object = resolve.obj, dims = 1:20, verbose = FALSE)
resolve.obj <- FindClusters(object = resolve.obj, verbose = FALSE,
  resolution = 0.1)

resolve.obj@active.assay <- "Spatial" # To show the real counts in the plots

###############
# Change the definitions of this function from the SeuratObject package
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
# From seurat to Giotto (only useful if we want to use Seurat before,
# if not we can just make a giotto object from scratch)
# See: https://github.com/RubD/Giotto/blob/suite/R/interoperability.R
# seuratToGiotto

# The SeuratToGiotto function supports only one Assay and image in the Seurat object
resolve.obj@assays <- resolve.obj@assays["Spatial"]
resolve.obj@images <- resolve.obj@images["cen"]
#
# Centroids are sufficient, giotto does not support ROIs
giotto_res <- seuratToGiotto(resolve.obj, subcellular_assay = "cen")

giotto_img <- createGiottoImage(giotto_res, spatial_locs = NULL,
	mg_object = DAPI_img, name = "image", xmax_adj = 0, xmin_adj = 0,
	ymax_adj = 0, ymin_adj = 0, scale_factor = DAPI_scale_factor,
	do_manual_adj = T)
# Even when do_manual_adj is TRUE Giotto will force some changes,
# Since the coordinates are already aligned we want to force them again to 0
giotto_img@boundaries <- DAPI_boundaries

giotto_res <- addGiottoImage(giotto_res, list(giotto_img))

############# Plotting test
thin_cluster_plot = spatPlot2D(giotto_res, show_image = T, image_name = "image",
	point_alpha = 0.1, point_size = 1, cell_color = "seurat_clusters",
	save_plot = T)

thick_cluster_plot = spatPlot2D(giotto_res, show_image = T, image_name = "image",
	cell_color = "seurat_clusters", save_plot = T)

############ Spatial analysis test
# Adapted From:
# https://rubd.github.io/Giotto_site/articles/tut12_giotto_cell_interaction.html
giotto_res <- normalizeGiotto(giotto_res)

giotto_res <- createSpatialNetwork(gobject = giotto_res, minimum_k = 2)

## identify genes with a spatial coherent expression profile
km_spatialgenes <- binSpect(giotto_res, bin_method = "kmeans")

cell_proximities <- cellProximityEnrichment(giotto_res,
	cluster_column = "seurat_clusters",
	spatial_network_name = "Delaunay_network",
	adjust_method = "fdr", number_of_simulations = 1000)

cellProximityNetwork(giotto_res, CPscore = cell_proximities,
	remove_self_edges = T, only_show_enrichment_edges = T)

spec_interaction = "1--6" # Random significant interaction

cellproxplot <- cellProximitySpatPlot2D(giotto_res,
	interaction_name = spec_interaction,
	show_network = F, cluster_column = "seurat_clusters",
	cell_color = "seurat_clusters",	cell_color_code = c("blue", "red"),
	point_size_select = 1)
cellproxplot <- addGiottoImageToSpatPlot(cellproxplot, giotto_img)
cowplot::save_plot(plot = cellproxplot, device = "pdf",
	filename = "cellproxplot.pdf", base_width = 10, base_height = 10)

