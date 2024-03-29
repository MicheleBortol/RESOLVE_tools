# ARCHIVED REPOSITORY:
## MOVED TO: https://gitlab.hzdr.de/michele.bortolomeazzi/RESOLVE_tools

# RESOLVE_tools
Tentative set of tools and scripts for analysing spatial transcriptomic data with the resolve platform

## Contents:
+ A) [Registration](#Registration)
+ B) [PreProcessing](#PreProcessing)
    + (B.1) [Nextflow pipeline](##Pipeline)
        + (B.1.1) [Parameters](###Parameters)
        + (B.1.2) [Input](###Input)
        + (B.1.3) [Output](###Output)
        + (B.1.4) [Example Run](###Example)
    + (B.2) [Scripts](##Scripts)
        + (B.2.1) [Segmentation](###Segmentation)
        + (B.2.2) [Expression assignment](###expression_assign)
+ C) [Analysis](#Analysis)
  + (C.1) [Seurat](##Seurat)
  + (C.2) [Giotto](##Giotto)

# A) Registration <a name="Registration"></a>
FiJi [macro](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/FIJI_REGISTRATION.ijm) for the registration of DAPI images from the RESOLVE and confocal microscopes.

**Usage**  
The script can be run interactively as a Fiji Macro. The script works as follows:
1) Ask the user for the pahts to open the `source` (Confocal) and `target` (Resolve) images.
2) Vertically mirrors the `source` to match the `target`.
3) Extract features with SIFT from both images and match them.
4) Tranform the `source` image using the extracted features to identify the transformation.
5) Extract only the `registered` image from the stack generated at point 4).
6) Change the `registered` image from 32 back to 16 bit (same bit depth of ´source´ and ´target´).
7) Create a RED (`registered`) and GREEN (`target`) image for checking the registration quality.
8) Save the `registered` image to a directory selected by the user as: `registered.tiff`.

**Warning:**
+ The script can require a lot of memory for large images.
+ If there are many empty areas in the target image, the registered image will be heavily distorted in their proximity.
+ The script assumes that the source image needs to be vertically mirrored.

**To Do**  
Convert to a script which can be run non-interactively. Hopefully with open-cv instead of FiJi.

# B) PreProcessing <a name="PreProcessing"></a>
[Nextflow](https://www.nextflow.io/) pipeline which runs image segmentation with [cellpose](https://github.com/MouseLand/cellpose) and then counts the transcripts in each cell. The pipeline uses two Python3 scripts for:
+ Segmentation
+ Expression assignment = counting the transcripts in each cell

These scripts can be used independently or as part of the Nextflow pipeline provided.

*Dependencies*
+ [Nextflow](https://www.nextflow.io/)
+ [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) 

The pipeline automatically fetches the following singularity container and uses it to run the scripts:

https://cloud.sylabs.io/library/michelebortol/resolve_tools/cellpose_skimage

The definition file is provided [here](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/singularity/cellpose.def).

## B.1) Nextflow pipeline <a name="##Pipeline"></a>

### B.1.1) Parameters <a name="##Parameters"></a>
For an example see the provided example config [file](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/example.config)
    
*Input/output Parameters:*
+ `params.input_path` = Path to the resolve folder with the Panoramas to be processed
+ `params.output_path` = Path for output

*cellpose Segmentation Parameters:*
+ `params.model_name` = "cyto" (recommended) or any model that uses 1 DNA channel.
+ `params.probability_threshold` = floating point number between -6 and +6 see [cellpose threshold documentation](https://cellpose.readthedocs.io/en/latest/settings.html#mask-threshold).
+ `params.cell_diameter` = Cell diameter or `None` for automatic estimation, see [cellpose diameter documentation](https://cellpose.readthedocs.io/en/latest/settings.html#diameter).
+ `params.do_zip` =	`true` or `false`.  Set to false to skip making ImageJ ROIs (faster)
+ `params.output_path` = "output/nextflow_test"

### B.1.2) Input <a name="##Input"></a>
Folder with the panoramas to be processed. All panoramas are expected to have:
+ DAPI image named: `Panorama_*_Channel3_R8_.tiff`
+ Transcripts coordinates named: `Panorama_*_results_withFP.txt`
    
### B.1.3) Output <a name="##Output"></a>
In `params.output_path`: 
+ `sample_metadata.csv`: .csv file with one row per sample and 3 columns: sample (sample name), dapi (path to the dapi image), counts (path to the transcript coordinates)
+ For each sample a folder: `SAMPLE_NAME` with: 
+ `SAMPLE_NAME-mask.tiff` = 16 bit segmentation mask (0 = background, N = pixels belonging to the Nth cell).
+ `SAMPLE_NAME-roi.zip` (optional) = ImageJ ROI file with the ROIs numbered according to the segmentation mask.
+ `SAMPLE_NAME-cell_data.csv` = Single cell data, numbered according to the semgentation mask.

### B.1.4) Example <a name="##Example"></a>
`nextflow run main.nf -profile cluster -c test.config`
Breakdown:
+ `-profile cluster` = For running on a PBS based cluster like the CURRY cluster. Default is local execution.
+ `-c test.config` = Use the parameters specified in the `test.config` file. ALternatively, parameters can be passed from the command line.


## B.2) Scripts <a name="#Scripts"></a>
Scripts used in the Nextflow pipeline, can also be run independently.

### B.2.1) Segmentation <a name="##Segmentation"></a>
[Segmentation script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/segmenter.py)
Just a wrapper around cellpose. It assumes the input is a single channel grayscale image with the nuclei. It requires the following positional arguments:
+ tiff_path = path to the image to segment
+ model_name = model to use for the segmentation			
+ prob_thresh = probability threshold
+ output_mask_file = path to the cell mask output
+ output_roi_file (optional) = path to the roi mask output or leave empty to skip (saves time).

**Example**  
`python3.9 segmenter.py DAPI_IMAGE cyto 0 70 OUTPUT_SEGMENTATION_MASK_NAME OUTPUT_ROI_ZIP_NAME`

### B.2.2) Expression assignment <a name="##expression_assign"></a>
[Expression assignment script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/segmenter.py)
Counts the transcripts in each cell from the segmentation mask. Equivalent to the Polylux counts unless:
+ Overlapping ROIs
+ Transcripts outside the border of the image or lying exactly on the ROI border (resolution is 1 pixel)
It requires the following positional arguments:
+ mask_file = Path to the input mask file
+ transcript_file = Path to the input transcript file
+ output_file = Path to the output single cell data file
Notes:
+ Removes all transcripts whose coordinates fall outside the size of the mask.

**Example**  
` python3.9 extracter.py SEGMENTATION_MASK TRANSCRIPT_COORDINATE_FILE OUTPUT_FILE_PATH.csv`

To Do:
+ Add option to filter by Z coordinate?
+ Add option to filter by transcript quality?
+ Add option to count transcripts from ROIs? 


# C) Analysis <a name="Analysis"></a>
Tentative scripts for data analyis.
## C.1) Seurat <a name="#Seurat"></a>
[Seurat example script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/resolveseurat.R)

The script is based on this [Seurat vignette](https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system).

The script currently requires a development version of Seurat from the [feat/imaging branch](https://github.com/satijalab/seurat/tree/feat/imaging)

To do:
+ Support for cells coming for more than one image in the same Seurat object.
+ Add image to Seurat Object? (Probably not worth it if we use Giotto)

## C.2) Giotto <a name="#Giotto"></a>

[Giotto from Seurat object example script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/resolvegiottoseurat.R)

The script is based on the [seuratToGiotto function](https://github.com/RubD/Giotto/blob/suite/R/interoperability.R).

The script currently requires a development version of Giotto from the [suite branch](https://github.com/RubD/Giotto/tree/suite)

To do:
+ Support for cells coming for more than one image in the same objet?
+ Add direct creation of a Giotto object without passing from Seurat?

Note:
+ The coordinate systems for cells are different between Seurat and Giotto (in Seurat X and Y are inverted). 
+ The transcript coordinates are also imported by `seuratToGiotto` but never used by the package.
