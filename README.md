# RESOLVE_tools
Tentative set of tools and scripts for analysing spatial transcriptomic data with the resolve platform

## Contents:
+ A) [PreProcessing](#PreProcessing)
    + (A.1) [Segmentation](##Segmentation)
    + (A.2) [Expression assignment](##expression_assign)  bin/segmenter.py = Just a wrapper around cellpose.
+ B) [Analysis](#Analysis)
  + (B.1) [Seurat](##Seurat)
  + (B.1) [Giotto](##Giotto)


# A) PreProcessingg <a name="PreProcessing"></a>
Python3 scripts for:
+ Segmentation
+ Expression assignment = counting the transcripts in each cell


## A.1) Segmentation <a name="#Segmentation"></a>
[Segmentation script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/segmenter.py)
Just a wrapper around cellpose. It assumes the input is a single channel grayscale image with the nuclei. It requires the following positional arguments:
+ sample_name = name of the sample (currently required but unused, would be useful if we make a Nextflow pipeline from these scripts)
+ tiff_path = path to the image to segment
+ model_name = model to use for the segmentation			
+ prob_thresh = probability threshold
+ output_mask_file = path to the cell mask output
+ output_roi_file = path to the roi mask output

To Do:
+ Add optional parameter for diameter selection. Currently the diameter is estimated from the image.

## A.2) Expression assignment <a name="#expression_assign"></a>
[Expression assignment script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/segmenter.py)
Counts the transcripts in each cell from the segmentation mask. Equivalent to the Polylux counts unless:
+ Overlapping ROIs
+ Transcripts outside the border of the image or lying exactly on the ROI border (resolution is 1 pixel)
It requires the following positional arguments:
+ sample_name = name of the sample (currently required but unused, would be useful if we make a Nextflow pipeline from these scripts)
+ mask_file = Path to the input mask file
+ transcript_file = Path to the input transcript file
+ output_file = Path to the output single cell data file
Notes:
+ Removes all transcripts whose coordinates fall outside the size of the mask.

To Do:
+ Add option to filter by Z coordinate?
+ Add option to filter by transcript quality?
+ Add option to count transcripts from ROIs? 


# B) Analysis <a name="Analysis"></a>
Tentative scripts for data analyis.
## B.1) Seurat <a name="#Seurat"></a>
[Seurat example script](https://github.com/MicheleBortol/RESOLVE_tools/blob/main/bin/resolveseurat.R)

The script is based on this [Seurat vignette](https://satijalab.org/seurat/articles/spatial_vignette_2.html#human-lymph-node-akoya-codex-system).

The script currently requires a development version of Seurat from the [feat/imaging branch](https://github.com/satijalab/seurat/tree/feat/imaging)

To do:
+ Support for cells coming for more than one image in the same Seurat object.
+ Add image to Seurat Object? (Probably not worth it if we use Giotto)

## B.2) Giotto <a name="#Giotto"></a>


