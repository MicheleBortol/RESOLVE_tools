nextflow.enable.dsl=2

script_folder = "$baseDir/bin"

include {data_collection} from "$script_folder/workflows.nf"
include {segmentation} from "$script_folder/workflows.nf"
include {sc_data_extraction} from "$script_folder/workflows.nf"

workflow {
	data_collection(params.input_path)
    sample_metadata = data_collection.out.data_csv \
	    | splitCsv(header:true) \
		| multiMap { row-> 
			sample: row.sample
			dapi: row.dapi
			counts: row.counts}

	segmentation(sample_metadata.sample, params.model_name, \
		params.probability_threshold, params.cell_diameter, \
		params.do_zip, sample_metadata.dapi)
	cell_masks = segmentation.out.mask_images
	roi_zips = segmentation.out.roi_zips

	sc_data_extraction(sample_metadata.sample, cell_masks, \
		sample_metadata.counts)
}
