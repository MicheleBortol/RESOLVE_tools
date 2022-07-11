script_folder = "$baseDir/bin"
process collect_data{

    memory { 1.GB }
    time '1h'
    
    publishDir "$params.output_path", mode:'copy', overwrite: true

    input:
		val(input_path)
	
    output:
        path("sample_metadata.csv", emit: metadata_csv)
    
    script:
    """
	echo "sample,dapi,counts" > sample_metadata.csv
	while IFS= read -d \$'\\0' -r DAPI
	do
		SAMPLE="\${DAPI##$input_path/Panorama_}"
		SAMPLE="\${SAMPLE%%_Channel3_R8_.tiff}"
		COUNTS="$input_path/Panorama_""\$SAMPLE""_results_withFP.txt"
		echo "\$SAMPLE,\$DAPI,\$COUNTS" >> sample_metadata.csv 
	done < <(find "$input_path/" -name "*_Channel3_R8_.tiff" -print0)
    """
}

process cellpose_segment{
    
    memory { 128.GB * task.attempt }
    time '72h'
    
    errorStrategy { task.exitStatus in 137..143 ? 'retry' : 'terminate' }
    maxRetries 3

    publishDir "$params.output_path/$sample_name", mode:'copy', overwrite: true
    container = "library://michelebortol/resolve_tools/cellpose_skimage:latest"

    input:
		val(sample_name)
		val(model_name)
		val(probability)
		val(diameter)
		val(do_zip)
		path(dapi_path)
	
    output:
        path("$sample_name-mask.tiff", emit: mask_image)
        path("$sample_name-roi.zip", emit: roi_zip, optional: true)

    script:
	  def zip_name = do_zip ? "$sample_name-roi.zip" : ''
    """
	python3.9 $script_folder/segmenter.py $dapi_path $model_name $probability \
		$diameter $sample_name-mask.tiff $zip_name

    """
}

process  extract_sc_data{

    memory { 16.GB * task.attempt }
    time '72h'

    publishDir "$params.output_path/$sample_name", mode:'copy', overwrite: true
    container = "library://michelebortol/resolve_tools/cellpose_skimage:latest"

    input:
        val(sample_name)
		path(mask_image_path)
		path(transcript_coord_path)

    output:
        path("$sample_name-cell_data.csv", emit: sc_data)

    script:

    """
	python3.9 $script_folder/extracter.py $mask_image_path $transcript_coord_path \
		${sample_name}-cell_data.tsv	
    """
}

