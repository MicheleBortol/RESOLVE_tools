import argparse
from cellpose import models, io, utils

def get_arguments():	
	"""
	Parses and checks command line arguments, and provides an help text.
	Assumes 6 and returns 6 positional command line arguments:
	sample_name = Name of the sample
	tiff_path = Path to the tiff file
	model_name = Model to use for the segmentation
	prob_thresh = Probability threshold
	output_mask_file = Path to output the cell mask
	output_roi_file = Path to output the ROI file
	"""
	parser = argparse.ArgumentParser(description = "Performs 2D segmentation with cellpose.")
	parser.add_argument("sample_name", help = "name of the sample")
	parser.add_argument("tiff_path", help = "path to the image to segment")
	parser.add_argument("model_name", help = "model to use for the segmentation")			
	parser.add_argument("prob_thresh", help = "probability threshold")
	parser.add_argument("output_mask_file", help = "path to the cell mask output")
	parser.add_argument("output_roi_file", help = "path to the roi mask output")
	args = parser.parse_args()
	return args.sample_name, args.tiff_path, args.model_name, 
			args.prob_thresh, args.output_mask_file, args.output_roi_file

if __name__ == "__main__":
		sample_name, tiff_path, model_name, prob_thresh, 
		output_mask_file, output_roi_file = get_arguments()

		# Define cellpose model
		model = models.Cellpose(gpu = False, model_type = model_name)

		channels = [0, 0] # We assume the input is a single grayscale image

		diameter = None
		# if diameter is set to None, the size of the cells is estimated on a per image basis
		# you can set the average cell `diameter` in pixels yourself (recommended) 
		# diameter can be a list or a single number for all images

		# Load the input image
		img = io.imread(tiff_path)
		
		# Apply the model
		mask, flows, style, diameter = model.eval(img, diameter = diameter,
				channels = channels) # recover diameter if it was set to None earlier

		# extract outlines = cell ROIs
		outlines = utils.outlines_list(mask)
		
		# save masks
		io.imsave(output_mask_file, mask)

		# save outlines
		io.outlines_to_text(output_roi_file, outlines)
