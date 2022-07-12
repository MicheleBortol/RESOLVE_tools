import argparse
import numpy as np
import cv2
from cellpose import models, io, utils
from roifile import ImagejRoi, roiwrite 

def get_arguments():	
	"""
	Parses and checks command line arguments, and provides an help text.
	Assumes 5 and returns 5 positional command line arguments:
	tiff_path = Path to the tiff file
	model_name = Model to use for the segmentation
	prob_thresh = Probability threshold
	cell_diameter = Expected cell diameter
	output_mask_file = Path to output the cell mask
	output_zip = path to the zip file for the ROIs (OPTIONAL)
	"""
	parser = argparse.ArgumentParser(description = "Performs 2D segmentation with cellpose.")
	parser.add_argument("tiff_path", help = "path to the image to segment")
	parser.add_argument("model_name", help = "model to use for the segmentation")			
	parser.add_argument("prob_thresh", help = "probability threshold")
	parser.add_argument("cell_diameter", help = "expected cell diameter")
	parser.add_argument("output_mask_file", help = "path to the cell mask output")
	parser.add_argument("output_zip", nargs = "?", help = "path to the zip file for the ROIs",
			default = None)
	args = parser.parse_args()
	return args.tiff_path, args.model_name, args.prob_thresh, args.cell_diameter, \
		args.output_mask_file, args.output_zip

if __name__ == "__main__":
		tiff_path, model_name, prob_thresh, cell_diameter, \
			output_mask_file, output_zip = get_arguments()

		# Define cellpose model
		model = models.Cellpose(gpu = False, model_type = model_name)

		channels = [0, 0] # We assume the input is a single grayscale image

		try:
			cell_diameter = float(cell_diameter)
		except ValueError:
			cell_diameter = None
		# if diameter is set to None, the size of the cells is estimated on a per image basis
		# you can set the average cell `diameter` in pixels yourself (recommended) 
		# diameter can be a list or a single number for all images

		# Load the input image
		img = io.imread(tiff_path)
		
		# Apply the model
		mask, flows, style, diameter = model.eval(img, diameter = cell_diameter,
				channels = channels) # recover diameter if it was set to None earlier

		# save masks
		io.imsave(output_mask_file, mask)
	
		# extract outlines and make ROIs
		if output_zip != None:
			mask = np.rot90(mask, k = -1)
			mask = np.fliplr(mask)
			mask_shape = mask.shape
			edges = np.where(utils.masks_to_outlines(mask), mask, 0)
			mapper = lambda n: np.column_stack(np.unravel_index(np.flatnonzero(edges == n),
				mask_shape))
			cells = map(mapper, range(1, np.amax(mask, axis=None) + 1))
			rois = list(map(ImagejRoi.frompoints, cells))
			roi_ids = list(map(str, range(0, len(rois))))
			roiwrite(output_zip, rois, roi_ids, mode = "w")
