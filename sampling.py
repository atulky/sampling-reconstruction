import random
import argparse
import time
import scipy.interpolate
from vtk import *
import numpy as np

# SIMPLE RANDOM SAMPLING 
def srs(data, percent):
	# Sampled output
	sample = []
	sample_values = []

	# selecting the 8 corner points
	xmax, ymax, zmax = data.GetDimensions()
	sample += [(0,0,0), 
		       (xmax-1,0,0), 
		       (0,ymax-1,0), 
		       (0,0,zmax-1), 
		       (xmax-1,ymax-1,0), 
		       (0,ymax-1,zmax-1), 
		       (xmax-1,0,zmax-1), 
		       (xmax-1,ymax-1,zmax-1)]

	# save the index of already sampled points
	sampled_indices = set([data.FindPoint(point) for point in sample])
	
	# extract scalar values of data
	scalar_values = data.GetPointData().GetScalars()
	
	# get values at the 8 corner points
	sample_values += [scalar_values.GetValue(i) for i in sampled_indices]

	# calculate number of samples to select
	total_points = data.GetNumberOfPoints()
	n_samples = int(total_points * (percent/100)) - 8

	# select points randomly
	while len(sampled_indices) < n_samples:
		sampled_indices.add(random.randint(0, total_points-1))
		
	# save point coordinates and corresponding values
	sample += [list(data.GetPoint(i)) for i in sampled_indices]
	sample_values += [scalar_values.GetValue(i) for i in sampled_indices]
	   
	return sample, sample_values

# COMPUTE SNR
def compute_SNR(arrgt, arr_recon):
	diff = arrgt - arr_recon
	sqd_max_diff = (np.max(arrgt)-np.min(arrgt))**2
	snr = 10*np.log10(sqd_max_diff/np.mean(diff**2))
	return snr

# READ VTI FILE
def readVTI(filepath):
	reader = vtkXMLImageDataReader()
	reader.SetFileName(filepath)
	reader.Update()
	data = reader.GetOutput()
	return data

# WRITE VTP FILE FOR SAMPLED POINTS
def writeVTP(sample, sample_values, filename):
	# Create VTK points and float array for points and values respectively
	points = vtkPoints()
	values = vtkFloatArray()
	values.SetNumberOfValues(len(sample_values))
	
	for point in sample:
		points.InsertNextPoint(point)
		
	for i, value in enumerate(sample_values):
		values.SetValue(i, value)

	pdata = vtkPolyData()
	pdata.SetPoints(points)
	pdata.GetPointData().SetScalars(values)

	# Write to .vtp file
	writer = vtkXMLPolyDataWriter()
	writer.SetInputData(pdata)
	writer.SetFileName(filename)
	writer.Write()

# WRITE RECONSTRUCTED DATA TO VTI FILE
def writeVTI(recon_values, filename):
	recon_data = vtkImageData()
	recon_data.SetDimensions(250,250,50)
	recon_data.SetSpacing(1.0,1.0,1.0)
	recon_data.SetOrigin(0.0,0.0,0.0)
	recon_data.AllocateScalars(VTK_FLOAT, 1)
	
	for i in range(250):
		for j in range(250):
			for k in range(50):
				recon_data.SetScalarComponentFromFloat(i, j, k, 0, recon_values[i][j][k])
		
	# Write to .vti file
	writer = vtkXMLImageDataWriter()
	writer.SetInputData(recon_data)
	writer.SetFileName(filename)
	writer.Write()

if __name__=='__main__':
	# Take inputs for sampling percentage and reconstruction method
	parser = argparse.ArgumentParser()
	parser.add_argument('--percent', required=True, type=float, help='sampling percentage')
	parser.add_argument('--method', required=True, type=str, help='reconstruction method')
	
	options = parser.parse_args()
	
	# Read data file
	data = readVTI("Isabel_3D.vti")
	
	# Perform sampling
	sample, sample_values = srs(data, options.percent)
	
	# Write the sampled points into sample.vtp file
	writeVTP(np.array(sample), np.array(sample_values), "sample.vtp")
	
	# Reconstruction using sampled points
	gx, gy, gz = np.mgrid[0:250, 0:250, 0:50] # create 3D grid
	print(f"starting {options.method} interpolation...")
	tic = time.time() # start clock
	recon_values = scipy.interpolate.griddata(sample, sample_values, (gx, gy, gz), method=options.method)
	print(f"finished {options.method} interpolation")
	
	# Nearest neighbor interpolation if nan values are present
	nan_vals = np.isnan(recon_values)
	nan_count = np.count_nonzero(nan_vals)
	if nan_count != 0:
		print("handling NaN values ...")
		nan_idx = np.argwhere(~nan_vals)
		recon_values = scipy.interpolate.griddata(nan_idx, recon_values[~nan_vals], (gx, gy, gz), method='nearest')
		print(f"{nan_count} NaN values were replaced with nearest neighbor value")
		
	toc = time.time() # stop clock
	print(f"reconstruction time: {round(toc - tic, 2)} seconds")
		
	
	# Data scalar values
	scalar_values = data.GetPointData().GetScalars()
	
	# Get ground truth values and reconstructed values
	arrgt, arr_recon = [], []
	for i in range(250):
		for j in range(250):
			for k in range(50):
				arrgt.append(scalar_values.GetValue(data.FindPoint(i, j, k)))
				arr_recon.append(recon_values[i][j][k])
	
	# compute SNR
	print(f"SNR: {round(compute_SNR(np.array(arrgt), np.array(arr_recon)), 2)}")
	
	# Write reconstructed values into {method}_recon.vti
	print(f"writing reconstructed data to {options.method}_recon.vti ...")
	writeVTI(recon_values, f"{options.method}_recon.vti")
	print("done")
