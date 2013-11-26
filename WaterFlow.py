import os, sys, time
from math import * 
from osgeo import gdal
from osgeo import gdalconst
from osgeo import ogr
from gdalconst import *

class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y

#
# Inputs:
#   raster - a handle to a raster from GDAL
#   start  - utm coordinate representing our start point, start.x and start.y should be filled
#
# Outputs:
#   flow - a list of points which make up our flow line
#
def WaterFlow(raster, start):
	
	# get raster information
	# width and height
	width = raster.RasterXSize
	height = raster.RasterYSize

	# get geographic information
	#top left x, pixel width, rotation, top left y, rotation, pixel height
	#notice pixel height is negative because it's measured from top!
	geoInfo = raster.GetGeoTransform()
	x0 = geoInfo[0]
	y0 = geoInfo[3]
	pWidth = geoInfo[1]
	pHeight = geoInfo[5]

	# load our data into an numpy array
	band = raster.GetRasterBand(1)
	data = band.ReadAsArray(0,0,width,height)

	# get the no data value for this raster
	NoData = band.GetNoDataValue()

	# determine our current point on the raster, based on the start point
	# convert the UTM coordinates to raster coordinates
	cX = int((start.x - x0) / pWidth)
	cY = int((start.y - y0) / pHeight)
	
	if(cX < 0 or cX > width):
		left = x0
		right = x0 + width * pWidth
		print(r'start.y/cX (' + str(start.x) + r'/' + str(cX) + r'is outside the bounds of the image: ' + str(left) + r' - ' + str(right))
		exit(1)

	if(cY < 0 or cY > height):
		top = y0
		bottom = y0 - height * pHeight
		print(r'start.y/cY (' + str(start.y) + r'/' + str(cY) + r'is outside the bounds of the image: ' + str(top) + r' - ' + str(bottom))
		exit(1)
	
	# declare a list to hold our flow points.
	# add the start point to the flow
	flow = []
	current = start
	flow.append(current)

	# Tuple holding vector and direction information
	#  7 0 1
	#  6   2
	#  5 4 3
	neighbors = ((0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1))
	direction = 0

	# precalculate the diagonal length
	diag = sqrt((pWidth * pWidth) + (pHeight * pHeight))

	# determine the start elevation
	cElevation = data[cX,cY]

	# loop while we are not on an edge, once we reach an edge we will quit
	while((cX > 0 and cX < width) or (cY > 0 and cY < height)):

		# initialize variables for this pass
		bestSlope = 0.0
		bestPixel = None

		# select the elevations from the set of neighbors
		nElevations = SelectNeighbors(data, neighbors, cX, cY, width, height)

		# this loop is to guarantee that we do not get stuck in depressions
		# that they fill as water enters them
		# this is a naive approach and possibly needs to be reconsidered.
		while(bestPixel is None):

			# loop through the neighbors, need a little trickery to more around the array
			for dir in range(direction, direction + 8, 1):

				# calculate some index fanciness
				# allows us to loop around our neighbors array
				i = dir % 8

				nElevation = nElevations[i]

				# if we have NoData in this cell, continue the loop
				if(nElevation == NoData): continue

				# calculate grade using rise/run
				# need to get a little more complex in case the image does not have square cells
				run = diag if (i % 2 == 1) else pWidth
				run = pHeight if (i % 4 == 0) else run

				# neighbor elevation - cElevation
				rise = nElevation - cElevation
				slope = float(rise / elevation)

				# we have a candidate from the neighbors
				# either the slope is lower than the best slope
				# or we have the same slope as the best, but we are heading in the last known direction, inertia?
				if(slope < bestSlope or (i == direction and slope == bestSlope)) : 
					bestSlope = slope
					
					# convert back to real world coordinates and store the pixel in teh flow
					bestPixel = Shape(current.x + float(neighbors[i][0] * pWidth), current.y + float(neighbors[i][1] * pHeight))
					bestPixel.x = current.x + float(neighbors[i][0] * pWidth)
					bestPixel.y = current.y + float(neighbors[i][1] * pHeight)

					# set the current direction
					direction = i

			# if we still have a null bestPixel, need to raise the elevation by one and continue.
			# we are slowly filling the depression with water, until it spills over
			if(bestPixel is null):
				data[cX, cY] += 1.0

		# we have our best pixel
		flow.append(bestPixel)

		#update cX and cY based on neighbor
		cX = cX + neighbors[direction][0]
		cY = cY + neighbors[direction][1]

	return flow


#
#	Inputs
#		data - elevation values from DEM
#		neighbors - tuples of vectors representing neighbor pixels
#		cX - current x position in raster
#		cY - current y position in raster
#		width - raster width
#		height - raster height
#
#	Outputs:
#		List of elevation values for the surrounding pixels
#
def SelectNeighbors(data, neighbors, cX, cY, width, height):

	elevations

	for i in range[0,7]:
		x = cX + neighbors[i][0]
		y = cY + neighbors[i][1]

		# make sure we aren't out of range
		if(x < 0 or x > width) : continue
		if(y < 0 or y > height) : continue

		elevations.append(data[x,y])

	return elevations

def PrintUsage():
	print(r'WaterFlow is a simulation which maps how water could flow across a DEM')
	print()
	print('Usage:')
	print('\tWaterFlow [pathToDEM] [UTM Easting] [UTM Northing]')


#
#   [1]: path to DEM to use
#   [2]: start x/easting/longitude value
#   [3]: stary y/northing/latitude value
#
def main():
		
	if(len(sys.argv) != 4): 
		PrintUsage()
		exit(-1)

	rasterHandle = None
	flowLine = None
	startPoint = Point(int(sys.argv[2]), int(sys.argv[3]))
		
	# if file exists
	if(os.path.exists(sys.argv[1])):

		# read the file into a 2 dimensional array
		rasterHandle = gdal.Open(sys.argv[1], GA_ReadOnly)

		geoInfo = rasterHandle.GetGeoTransform()
		if not geoInfo is None:
			print 'Origin = (', geoInfo[0], ',',geoInfo[3],')'
			print 'Pixel Size = (',geoInfo[1],',',geoInfo[5],')'
	
	else:        
		print(sys.argv[1] + ' does not exist! Exiting..')
		sys.exit(-1)

	# we now have our elevation and start point, call water flow
	flowLine = WaterFlow(rasterHandle, startPoint)

	# output the flowline to some format...
	

if __name__ == '__main__':
	main()



