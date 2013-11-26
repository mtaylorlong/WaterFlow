import os, sys, time, math
from math import * 
from osgeo import gdal
from osgeo import gdalconst
from osgeo import ogr
from gdalconst import *

class Point:
	def __init__(self, x, y):
		self.x = float(x)
		self.y = float(y)

#
# Converts a pixel coordinate to world coordinate
# 
# Inputs:
#	pixel - a point object with pixel coordinates
#	origin - the origin coordinate for raster (0,0)
#	pWidth - the width of the pixel in world 
#	pHeight - the height of the pixel in world
#
# Output:
#	A point object in world coordinates
#
def PixelToWorld(pixel, origin, pWidth, pHeight):

	# calculate the X, Y
	x = float(pixel.x * pWidth)
	y = float(pixel.y * abs(pHeight))

	# add the X to the origin and subtract the Y from the origin
	world = Point(origin.x + x, origin.y - y)

	return world

#
# Converts a world coordinate to pixel coordinate
# 
# Inputs:
#	world - a point object with world coordinates
#	origin - the origin coordinate for raster (0,0)
#	pWidth - the width of the pixel in world 
#	pHeight - the height of the pixel in world
#
# Output:
#	A point object in pixel coordinates
#
def WorldToPixel(world, origin, pWidth, pHeight):
	
	# calculate the x, y differences
	x = float(world.x - origin.x)
	y = float(origin.y - world.y)

	# determine how many pixels we are
	xPixels = x / pWidth
	yPixels = y / abs(pHeight)

	# return the point, use int() to truncate to pixel coordinate
	return Point(int(xPixels), int(yPixels))

#
# Converts degree measurements to meters
# Inputs
#   pWidth - The pixel width in degrees
#   pHeight - The pixel height in degrees
#   lat - The latitude to perform this conversion
#
# Outputs
#   The measurement in meters [x,y]  
#
# Notes: Uses Spherical Law of Cosines found at http://www.movable-type.co.uk/scripts/latlong.html
def DegreesToMeters(pWidth, pHeight, lat):
	
	# http://en.wikipedia.org/wiki/World_Geodetic_System
	radiusInMeters = 6378137.0	

	# convert the lat/width/height to radians
	rLat = float(math.pi * lat / 180)
	rWidth = float(math.pi * pWidth /180)
	rHeight = abs(float(math.pi * pHeight / 180))

	# convert the radians to meters
	x = rWidth * math.cos(rLat) * radiusInMeters
	y = rHeight * radiusInMeters

	return [x,y]
	
#
# Outputs the flowline to KML
#	inputs:
#		flow - a list representing the coordinates of the flowline
#
#	outputs:
#		none, but writes out the geometry to a KML file
#
def OutputFlowline(flow):

	# the generic KML document, we modify by replacing "MyGeometry" and "myKML" with teh appropriate data
	kmlHeader = r'<?xml version="1.0" encoding="UTF-8"?><kml xmlns="http://www.opengis.net/kml/2.2"><Placemark><name>MyGeometry</name>myKML</Placemark></kml>'
	
	# provide a name for this geometry, then replace
	name = str(len(flow))+'points'
	kmlFile = kmlHeader.replace('MyGeometry', name)
	
	# create the linestring using OGR, then add the points
	line = ogr.Geometry(ogr.wkbLineString)

	for point in flow:
		line.AddPoint(point.x, point.y)

	# convert to KML, then put in the KML document
	kml = line.ExportToKML()
	kmlFile = kmlFile.replace('myKML', kml)

	# output the text to a ".kml" file, the name is # of points
	name = name + '.kml'
	filename = 'H:\\GEOG 5562\\output\\' + name

	if(os.path.exists(filename)):
		os.remove(filename)

	f = open(filename, 'w')
	f.write(kmlFile)
	f.close()	

#
# Inputs:
#   raster - a handle to a raster from GDAL
#   start  - pixel coordinate representing our start point, start.x and start.y should be filled
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
	origin = Point(geoInfo[0], geoInfo[3])
	pWidth = float(geoInfo[1])
	pHeight = float(geoInfo[5])

	# load our data into an numpy array
	band = raster.GetRasterBand(1)
	data = band.ReadAsArray(0,0,width,height)

	# get the no data value for this raster
	NoData = band.GetNoDataValue()

	# determine our current point on the raster, based on the start point
	# convert the world coordinates to raster coordinates
	cX = int((start.x - x0) / pWidth)
	cY = int((start.y - y0) / pHeight)
	
	if(start.x < 0 or start.x > width):
		left = x0
		right = x0 + width * pWidth
		print(r'start.y/cX (' + str(start.x) + r'/' + str(cX) + r')is outside the bounds of the image: ' + str(left) + r' - ' + str(right))
		exit(1)

	if(start.y < 0 or start.y > height):
		top = y0
		bottom = y0 - height * pHeight
		print(r'start.y/cY (' + str(start.y) + r'/' + str(cY) + r')is outside the bounds of the image: ' + str(top) + r' - ' + str(bottom))
		exit(1)
	
	# need to determine the pixel size in meters, since we are using lat/lon
	# convert the current point to a world coordinate
	currentWorld = PixelToWorld(start, origin, pWidth, pHeight)
	meters = DegreesToMeters(pWidth, pHeight, currentWorld.y)

	# declare a list to hold our flow points.
	# add the start point to the flow
	
	flow = []
	flow.append(currentWorld)

	# Tuple holding vector and direction information
	#  7 0 1
	#  6   2
	#  5 4 3
	neighbors = ((0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1))
	direction = 0

	# precalculate the diagonal length and pixel size in meters
	mWidth = meters[0]
	mHeight = meters[1]
	diag = sqrt((mWidth * mWidth) + (mHeight * mHeight))	
	
	maxX = width - 1
	maxY = height - 1

	# set cX and cY to the current
	cX = int(start.x)
	cY = int(start.y)
		
	# start a dictionary to ensure we do not double back on pixels
	visitedPixels = dict()	
	
	# debugging
	#pixelCount = 0
	
	# loop while we are not on an edge, once we reach an edge we will quit
	while((cX > 0 and cX < maxX) or (cY > 0 and cY < maxY)):

		# debugging
		#pixelCount+= 1
		#if(pixelCount >=58):
		#	OutputFlowline(flow)			
		#print 'Searching at pixel ', cX, ',', cY

		# update this pixel as visited
		visitedPixels[(cX,cY)] = True

		# initialize variables for this pass
		bestSlope = 0.0
		bestPixel = None

		# determine the start elevation
		cElevation = float(data.item(cY,cX).real)

		# select the pixels from the set of neighbors
		neighborPixels = SelectNeighbors(neighbors, cX, cY, maxX, maxY)

		# this loop is to guarantee that we do not get stuck in depressions
		# that they fill as water enters them
		# this is a naive approach and possibly needs to be reconsidered.
		while(bestPixel is None):

			# loop through the neighbors, need a little trickery to move around the array
			directions = range(direction, direction + 8, 1)
			
			for dir in directions:

				# calculate some index fanciness
				# allows us to loop around our neighbors array
				i = dir % 8

				# pull the neighboring pixel and its elevation data
				neighbor = neighborPixels[i]                
				nElevation = float(data.item(neighbor[1], neighbor[0]).real)

				# if we have NoData in this cell or we have (-1,-1) or we have visited the pixel, continue the loop
				if((neighbor in visitedPixels) or nElevation == NoData or neighbor == (-1, -1)): continue
								
				# calculate grade using rise/run
				# need to get a little more complex in case the image does not have square cells
				run = diag if (i % 2 == 1) else mWidth
				run = mHeight if (i % 4 == 0) else run

				# neighbor elevation - cElevation
				rise = float(nElevation - cElevation)
				slope = float(rise / run)

				# we have a candidate from the neighbors
				# either the slope is lower than the best slope
				# or we have the same slope as the best, but we are heading in the last known direction, inertia?
				if(slope <= bestSlope or (i == direction and slope <= bestSlope)) : 
					bestSlope = slope
					
					# update the best pixel
					bestPixel = Point(cX + neighbors[i][0], cY + neighbors[i][1])
					print '\t Best Pixel: ', bestPixel.x, ',', bestPixel.y, ' in direction: ', i

					# set the current direction
					direction = i

			# if we still have a null bestPixel, we are stuck in a depresson
			# Find the next lowest pixel that is and set that regardless of visitation
			if(bestPixel is None):
				
				# debugging
				#dWorld = PixelToWorld(Point(cX, cY), origin, pWidth, pHeight)
				#print '\tDepression at ', dWorld.x, ' (',cX,'), ', dWorld.y, ' (',cY,'), Elev: ', cElevation
				
				# find the next lowest elevation
				nextBestIndex = -1
				nextLowestHeight = 500000

				for d in range(0,8):
					
					neighbor = neighborPixels[d]
					elev = data[neighbor[1], neighbor[0]]

					if(elev >= cElevation and elev < nextLowestHeight):
						nextBestIndex = d
						nextLowestHeight = elev
						
				# create the best pixel from the next best.
				bestPixel = Point(neighborPixels[nextBestIndex][0], neighborPixels[nextBestIndex][1])					

				# debugging
				#dWorld = PixelToWorld(bestPixel, origin, pWidth, pHeight)
				#print '\tNext Best : ', dWorld.x, ' (',bestPixel.x,'), ', dWorld.y, ' (',bestPixel.y,'), Elev: ', nextLowestHeight				
				

		# we have our best pixel, convert to world and add
		bestWorld = PixelToWorld(bestPixel, origin, pWidth, pHeight)
		flow.append(bestWorld)
		
		#debugging
		#print '\tAdded ', bestWorld.x, ',', bestWorld.y , ' to the flow line in direction ', direction 		

		#update cX and cY based on neighbor
		cX = int(neighborPixels[direction][0])
		cY = int(neighborPixels[direction][1])

		# limit how often we update/output lines
		if(len(flow) % 500 == 0):
			OutputFlowline(flow)

	return flow

#
#	Inputs
#		neighbors - tuples of vectors representing neighbor pixels
#		cX - current x position in raster
#		cY - current y position in raster
#		width - raster width
#		height - raster height
#
#	Outputs:
#		List of pixel coordinates for the surrounding pixels
#
def SelectNeighbors(neighbors, cX, cY, width, height):

	elevations = []

	for i in range(0,8):
		x = int(cX + neighbors[i][0])
		y = int(cY + neighbors[i][1])

		# make sure we aren't out of range
		if(x < 0 or x > width or y < 0 or y > height) : 
			elevations.append((-1, -1))
			continue	
		
		pixel = (x,y)
		elevations.append(pixel)

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
	startPoint = Point(float(sys.argv[2]), float(sys.argv[3]))
		
	origin = None
	pWidth = None
	pHeight = None

	# if file exists
	if(os.path.exists(sys.argv[1])):

		# read the file into a 2 dimensional array
		rasterHandle = gdal.Open(sys.argv[1], GA_ReadOnly)

		geoInfo = rasterHandle.GetGeoTransform()
		origin = Point(geoInfo[0], geoInfo[3])
		pWidth = geoInfo[1]
		pHeight = geoInfo[5]

		if not geoInfo is None:
			print 'Origin = (', geoInfo[0], ',',geoInfo[3],')'
			print 'Pixel Size = (',geoInfo[1],',',geoInfo[5],')'
	
	else:        
		print(sys.argv[1] + ' does not exist! Exiting..')
		sys.exit(-1)
	
	startPixel = WorldToPixel(startPoint, origin, pWidth, pHeight)
	print r'Beginning analysis at ', startPoint.x, ' ', startPoint.y
	print '\tPixel coordinate: ', startPixel.x, ' ', startPixel.y

	# we now have our elevation and start point, call water flow
	try:
		flowLine = WaterFlow(rasterHandle, startPixel)
	except ValueError:
		print ValueError

	# output the flowline to some format...
	OutputFlowline(flowLine)
	

if __name__ == '__main__':
	main()



