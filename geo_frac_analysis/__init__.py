""" GeoFracAnalysis
	===============
	
	A module for spatial analysis of geological fractures.

	Details on the type of analysis can be found here:	

	Classes Included
	-----------------------------------------------------------------
		FracAnalysisPoly
		- Class for spatially analysing geological fracture data from shapefiles
		 if the data is in polylines.
		
		FracAnalysisPoint
		- Like FracAnalysisPoly except for shapefiles containing point data.		

		Both of these classes analyse fracture data by spatial distribution, 
		lengths (in the case of FracAnalysisPoly)  and angle, all based on user 
		defined cell size and angle bins.

		Can give Number density, Length Density, Number Anisotropy, Group Dominance
		Frequency, and can write out to a shapefile.

	Functions Included
	-----------------------------------------------------------------

		FancyPlot
		- For plotting fractures on a user-defined grid.

		FancyPlotTotals
		- For plotting groups of separate shapefiles.

	Usage
	-----------------------------------------------------------------
	Example of use for single shapefile containing polylines:
	
	>>> import GeoFracAnalysis as GFA

	>>> shapefile_address = "wexfaults.dbf"

	For cell size of 5km, and angle bins of 60 degrees (180/3):

	>>> cell_size, angle_bins = 5, 3

	>>> a = FracAnalysisPoly(shapefile_address, cell_size, angle_bins)

	>>> FancyPlot(a, Rose = True, Fractures = True, Patches = "NumberAnisotropy", Circles = True, SquareNumbers = True, FigureNumber = 1)

		To save the output to a shapefile:

	>>> a.save_output("address_for_output")

	To analyse and plot an additional region contained in a shapefile:

	>>> shapefile_address = "Galwayfaults.dbf"

	>>> b = FracAnalysisPoly(shapefile_address2, cell_size, angle_bins)

	>>> FancyPlot([a, b], Rose = True, Fractures = True, Patches = "NumberAnisotropy", Circles = True, SquareNumbers = True, FigureNumber = 1)


	Author
	-------------------

	Written by Sean Blake in Trinity College Dublin, 2014-2016.

	Email: blakese@tcd.ie

	GITHUB: https://github.com/TerminusEst

	Uses the MIT license.

"""

import shapefile

from matplotlib import pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

import numpy as np

################################################################################
################################################################################
################################################################################
################################################################################
# clockwise angle between two points

def angle_calc(x1, y1, x2, y2):
	""" clockwise angle between two points """
	x = x1 - x2
	y = y1 - y2

	anglez = np.degrees(np.arctan2(x, y))
	if anglez < 0:
		anglez += + 180

	if anglez == 180.0:
		anglez -= 0.01
	return anglez

################################################################################
################################################################################
# LIANGBARSKY THING

def liangbarsky(left, top, right, bottom, x1, y1, x2, y2):
	""" Calculate length of polyline in a grid

		Parameters
		-----------
		left, top, right, bottom 	= bounding co-ordinates of box
		x1, y1, x2, y2		= number of points in the grid

		Returns
		-----------
		length of polyline in box
	
		-----------------------------------------------------------------
	"""
	dx = x2 - x1
	dy = y2 - y1
	dt0, dt1 = 0, 1

	checks = ((-dx, -(left - x1)), (dx, right - x1), (-dy, -(bottom - y1)), (dy, top - y1))

	for p, q in checks:
		if p == 0:
			p = 0.01
			if q < 0:
				return 0, 0, 0, 0, 0
		dt = q / (p * 1.0)
		if p < 0:
			if dt > dt1:
				return 0, 0, 0, 0, 0
			dt0 = max(dt0, dt)
		else:
			if dt < dt0:
				return 0, 0, 0, 0, 0
			dt1 = min(dt1, dt)

	x2 = x1 + (dt1 * dx)
	y2 = y1 + (dt1 * dy)
	x1 = x1 + (dt0 * dx)
	y1 = y1 + (dt0 * dy)

	distance = np.hypot(x2 - x1, y2 - y1)
	return distance, x1, y1, x2, y2

################################################################################
################################################################################
################################################################################

class FracAnalysisPoly:
	"""
		Return analyzed shapefile in a grid.

		Parameters
		-----------
		address : string
			location of shapefile with only polylines

		cell_size : int, float
			size of grid square sides in km

		angle_divs : int
			number of angle divisions in 180 degrees

		Attributes
		-----------
		X, Y	: numpy.ndarray
			arrays of x and y-coordinates of grids with non-zero points

		N : numpy.ndarray
			array of number of counts in square, divided by angle bin

		N_total : numpy.ndarray
			array of true number of fractures per square (excludes double counting
			of fractures which are in >1 angle bin)

		L : numpy.ndarray
			array of length of fractures in square in km, divided by angle bin

		Number_Anisotropy : numpy.ndarray
			array of number Anisotropy (1 value per square)

		Functions
		-----------
		save_output:
			saves attributes to a shapefile
		-----------------------------------------------------------------
	"""
	def __init__(self, address, cell_size, angle_divs):
		self.address = address
		self.cell_size = cell_size
		self.angle_divs = angle_divs

		# Read in the data
		sf = shapefile.Reader(address)
		shapes = sf.shapes()

		x_coord, y_coord = [], []
		for i in shapes[::]:
			a = np.array(i.points)
			x_coord.append(a[:,0])
			y_coord.append(a[:,1])

		if type(angle_divs) != int:
			print "\n The argument 'angle_divs' needs to be int\n"
			return

		km = 1000

		# get bounds for plot
		maxx = max([item for sublist in x_coord for item in sublist])
		minx = min([item for sublist in x_coord for item in sublist])
		maxy = max([item for sublist in y_coord for item in sublist])
		miny = min([item for sublist in y_coord for item in sublist])

		x1 = 0 
		x2 = (int((maxx/km)/cell_size)*float(cell_size)) + cell_size
		y1 = 0
		y2 = (int((maxy/km)/cell_size)*float(cell_size)) + cell_size

		x_vals = (np.arange(x1, x2+cell_size, cell_size)*km)
		y_vals = (np.arange(y1, y2+cell_size, cell_size)*km)

		meshedx, meshedy = np.meshgrid(x_vals, y_vals)

		meshedxx = np.reshape(meshedx, (len(x_vals) * len(y_vals), 1), order = "C")
		meshedyy = np.reshape(meshedy, (len(x_vals) * len(y_vals), 1), order = "C")

		a = np.zeros((len(x_vals) * len(y_vals), 4))			# a will contain X and Y coordinates of squares
		b = np.zeros((len(x_vals) * len(y_vals), angle_divs))	# b will be 2d list of lengths by angle bin
		c = np.zeros((len(x_vals) * len(y_vals), angle_divs))	# c will be 2d list of number by angle bin
		d = np.zeros(len(a))	# d will be true number density ( to prevent double counting)

		
		a[:,0] = meshedxx[:,0]
		a[:,1] = meshedyy[:,0]
		temp_num = []	# temp_num will be populated with tuples of the form (square index, polyline id)

		for x, y, polyline_id in zip(x_coord, y_coord, range(len(x_coord))):
			temp_dens = []	# temp_dens is used to get density per angle_bin

			for i in np.arange(0, len(x)-1, 1):	# for each polyline:
				xi = int((x[i]/km)/cell_size)*float(cell_size) # starting square x
				xi_1 = int((x[i + 1]/km)/cell_size)*float(cell_size) # ending square y
				yi = int((y[i]/km)/cell_size)*float(cell_size) # starting square y
				yi_1 = int((y[i + 1]/km)/cell_size)*float(cell_size) # ending square y

				if (xi == xi_1) and (yi == yi_1):	# line stays in single square

					x_index = ((xi - x1)/cell_size)
					y_index = ((yi - y1)/cell_size)
					index = int((x_index + (y_index * len(x_vals)))) # index of square

					distance = np.hypot(x[i] - x[i + 1], y[i] - y[i + 1])	# distance

					frac_angle = angle_calc(x[i], y[i], x[i+1], y[i+1])	# angle
			
					# add to b in correct square, angle_index:
					angle_bin = int(frac_angle)/int(180/angle_divs)
					b[index][angle_bin] += distance	

					a[index][3] += distance # add to distance
					a[index][2] += 1		# add to row 3 (used to get rid of empty squares

					# add the square index and the angle bin to temp_dens:
					temp_dens.append((index, int(frac_angle)/int(180/angle_divs)))
					# add the square index and the polyline_id to temp_num:
					temp_num.append((index, polyline_id))

				else:	# if the polyline starts in a square and ends in a different square:
					# loop over each square which the line could inhabit
					for l in np.arange(min(xi, xi_1), max(xi, xi_1) + cell_size, cell_size):
						for j in np.arange(min(yi, yi_1), max(yi, yi_1) + cell_size, cell_size):

							x_index = (l - x1)/cell_size
							y_index = (j - y1)/cell_size

							index = int((x_index + (y_index * len(x_vals))))

							# use liangbarsky formula to work out distance in each square
							distance, xx1, yy1, xx2, yy2 = liangbarsky(l * 1000, (j + cell_size) * 1000, (l + cell_size) * 1000, j * 1000, x[i], y[i], x[i + 1], y[i + 1])
							if distance != 0:
								# if distance = 0, save relevant data
								frac_angle = angle_calc(xx1, yy1, xx2, yy2)

								angle_bin = int(frac_angle)/int(180/angle_divs)
								b[index][angle_bin] += distance	



								temp_dens.append((index, int(frac_angle)/int(180/angle_divs)))
								temp_num.append((index, polyline_id))

							a[index][3] += distance
							a[index][2] += 1
							
			# loop over temp_dens, extract data to c -> number by angle bin
			for h in set(temp_dens):
				c[h[0]][h[1]] += 1
			
		true_list = [] # list of bools to get rid of empty squares
		for row in a:
			if row[2] == 0 or row[3] == 0:
				true_list.append(False)
			else:
				true_list.append(True)

		true_list = np.array(true_list)

		temp_num =  list(set(temp_num))
		for number in temp_num:
			d[number[0]] += 1

		a = a[true_list]
		b = b[true_list]
		c = c[true_list]
		d = d[true_list]

		X, Y = a[:,0], a[:,1]

		self.X = X
		self.Y = Y
		self.N = c
		self.L = b / 1000.
		self.N_total = d
		self.__name__ = "FracAnalysisPoly"
		self.GDF = np.sum(self.N, axis = 0)/np.sum(self.N)
		self.Number_Anisotropy = np.max(self.N, axis = 1) / np.array([min([x for x in i if x > 0]) for i in self.N])

	def save_output(self, address_out):
		"""
			Creates shapefile of attributes

			Parameters
			-----------
			address_out : string
				location of shapefile to write out to

			-----------------------------------------------------------------
		"""
		#Set up shapefile writer and create empty fields
		w = shapefile.Writer(shapefile.POINT)
		w.autoBalance = 1 #ensures gemoetry and attributes match
		w.field('X','F',10,8)
		w.field('Y','F',10,8)

		for i in range(len(self.N.T)):
			header_string = "N" + str(i)
			w.field(header_string, "F", 10, 8)

		for i in range(len(self.N.T)):
			header_string = "L" + str(i)
			w.field(header_string, "F", 10, 8)

		w.field('N_tot','F',10,8)
		w.field("N_Anisotropy", "F", 10, 8)

		for index, value in enumerate(self.X):
			w.point(self.X[index], self.Y[index])

			zzz = [self.X[index], self.Y[index]]
			for j in self.N[index]:
				zzz.append(j)
			for j in self.L[index]:
				zzz.append(j/1000.)

			zzz.append(self.N_total[index])
			zzz.append(self.Number_Anisotropy[index])

			w.record(*zzz)

		#Save shapefile
		w.save(address_out)
		print "Saved output to " + address_out

################################################################################
################################################################################
################################################################################

class FracAnalysisPoint:
	"""
		Return analyzed shapefile in a grid.

		Parameters
		-----------
		address : string
			location of shapefile with only points

		cell_size : int, float
			size of grid square sides in km

		angle_divs : int
			number of angle divisions in 180 degrees

		Attributes
		-----------
		X, Y	: numpy.ndarray
			arrays of x and y-coordinates of grids with non-zero points

		N : numpy.ndarray
			array of number of counts in square, divided by angle bin

		N_total : numpy.ndarray
			array of true number of fractures per square (excludes double counting
			of fractures which are in >1 angle bin)

		Number_Anisotropy : numpy.ndarray
			array of number Anisotropy (1 value per square)

		Functions
		-----------
		save_output:
			saves attributes to a shapefile
		-----------------------------------------------------------------
	"""
	def __init__(self, address, cell_size, angle_divs):

		self.address = address
		self.cell_size = cell_size
		self.angle_divs = angle_divs

		sf = shapefile.Reader(address)

		data = list(sf.records())
		x_input = [x[0] for x in data]
		y_input = [y[1] for y in data]
		strike_input = [z[2] for z in data]

		if type(angle_divs) != int:
			print "\n The argument 'angle_divs' needs to be int\n"
			return

		km = 1000.

		# get bounds for plot
		maxx = max(x_input)
		minx = min(x_input)
		maxy = max(y_input)
		miny = min(y_input)

		x1 = 0 
		x2 = (int((maxx/km)/cell_size)*float(cell_size)) + cell_size
		y1 = 0
		y2 = (int((maxy/km)/cell_size)*float(cell_size)) + cell_size

		x_vals = (np.arange(x1, x2+cell_size, cell_size)*km)
		y_vals = (np.arange(y1, y2+cell_size, cell_size)*km)

		meshedx, meshedy = np.meshgrid(x_vals, y_vals)

		meshedxx = np.reshape(meshedx, (len(x_vals) * len(y_vals), 1), order = "C")
		meshedyy = np.reshape(meshedy, (len(x_vals) * len(y_vals), 1), order = "C")

		a = np.zeros((len(x_vals) * len(y_vals), 3))	# a will contain X and Y coordinates of squares
		c = np.zeros((len(x_vals) * len(y_vals), angle_divs))	# c will be 2d list of number by angle bin

		a[:,0] = meshedxx[:,0]
		a[:,1] = meshedyy[:,0]

		for x, y, strike in zip(x_input, y_input, strike_input):
			temp_dens = []	# temp_dens is used to get density per angle_bin

			xi = int((x/km)/cell_size)*float(cell_size) # starting square x
			yi = int((y/km)/cell_size)*float(cell_size) # starting square y

			x_index = ((xi - x1)/cell_size)
			y_index = ((yi - y1)/cell_size)
			index = int((x_index + (y_index * len(x_vals)))) # index of square

			a[index][2] += 1	# add to row 3 (used to get rid of empty squares
			try:
				c[index][int(strike)/int(180/angle_divs)] += 1
			except:
				c[index][(int(strike)-1)/int(180/angle_divs)] += 1

		true_list = [] # list of 
		for row in a:
			if row[2] == 0:
				true_list.append(False)
			else:
				true_list.append(True)

		true_list = np.array(true_list)


		a = a[true_list]
		c = c[true_list]
		X, Y = a[:,0], a[:,1]

		self.A = a
		self.X = X
		self.Y = Y
		self.N = c
		self.GDF = np.sum(self.N, axis = 0)/np.sum(self.N)
		self.N_total = np.max(self.N, axis = 1)
		self.__name__ = "FracAnalysisPoint"
		self.Number_Anisotropy = np.max(self.N, axis = 1) / np.array([min([x for x in i if x > 0]) for i in self.N])

	def save_output(self, address_out):

		"""
			Creates shapefile of attributes

			Parameters
			-----------
			address_out : string
				location of shapefile with only polylines

			-----------------------------------------------------------------
		"""
		#Set up shapefile writer and create empty fields
		w = shapefile.Writer(shapefile.POINT)
		w.autoBalance = 1 #ensures gemoetry and attributes match
		w.field('X','F',10,8)
		w.field('Y','F',10,8)


		for i in range(len(self.N.T)):
			header_string = "N" + str(i)
			w.field(header_string, "F", 10, 8)

		w.field('N_tot','F',10,8)
		w.field("N_Anisotropy", "F", 10, 8)

		for index, value in enumerate(self.X):
			w.point(self.X[index], self.Y[index])

			zzz = [self.X[index], self.Y[index]]
			for j in self.N[index]:
				zzz.append(j)

			zzz.append(self.N_total[index])
			zzz.append(self.Number_Anisotropy[index])

			w.record(*zzz)

		#Save shapefile
		w.save(address_out)
		print "Saved output to " + address_out

################################################################################
################################################################################
################################################################################

def FancyPlot(FracAnalyzed, Rose = True, Fractures = True, Patches = False, Circles = False, SquareNumbers = False, Title = "", FigureNumber = 1):

	"""Plots the point_number_density 

		Parameters
		-----------
		FracAnalyzed : list or single instance of FracAnalysisPoly/Point object
			Can be a list of different FracAnalysisPoint and/or FracAnalysisPoly
			objects. These objects must have the same cell_size and anlge_bins

		Rose : Boolean
			True by default
			Includes rose plots for each square

		Fractures : Boolean
			True by default
			Includes the 'raw' polylines and points in the plot

		Patches : "Number", "Length", "NumberAnisotropy" or False
			False by default
			Includes squares which indicate length density, number density, or 
			number anisotropy with colorbar

		Circles : Boolean
			False by default
			Includes circles around each Rose diagram, divided into sectors by
			angle bin. The circles indicate 25, 50, 75 and 100 % proportion.

		SquareNumbers : Boolean
			False by default
			Includes numbers in each square.

		Title: string
			Blank by default
			Title of the plot
		Returns
		--------
		A fancy plot.

		-----------------------------------------------------------------
	"""
	if type(FracAnalyzed) != list:
		FracAnalyzed = [FracAnalyzed]

	cell_size = FracAnalyzed[0].cell_size

	fig = plt.figure(FigureNumber)
	fig.clf()
	ax = fig.add_subplot(111)

	minx = min([min(temp.X) for temp in FracAnalyzed])
	maxx = max([max(temp.X) for temp in FracAnalyzed])

	miny = min([min(temp.Y) for temp in FracAnalyzed])
	maxy = max([max(temp.Y) for temp in FracAnalyzed])


	plt.xlim([minx, maxx + (cell_size * 1000)])
	plt.ylim([miny, maxy + (cell_size * 1000)])
	plt.grid(True)

	############################################################################
	# Drawing Fractures
	if Fractures == True:
		for classs in FracAnalyzed:

			if classs.__name__ == "FracAnalysisPoly": 
				sf = shapefile.Reader(classs.address)
				shapes = sf.shapes()

				for i in shapes[::]:
					a = np.array(i.points)
					plt.plot(a[:,0], a[:,1], zorder = 2)

			if classs.__name__ == "FracAnalysisPoint":
				sf = shapefile.Reader(classs.address)

				data = list(sf.records())
				x_input = [x[0] for x in data]
				y_input = [y[1] for y in data]

				plt.scatter(x_input, y_input, zorder = 2, color = "white", edgecolor = "black")		


	############################################################################
	# Drawing Patches Squares
	cmapp = plt.cm.jet

	#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
	# Number Density 
	if Patches == "Number":
		max_value = max([np.max(classs.N_total) for classs in FracAnalyzed])

		for classs in FracAnalyzed:
			X = classs.X
			Y = classs.Y
			ZZ = classs.N_total

			s = plt.scatter(X + (500 * cell_size), Y + (500 * cell_size), c = ZZ, cmap = cmapp, s = 1)

			for x, y, c, number in zip(X, Y, ZZ, range(len(ZZ))):
				ax.add_artist(plt.Rectangle(xy=(x, y), color = cmapp(c/max_value), width= cell_size * 1000, height=cell_size * 1000, alpha = 0.8))

		cbar = plt.colorbar(s)
		cbar.set_label("# of Fractures per\n   {} km squared".format(cell_size), fontsize = 16, rotation = 90)
		plt.clim(0, max_value)

	#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
	# Length Density
	elif Patches == "Length":
		max_value = max([np.max(np.sum(classs.L, axis = 1)) for classs in FracAnalyzed])

		for classs in FracAnalyzed:
			X = classs.X
			Y = classs.Y
			ZZ = np.sum(classs.L, axis = 1)

			s = plt.scatter(X + (500 * cell_size), Y + (500 * cell_size), c = ZZ, cmap = cmapp, s = 1)

			for x, y, c, number in zip(X, Y, ZZ, range(len(ZZ))):
				ax.add_artist(plt.Rectangle(xy=(x, y), color = cmapp(c/max_value), width= cell_size * 1000, height=cell_size * 1000, alpha = 0.8))

		cbar = plt.colorbar(s)
		cbar.set_label("Length of Fractures (km) per\n   {} km squared".format(cell_size), fontsize = 16, rotation = 90)
		plt.clim(0, max_value)

	#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#--#
	# Number Anisotropy 
	elif Patches == "NumberAnisotropy":
		max_value = max([np.max(classs.Number_Anisotropy) for classs in FracAnalyzed])

		for classs in FracAnalyzed:
			X = classs.X
			Y = classs.Y
			ZZ = classs.Number_Anisotropy

			s = plt.scatter(X + (500 * cell_size), Y + (500 * cell_size), c = ZZ, cmap = cmapp, s = 1)

			for x, y, c, number in zip(X, Y, ZZ, range(len(ZZ))):
				ax.add_artist(plt.Rectangle(xy=(x, y), color = cmapp(c/max_value), width= cell_size * 1000, height=cell_size * 1000, alpha = 0.8))

		cbar = plt.colorbar(s)
		cbar.set_label("Number Anisotropy per\n   {} km squared".format(cell_size), fontsize = 16, rotation = 90)
		plt.clim(0, max_value)

	############################################################################
	# Add numbers to the squares
	if SquareNumbers == True:
		for classs in FracAnalyzed:
			X = classs.X
			Y = classs.Y

			for x, y, number in zip(X, Y, range(len(X))):
				plt.text(x, y + ((0.7*cell_size)*1000), str(number), color = "white", zorder = 100000, fontsize = 10)

	############################################################################
	# plot rose diagram
	if Rose == True:
		for classs in FracAnalyzed:

			X = classs.X
			Y = classs.Y
			Z = classs.N
	
			patches = []
			for x, y, l in zip(X, Y, Z):
				L = np.sqrt(l/sum(l)) * cell_size * 500
				increment = 180./len(L)
				start = 0

				startx, starty = (x + (0.5 * cell_size)*1000), (y + (0.5 * cell_size)*1000)
				for i in L:
					end = start + increment
					wedgez = Wedge((startx, starty), i, 90 - end, 90 - start)
					patches.append(wedgez)

					wedgez = Wedge((startx, starty), i, (270 - end), (270 - start))
					patches.append(wedgez)
					start = start + increment

				if Circles == True:
					for classs in FracAnalyzed:
						radius = cell_size * 0.5
						# Circles around rose plot
						for i in (.25, .5, .75, 1):
							circle1 = plt.Circle((startx, starty), np.sqrt(i)*1000.0*radius, color = "black", fill = False, lw = 1.5, alpha = 0.25)
							fig.gca().add_artist(circle1)

						my_angle = 0.

						# Spokes in these circles
						for i in range(len(L)*2):
							plt.plot([startx, startx + (1000 * radius*np.sin(np.deg2rad(90 - 180-90-my_angle)))], [starty, starty + (1000. * radius*np.sin(np.deg2rad(90 - my_angle)))], color = "k", lw = 1.5, alpha = 0.25)
							my_angle += increment

			p = PatchCollection(patches, alpha=0.5, color = "black", zorder = 1000)
			ax.add_collection(p)
	plt.show()

################################################################################
################################################################################
################################################################################


def FancyPlotTotals(FracAnalyzed, Fractures = True, Circles = True, FigureNumber = 1):

	"""Plots the point_number_density 

		Parameters
		-----------
		FracAnalyzed : list or single instance of FracAnalysisPoly/Point object
			Can be a list of different FracAnalysisPoint and/or FracAnalysisPoly
			objects. These objects must have the same cell_size and anlge_bins

		Fractures : Boolean
			True by default
			Includes the 'raw' polylines and points in the plot

		Circles : Boolean
			False by default
			Includes circles around each Rose diagram, divided into sectors by
			angle bin. The circles indicate 25, 50, 75 and 100 % proportion.

		Title : string
			Blank by default
			Title of the plot

		FigureNumber : int
			Number of figure created

		Returns
		--------
		Another fancy plot.

		-----------------------------------------------------------------
	"""

	if type(FracAnalyzed) != list:
		FracAnalyzed = [FracAnalyzed]

	
	fig = plt.figure(FigureNumber)
	fig.clf()
	ax = fig.add_subplot(111)

	minx = min([min(temp.X) for temp in FracAnalyzed])
	maxx = max([max(temp.X) for temp in FracAnalyzed])

	miny = min([min(temp.Y) for temp in FracAnalyzed])
	maxy = max([max(temp.Y) for temp in FracAnalyzed])

	cell_size = FracAnalyzed[0].cell_size
	plt.xlim([minx, maxx + (cell_size * 1000)])
	plt.ylim([miny, maxy + (cell_size * 1000)])
	plt.grid(True)

	############################################################################
	# Drawing Fractures
	if Fractures == True:
		for classs in FracAnalyzed:

			if classs.__name__ == "FracAnalysisPoly": 
				sf = shapefile.Reader(classs.address)
				shapes = sf.shapes()

				for i in shapes[::]:
					a = np.array(i.points)
					plt.plot(a[:,0], a[:,1], zorder = 100)

			if classs.__name__ == "FracAnalysisPoint":
				sf = shapefile.Reader(classs.address)

				data = list(sf.records())
				x_input = [x[0] for x in data]
				y_input = [y[1] for y in data]

				plt.scatter(x_input, y_input, zorder = 100)			

	#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
	# Drawing the Rose plot
	for classs in FracAnalyzed:

		pointsx, pointsy = [], []
		sf = shapefile.Reader(classs.address)
		shapes = sf.shapes()

		for i in shapes[::]:
			a = np.array(i.points)
			pointsx.extend(a[:,0])
			pointsy.extend(a[:,1])

		X = classs.X
		X_max, X_min = max(X), min(X)
		X_width =  X_max - X_min

		Y = classs.Y
		Y_max, Y_min = max(Y), min(Y)
		Y_width = Y_max - Y_min

		cell_size = min([X_width, Y_width])/1000.0
		startx, starty = np.mean(pointsx), np.mean(pointsy)

		Z = np.sum(classs.N, axis = 0)

		radius = cell_size/2.
		a = Z/np.sum(Z)
		L = np.sqrt(a) * radius * 1000.
		increment = 180./len(L)
		start = 0
		patches = []

		for i in L:
			end = start + increment
			wedgez = Wedge((startx, starty), i, 90 - end, 90 - start)
			patches.append(wedgez)

			wedgez = Wedge((startx, starty), i, (270 - end), (270 - start))
			patches.append(wedgez)
			start = start + increment

		p = PatchCollection(patches, alpha=0.5, color = "black", zorder = 1000)
		ax.add_collection(p)

		#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
		# Add Circles
		if Circles == True:
			# Circles around rose plot
			for i in (.25, .5, .75, 1):
				circle1 = plt.Circle((startx, starty), np.sqrt(i)*1000.0*radius, color = "black", fill = False, lw = 1.5, alpha = 0.5)
				fig.gca().add_artist(circle1)

			my_angle = 0.

			# Spokes in these circles
			for i in range(len(L)*2):
				plt.plot([startx, startx + (1000 * radius*np.sin(np.deg2rad(90 - 180-90-my_angle)))], [starty, starty + (1000. * radius*np.sin(np.deg2rad(90 - my_angle)))], color = "k", lw = 1.5, alpha = 0.5)
				my_angle += increment

	plt.show()

################################################################################
