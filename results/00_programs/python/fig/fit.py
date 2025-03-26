## Matthew A. Dorsey
## @sunprancekid
## contains Line class, and methods to fit data to certain relationships

##############
## PACKAGES ##
##############
# from conda
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
# local
from fig.Figure import Label
from util.smoothie import lin2log, log2lin

################
## PARAMETERS ##
################
# line types
logarithmic_type = "log"
linear_type = "linear"
langevin_type = "langevin"
# defaults
default_label_string = ""
default_label_size = 6
default_linestyle = '--'
default_linewidth = 2
default_linecolor = 'k'
default_marker = ''
default_markersize = 2


#############
## METHODS ##
#############
# method for parsing x and y data for fitting from starting and ending x-vales
def parse_data(x, y, x_start = None, x_end = None, log = False):
	
	# initialize arrays
	x_fit = []
	y_fit = []

	# if x_start has been specified, determine while to start adding the data which needs to be fit
	if x_start is None or x_start > max(x):
		idx_x_start = 0
	else:
		for i in range(len(x)):
			if x[i] >= x_start:
				idx_x_start = i
				break

	# cumulate the data as needed
	if x_end is None or x_end > max(x):
		for i in range(idx_x_start, len(x)):
			print(i, x[i], y[i])
			x_fit.append(x[i])
			y_fit.append(y[i])
	else:
		for i in range(idx_x_start, len(x)):
			if x[i] <= x_end:
				x_fit.append(x[i])
				y_fit.append(y[i])
			else:
				break

	if log:
		for i in range(len(x_fit)):
			x_fit[i] = lin2log(x_fit[i])
			y_fit[i] = lin2log(y_fit[i])

	return x_fit, y_fit

# method for fitting data to expotentional
def power_fit(x, p):
	return math.pow(x, p)

# uses paraemters which describe a linear relationship
def linear_fit(x, m, b):
	return m * x + b

# linear fit without y-intercept
def linear_fit_no_intercept (x, m):
	return linear_fit (x, m, 0.)

# linear fit without slope
def linear_fit_no_slope (x, b):
	return linear_fit (x, 0., b)

# methods for fitting data to line
# returns Line object which contains the parameters for linear relationship for data passed to method
def fit_line (x, y, x_start = None, x_end = None, log = False):
	
	# find the data that corresponds to the starting and ending regime, if specified
	x_fit, y_fit = parse_data(x, y, x_start = x_start, x_end = x_end, log = log)
	popt, pcov = curve_fit(f = linear_fit, xdata =  x_fit, ydata = y_fit)
	
	# create line object which contains parameters
	l = Line()
	l.set_label(l = "Linear Fit", s = None)
	l.set_log(log)
	l.set_linear_type(opt = popt, cov = pcov)
	return l

# returns Line object which fits data to line with out y-intercept
def fit_line_no_intercept(x, y, x_start = None, x_end = None):

	# fit the data to a line without a y-intercept
	x_fit, y_fit = parse_data(x, y, x_start = x_start, x_end = x_end)
	popt, pcov = curve_fit(f = linear_fit_no_intercept, xdata =  x_fit, ydata = y_fit)

	# create line with parameters fit to linear relationship
	l = Line()
	l.set_label(l = "Linear Fit", s = None)
	l.set_linear_type(opt = [popt[0], 0.], cov = pcov)
	return l

# returns Line object which fits data to line without slope
def fit_line_no_slope(x, y, x_start = None, x_end = None):

	# fit the data to a line without a y-intercept
	x_fit, y_fit = parse_data(x, y, x_start = x_start, x_end = x_end)
	popt, pcov = curve_fit(f = linear_fit_no_slope, xdata =  x_fit, ydata = y_fit)

	# create line with parameters fit to linear relationship
	l = Line()
	l.set_label(l = "Linear Fit", s = None)
	l.set_linear_type(opt = [0., popt[0]], cov = pcov)
	return l

# method for fitting data to langevin function
def langevin_fit():
	# initialize line object, assign label and type
	l = Line()
	l.set_label(l = "Langevin Function", s = None)
	l.set_langevin_type()
	return l


#############
## CLASSES ##
#############
# Line class
class Line(object):
	"""docstring for Line"""
	def __init__(self):
		self.set_label()
		self.reset_type()
		self.set_parameters() 
		self.set_log()
		self.set_linestyle() # assign default linestyle
		self.set_marker() # assign the default marker
		self.set_markersize() # assign the default markersize
		self.set_linewidth() # assign the default line width
		self.set_linecolor() # assign the default line color

	## LABEL ## 

	# method that sets the label used for the string
	""" assign label to Line object. """
	def set_label(self, l = default_label_string, s = default_label_size):
		self.label = Label(l, s)

	# returns label assigned to Line object
	""" returns Label object describing the Line object. """
	def get_label(self):
		return self.label

	## LOG SCALE ##
	# sets logscale for fit
	def set_log(self, l = False):
		self.log = l

	# return logscale
	def get_log (self):
		return self.l

	## TYPE ## 
	# used to describe the type of line

	# reset line type
	""" resets line type as none. """
	def reset_type (self):
		self.type = None

	# assign the line as having langevin type fit
	""" assign langevin type to line object. """
	def set_langevin_type(self, opt = None, cov = None):
		self.type = langevin_type

	def set_linear_type(self, opt = None, cov = None):
		self.type = linear_type
		self.set_parameters(opt)

	""" returns string describing the type assigned to Line object. """
	def get_type(self):
		return self.type

	## PARAMETERS ## 
	# set parameters
	def set_parameters (self, opt = None):
		self.parameters = opt

	# get parameters
	def get_parameters (self):
		return self.parameters

	## LINESTYLE ## 

	# set the linestyle used for the line
	""" set the line style used for the line. if not linestyle is specified, default is used. """
	def set_linestyle(self, ls = default_linestyle):
		self.linestyle = ls

	# get the linestyle used for the line
	def get_linestyle(self):
		return self.linestyle

	## MARKER ## 

	# set the marker used for plotting the line
	""" set the marker used when plotting the line. if no marker is specified, the default marker is assigned to the line. """
	def set_marker (self, m = default_marker):
		self.marker = m

	# returns the marker used for the line
	""" return the marker assigned to the line object. """
	def get_marker(self):
		return self.marker

	## MARKERSIZE ## 

	# set the markersize used for plotting the line
	""" set the marker size used for plotting the line. if not markersize as specified, assign the default."""
	def set_markersize(self, ms = default_markersize):
		self.markersize = ms 

	# get the markersize used for plotting the line
	""" return the markersize used for plotting the line assigned to the object. """
	def get_markersize(self):
		return self.markersize

	## LINEWIDTH ## 

	# set the linestyle used for plotting the line
	""" assign the linestyle used when plotting the object. if a linestyle is not passed to the method, assign the default. """
	def set_linewidth(self, ls = default_linewidth):
		self.linewidth = default_linewidth

	# get the line style used for plotting the object
	""" return the linestyle used for plotting the object. """
	def get_linewidth(self):
		return self.linewidth

	## LINECOLOR ##

	# set the line color assigned to the object
	""" assign the line color used when plotting the object. if a linestyle is not assigned to the method, assign the default. """
	def set_linecolor(self, lc = default_linecolor):
		self.linecolor = lc

	# return the line color assigned to the object
	def get_linecolor(self):
		return self.linecolor


	## X AND Y VALS ##

	# returns list containing x values used for line
	def get_xval_list(self, lims = None, n = None, log = False):
		# unpack touple containing min and max
		xmin, xmax = lims
		if log:
			# convert to log scale if requested
			xmin = lin2log(xmin)
			xmax = lin2log(xmax)

		# create linear spacing 
		xvals = []
		for i in range(n):
			x = (i / (n - 1)) * (xmax - xmin) + xmin
			if log:
				# convert back to linear scale 
				xvals.append(log2lin(x))
			else:
				xvals.append(x)

		return xvals

	# return list containing y value list
	def get_yval_list(self, lims = None, n = None, log = False):

		# get xvals
		xvals = self.get_xval_list(lims = lims, n = n, log = log)

		# use xvals and line type to get yvals
		yvals = []
		if self.type == langevin_type:
			for i in range(len(xvals)):
				x = xvals[i]
				if not (x == 0.):
					yvals.append((np.cosh(x) / np.sinh(x)) - (1. / x))
				else:
					yvals.append(0.)
		elif self.type == linear_type:
			linparms = self.get_parameters()
			for i in range(len(xvals)):
				if linparms is not None:
					if self.log:
						yvals.append(log2lin(linear_fit(lin2log(xvals[i]), *linparms)))
					else:
						yvals.append(linear_fit(xvals[i], *linparms))
		else:
			print("Line :: get_yval_list :: ERROR :: line type \'{}\' does not exists.".format(self.type))
			exit()

		return yvals

###############
## ARGUMENTS ##
###############
# none


############
## SCRIPT ##
############
# none