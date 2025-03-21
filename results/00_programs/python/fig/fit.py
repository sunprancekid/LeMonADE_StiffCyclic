## Matthew A. Dorsey
## @sunprancekid
## contains Line class, and methods to fit data to certain relationships

##############
## PACKAGES ##
##############
# from conda
import pandas as pd
import numpy as np
# local
from fig.Figure import Label

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
# method for fitting data to expotentional


# method for fitting data to line


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
		# self.reset_parameters()
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

	## TYPE ## 
	# used to describe the type of line

	# reset line type
	""" resets line type as none. """
	def reset_type (self):
		self.type = None

	# assign the line as having langevin type fit
	""" assign langevin type to line object. """
	def set_langevin_type(self):
		self.type = langevin_type

	""" returns string describing the type assigned to Line object. """
	def get_type(self):
		return self.type

	## PARAMETERS ## 

	# reset parameteres

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
		if not log:
			# create linear spacing 
			xvals = []
			for i in range(n):
				xvals.append((i / (n - 1)) * (xmax - xmin) + xmin)
		else:
			# create logarithimic spacing
			pass

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