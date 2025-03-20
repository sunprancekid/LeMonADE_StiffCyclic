## Matthew A. Dorsey
## @sunprancekid
## contains Line class, and methods to fit data to certain relationships

##############
## PACKAGES ##
##############
# from conda
import pandas as pd
# local
from fig.Figure import Label

################
## PARAMETERS ##
################
# line types
logarithmic_type = "log"
linear_type = "linear"
# defaults
default_label_string = ""
default_label_size = 6
default_linestyle = '--'
default_linewidth = 2
default_linecolor = 'k'


#############
## METHODS ##
#############
# method for fitting data to expotentional


# method for fitting data to line


#############
## CLASSES ##
#############
# Line class
class Line(object):
	"""docstring for Line"""
	def __init__(self):
		self.set_label()
		self.reset_type()
		self.reset_parameters()
		self.line_style = default_linestyle # line style
		self.line_width = default_linewidth # width of line
		self.line_color = default_linecolor # key describing color of line

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

	""" returns string describing the type assigned to Line object. """
	def get_type(self):
		return self.type

	## PARAMETERS ## 

	# reset parameteres

###############
## ARGUMENTS ##
###############
# none


############
## SCRIPT ##
############
# none