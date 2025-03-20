## Matthew A. Dorsey
## @sunprancekid
## contains Line class, and methods to fit data to certain relationships

##############
## PACKAGES ##
##############
# from conda
import pandas as pd
# local


################
## PARAMETERS ##
################
# defaults
default_linestyle = '--'


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
		self.linear = False # boolean determining if the line is of a linear type
		self.log = False 	# boolean determining if the line is of a logarithmic type
		self.label = None	# string used to label the line
		self.parameters = None # parameters describing the x -> relationship of the line
		self.line_style = None # line style
		self.line_width = None # width of line
		self.line_color = None # key describing color of line

	# method that sets the label used for the string


###############
## ARGUMENTS ##
###############
# none


############
## SCRIPT ##
############
# none