# Matthew Dorsey
# @sunprancekid
# 19.02.2025
# program for generating scatter plots within a common Figure framework

##############
## PACKAGES ##
##############

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from fig.Figure import Figure

################
## PARAMETERS ##
################

# defaults associated with scatter plot
default_markersize = 42
default_edgecolor = 'black'
default_colormap = 'coolwarm'
default_legendloc = 'best'

#############
## METHODS ##
#############

# generate scatter plot
def gen_scatter(fig = None, edgecolor = default_edgecolor, markersize = default_markersize, legendloc = default_legendloc):

    # if no figure was provided ..
    if fig is None:
        # if a figure has not been specified, we have a problem
        # if figure and data are seperate objects, then maybe figure can be a deafult
        exit()

    ## CHECK FIGURE
    # figure must have x
    # figure must have y
    # if figure does not have i, then an ival must be assigned to the entire df
    # if figure does not have c, then an cval must be assigned to the entire df (default is already tab)

    ## START FIGURE
    plt.figure()

    # plot scatter
    leg = [] # empty list used for legend
    for i in fig.get_unique_ivals():
        leg.append(mlines.Line2D([], [], color = edgecolor, marker = fig.get_marker(i), ls = '', label = fig.get_label(i)))
        plt.scatter(fig.get_xval_list(i), fig.get_yval_list(i), marker = fig.get_marker(i), s = markersize, edgecolor = edgecolor) # == colorval, label == label, cmap == cmap

    # add min and max, labels
    plt.legend(handles = leg, loc = legendloc) # TODO increase size of legend labels
    plt.xlabel(fig.get_xaxis_label().get_label(), fontsize = fig.get_xaxis_label().get_size())
    plt.ylabel(fig.get_yaxis_label().get_label(), fontsize = fig.get_yaxis_label().get_size())

    plt.show()
    # fig, ax = plt.subplots()
    # ax.set_xlim(fig.get_yaxis_min(), fig.get_xaxis_min())
    # ax.set_ylim(fig.get_yaxis_min(), fig.get_yaxis_max())


#############
## CLASSES ##
#############

# none

############
## SCRIPT ##
############

# none
