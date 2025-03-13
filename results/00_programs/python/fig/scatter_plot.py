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
def gen_scatter(fig = None, edgecolor = default_edgecolor, markersize = default_markersize, legendloc = default_legendloc, show = True, save = True):

    # if no figure was provided ..
    if fig is None:
        # if a figure has not been specified, we have a problem: cannot generate a default
        # if figure and data are seperate objects, then maybe figure can be default while data would be mandatory
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
    plt.xlim(fig.get_xaxis_min(), fig.get_xaxis_max())
    plt.ylim(fig.get_yaxis_min(), fig.get_yaxis_max())
    if fig.get_title_label() is not None:
        plt.suptitle(fig.get_title_label().get_label(), fontsize = fig.get_title_label().get_size())
    if fig.get_subtitle_label() is not None:
        plt.title(fig.get_subtitle_label().get_label(), fontsize = fig.get_subtitle_label().get_size())
    plt.xlabel(fig.get_xaxis_label().get_label(), fontsize = fig.get_xaxis_label().get_size())
    plt.ylabel(fig.get_yaxis_label().get_label(), fontsize = fig.get_yaxis_label().get_size())

    # add logscale
    if fig.xaxis_is_logscale():
        plt.xscale(fig.get_xaxis_scale(), base = fig.get_xaxis_scale_base())
    if fig.yaxis_is_logscale():
        plt.yscale(fig.get_yaxis_scale(), base = fig.get_yaxis_scale_base())

    if save:
        plt.savefig(fig.get_saveas(), dpi = fig.get_dpi(), bbox_inches='tight')
    if show:
        plt.show()

    plt.close()


#############
## CLASSES ##
#############

# none

############
## SCRIPT ##
############

# none
