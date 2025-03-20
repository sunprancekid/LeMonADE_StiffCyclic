## Matthew A. Dorsey
## @sunprancekid
## constaints methods for generating plots using figure class


##############
## PACKAGES ##
##############
# from conda
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
# local
from fig.Figure import Figure


################
## PARAMETERS ##
################
# defaults associated with scatter plot
default_plot_markersize = 10
default_plot_linewidth = 2
default_edgecolor = 'black'
default_colormap = 'coolwarm'
default_legendloc = 'best'


#############
## METHODS ##
#############
# generate plot
def gen_plot (fig = None, linewidth = default_plot_linewidth, markersize = default_plot_markersize, legendloc = default_legendloc, show = True, save = True):

    # if no figure was provided ..
    if fig is None:
        # if a figure has not been specified, we have a problem: cannot generate a default
        # if figure and data are seperate objects, then maybe figure can be default while data would be mandatory
        exit()
        
    ## TODO :: check figure 
    

    # plot scatter
    leg = [] # empty list used for legend
    if fig.has_ivals():
        # if the figure has unique isolated values
        for i in fig.get_unique_ivals():
            leg.append(mlines.Line2D([], [], color = edgecolor, marker = fig.get_marker(i), ls = '', label = fig.get_label(i)))
            plt.plot(fig.get_xval_list(i), fig.get_yval_list(i), s = markersize) # == colorval, label == label, cmap == cmap, #  
        # add the legend
        plt.legend(handles = leg, loc = legendloc) # TODO increase size of legend labels
    else:
        # otherwise the figure does not have isolated values, so just create one plot
        plt.plot(fig.get_xval_list(), fig.get_yval_list(), linewidth = linewidth, marker = fig.get_marker(), markersize = markersize)

    # add min and max, labels
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


###############
## ARGUMENTS ##
###############
# none


############
## SCRIPT ##
############
# none