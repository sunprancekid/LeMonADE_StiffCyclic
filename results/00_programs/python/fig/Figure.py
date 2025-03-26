## Matthew Dorsey
## @sunprancekid
## 19.02.2025
## class for containing information for figure properties which is consisten across all plots

##############
## PACKAGES ##
##############

import pandas as pd
import numpy as np
import os # used to check path
import itertools # used for iterating over markers


################
## PARAMETERS ##
################

## constants, defaults for Label class
default_label_size = 12
minimum_label_size = 6

## constants, defaults used for Axis class
default_number_major_ticks = 3
default_number_minor_ticks = 3
default_major_format_string = "{:.2f}" # floating point number with two decimal places

## constants, defaults for Figure class
default_title_label = None
default_subtitle_label = None
default_xaxis_label = None
default_yaxis_label = None
default_dpi = 200
minimum_dpi = 100
default_file_location = "./"
default_file_name = "figure"
default_file_type = ".png"
default_markerset = ["D", "^", "v", "<", "o", "s", "p", "*"]
default_logscale_base = 10
minimum_logscale_base = 0.1
accepted_filetypes = [".png", ".tif"]
# defaults used for figures used in publications
pubdefault_dpi = 300
pubdefault_label_size = 18
# pubdefault_tick_size


#############
## METHODS ##
#############

# none


#############
## CLASSES ##
#############

## Label class
class Label (object):
    """ stand initialization routine for Label object. """
    def __init__ (self, l = None, s = None):
        self.set_label(l)
        self.set_size(s)

    """ output Label object as string. """
    def __str__ (self):

        out = ""
        if self.label is None:
            out += "Empty Label"
        else:
            out += self.get_label()
        out += " (font size: {})".format(self.get_size())

        return out


    ## GETTERS AND SETTERS ##

    # set label
    """ set label as string. """
    def set_label (self, l = None):
        if l is None or not isinstance(l, str):
            # if l is not an object or not a string
            self.label = None
        else:
            # otherwise l is a string
            self.label = l

    # return label
    """ return label as string. """
    def get_label (self):
        if self.label is None:
            return ""
        else:
            return self.label

    # set label fontsize
    """ get label size as interger no less than minimum. """
    def set_size (self, s = None):
        if s is None or not isinstance(s, int):
            # if s is not an object or not an integer
            self.size = default_label_size
        else:
            # if an integer was passed to the method, check the value is not less than the minimum
            if (s < minimum_label_size):
                self.size = minimum_label_size
            else:
                # otherwise assign the value passed to the method
                self.size = s

    # return label font size
    """ return label font size as integer. """
    def get_size (self):
        return self.size

## Axis class
class Axis (object):

    """ standard initialization for Axis object. """
    def __init__ (self):
        self.set_label()
        self.reset_limits()
        self.set_scale()
        self.set_major_ticks()
        self.set_minor_ticks()

## Figure class
class Figure (object):

    """ standard initialization routine for Figure object. """
    def __init__ (self):

        ## related to figure formatting, labelling
        self.set_title_label()
        self.set_subtitle_label()
        self.reset_xaxis_label()
        self.reset_xaxis_limits()
        # self.set_xaxis_major_ticks()
        # self.set_xaxis_minor_ticks()
        self.reset_yaxis_label()
        self.reset_yaxis_limits()
        # self.set_yaxis_major_ticks()
        # self.set_yaxis_minor_ticks()
        # self.set_cbar_label()
        # self.set_cscheme()
        # self.set_max_cbar()
        # self.set_min_cbar()
        self.set_linscale()
        self.set_dpi()
        self.set_saveas()
        # self.set_file_name() # TODO check file type
        # self.set_file_location() # TODO check that the location exists

        ## related to data and specification
        self.reset_data()

    """ output Figure object as string."""
    def __str__ (self):
        out = ""

        # add title string
        out += "\nTitle: {}".format(self.get_title_label_str())

        # add subtitle string
        out += "\nSub-Title: {}".format(self.get_subtitle_label_str())

        # add xaxis label string
        out += "\n\nX-Axis"
        out += "\nLabel: {}".format(self.get_xaxis_label_str())
        out += "\nMinimum Axis Limit: {}".format(self.get_xaxis_min())
        out += "\nMaximum Axis Limit: {}".format(self.get_xaxis_max())

        # add yaxis label string
        out += "\n\nY-Axis"
        out += "\nLabel: {}".format(self.get_yaxis_label_str())
        out += "\nMinimum Axis Limit: {}".format(self.get_yaxis_min())
        out += "\nMaximum Axis Limit: {}".format(self.get_yaxis_max())

        # add figure dpi
        out += "\nDPI: {}".format(self.dpi)

        ## TODO add data to string output

        # return to user
        return out

    """ adjust figure limits according to data passed to method. """
    def adjust_limits(self, df = None, xcol = None, ycol = None, ccol = None):

        # check for df
        if df is None or not isinstance(df, pd.DataFrame):
            # if df is not data frame, exit without adjusting axis limits
            return

        # set xaxis limits
        self.set_xaxis_limits(l = df[xcol].to_list())

        # set yaxis limits
        self.set_yaxis_limits(l = df[ycol].to_list())

        # set cbar limits

    """ method sets fontsizes for labels, ticks, etc. to presets for figures in publications. """
    def set_publication(self):
        self.filename = self.filename + "_pub"
        # adjust axis labels size
        self.set_xaxis_label(s = pubdefault_label_size)
        self.set_yaxis_label(s = pubdefault_label_size)
        # self.set_caxis_label(s = pubdefault_label_size)
        # adjust tick size
        # adjust legend label size
        # adjust dpi
        self.set_dpi (d = pubdefault_dpi)

    """ assign location, name, and type of save file. """
    def set_saveas(self, savedir = default_file_location, filename = default_file_name, filetype = default_file_type):

        # check the save directory passed to the method
        if not os.path.exists(savedir):
            # if the directory does not exist, make it
            os.makedirs(savedir)
            self.savedir = savedir
        else:
            self.savedir = savedir

        # check the filename passed to the method
        if type(filename) is not str:
            self.filename = default_file_name
        else:
            self.filename = filename

        # check the filetype passed to the method
        if filetype not in accepted_filetypes:
            self.filetype = default_file_type
        else:
            self.filetype = filetype

    """ return path to location of saved file. """
    def get_saveas(self):
        saveas = self.savedir + self.filename + self.filetype
        return saveas

    """ save data frame data to directory and with name associated with object."""
    def save_data(self):
        self.df.to_csv(self.savedir + self.filename + ".csv", index = False)

    ## GETTERS AND SETTERS ##

    ## DATA ##

    # method used to initialize data stored withing figure object
    """ initializes data stored withing figure object. dataframe is removed, x, y, c, and i columns are reset. """
    def reset_data (self):
        self.df = None
        self.xcol = None
        self.ycol = None
        self.ccol = None
        self.icol = None
        self.icol_marker_dict = None
        self.icol_label_dict = None

    # method that loads data from datafram into figure object
    """ loads data from data frame into Figure object. xcol specifies the xaxis data, ycol specifies the yaxis data, ccol specifies the color column data, icol specifies the isolation column data. """
    def load_data (self, d = None, xcol = None, ycol = None, ccol = None, icol = None):

        if d is None:
            # df has not been specified, cannot load data
            return

        # reset data and load
        self.reset_data()
        df_dict = {}

        # add xcol if it is specified
        if xcol is not None:
            self.xcol = xcol
            df_dict.update({xcol: d[xcol].to_list()})

        # add ycol if it is specified
        if ycol is not None:
            self.ycol = ycol
            df_dict.update({ycol: d[ycol].to_list()})

        # add ccol if it is specified
        if ccol is not None:
            self.ccol = ccol
            df_dict.update({ccol: d[ccol].to_list()})

        # add icol if it is specified
        if icol is not None:
            self.icol = icol
            df_dict.update({icol: d[icol].to_list()})

        self.df = pd.DataFrame(df_dict)
        self.reset_markers()
        self.reset_labels()

    # initialize set of random set of markers that can be used for each unique ival in icol
    """ method generates a random set of markers than can be used with matplotlib. """
    def reset_markers(self, markerset = None):

        # create interable list
        if markerset is None or type(markerset) is not list:
            # if the marker set is not passed to the method, or is not a list
            markerset = default_markerset
        else:
            # that marker set passed to the method is a list
            # check that the marker set contains the correct number of markers
            # if it does not, add to the marker list from the default until it contains the appropriate number
            if len(markerset) < len(self.get_unique_ivals()):
                n_default = len(default_markerset)
                i = 0
                while True:
                    if default_markerset[i] not in markerset or (i >= len(default_markerset)):
                        markerset.append(default_markerset[i])
                    i+=1
                    if len(markerset) >= len(self.get_unique_ivals()):
                        break

        # create empty dictionary
        marks = itertools.cycle(markerset)
        self.marker_dict = {} # empty dictionary
        for i in self.get_unique_ivals():
            self.marker_dict.update({i: next(marks)}) # assign random marker to each ival

    # adjusts one marker in marker dictionary
    """ method changes one marker in the marker dictionary to a new marker type. the marker that is changed is the one that corresponds to the ival used as a key in the marker dictionary. """
    def set_marker(self, ival = None, marker = None):
        if ival in self.marker_dict:
            self.marker_dict[ival] = marker

    # return marker corresponding to ival
    """ method returns marker that correspons to ival in marker dictionary. """
    def get_marker (self, ival = None):
        if ival is None:
            return default_markerset[0]
        else:
            return self.marker_dict[ival]

    # initialize list of labels that correspons to each unique ival in icol
    """ method initializes labels used to describe each unique ival in plots as that ival stored within that Figure dataframe. """
    def reset_labels(self):
        # create empty dictionary
        self.label_dict = {}
        for i in self.get_unique_ivals():
            self.label_dict.update({i: i}) # assign random marker to each ival

    # adjusts one label in label dictionary
    """ method changes one label in label dictionary to new string (not Label class). the label that is changed is the one that correspons to the ival used as a key in the label dictionary. """
    def set_label (self, ival = None, label = None):
        if ival in self.label_dict:
            self.label_dict[ival] = label

    # returns one label in label dictionary
    """ method returns label that correspons to ival in label dictionary. """
    def get_label (self, ival = None):
        return self.label_dict[ival]

    # method that determines if isolation column has been specified within the dataframe
    """ returns boolean the determines if isolation column has been specified within dataframe. """
    def has_ivals(self):
        return self.icol is not None

    # method that returns unique values for the isolation column
    """ returns list of all unique values contained within icol. """
    def get_unique_ivals (self, rev = False):
        if self.icol is not None:
            # if an icol has been specified return all unique items
            l = self.df[self.icol].unique()
            if rev:
                return np.flip(l)
            else:
                return l
        else:
            # otherwise, if an icol has not been specified, return a list with empty string
            return [""]

    # method that returns list of xvalues
    """ returns list of xvals contained within xcol. if ival is specified, xvals returned are those which share the same ival in icol (if any). """
    def get_xval_list (self, ival = None):
        if ival is None:
            # return all xvals as a list
            return self.df[self.xcol].to_list()
        else:
            # return xvals which share the same ival (if any)
            tmpdf = self.df[self.df[self.icol] == ival]
            return tmpdf[self.xcol].to_list()

    # method that returns list of yvalues
    """ returns list of yvals contained within ycol. if ival is specified, yvals returned are those which share the same ival in icol (if any)."""
    def get_yval_list (self, ival = None):
        if ival is None:
            # return all yvals as a list
            return self.df[self.ycol].to_list()
        else:
            # return yvales which share the same ival (if any)
            tmpdf = self.df[self.df[self.icol] == ival]
            return tmpdf[self.ycol].to_list()

    ## TITLE ##

    # sets the title label
    """ method for setting Figure title label. """
    def set_title_label(self, l = default_title_label, s = None):
        self.title_label = Label (l, s)

    # gets the title label as Label object
    """ method for returning Figure title as Label object. """
    def get_title_label (self):
        return self.title_label

    # gets the title label as string
    """ method for returning the Figure title as string. """
    def get_title_label_str (self):
        return str(self.title_label)

    ## SUBTITLE ##

    # sets the subtitle label
    """ method for setting Figure subtitle label. """
    def set_subtitle_label (self, l = default_subtitle_label, s = None):
        self.subtitle_label = Label (l , s)

    # gets the subtitle label as Label object
    """ method for getting the Figure subtitle label as Label object. """
    def get_subtitle_label (self):
        return self.subtitle_label

    # get the subtitle label as string
    """ method for getting the Figure subtitle label as string. """
    def get_subtitle_label_str (self):
        return str(self.subtitle_label)

    ## XAXIS ##

    # resets the xaxis label
    """ method for initializing Figure x-axis label as empty label. """
    def reset_xaxis_label (self, l = default_xaxis_label, s = None):
        self.xaxis_label = Label(l, s)

    # adjust xaxis label properties
    def set_xaxis_label (self, l = None, s = None):
        if l is not None:
            self.xaxis_label.set_label(l)
        if s is not None:
            self.xaxis_label.set_size(s)

    # gets xaxis label as object
    """ returns x-axis label as Label object. """
    def get_xaxis_label (self):
        return self.xaxis_label

    # gets xaxis label as string
    """ returns x-axis label as string. """
    def get_xaxis_label_str (self):
        return str(self.xaxis_label)

    # sets the minimum and maximum limits for the xaxis
    """ method used to assign the xaxis minimum and maximum values at the same time. """
    def set_xaxis_limits (self, l = None, min_val = None, max_val = None):

        # check for a list
        if l is not None:
            # pull min and max values from the list if they were not passed to the method
            if min_val is None:
                min_val = min(l)
            if max_val is None:
                max_val = max(l)

        # assign the minimum and maximum values
        self.set_xaxis_min(min_val)
        self.set_xaxis_max(max_val)

    # resets the minimum and maximum limits for the xaxis to None (used for initialization)
    """ method used to initialize the xaxis limits to None. """
    def reset_xaxis_limits (self):
        self.xaxis_min = None
        self.xaxis_max = None

    # sets the minimum value for the xaxis limit
    """ method used the assign the xaxis minimum limit as double. """
    def set_xaxis_min (self, val = None):
        # check the type passed as the minimum
        if isinstance(val, float):
            # if the minimum value is a float, assign it
            self.xaxis_min = val
        elif isinstance(val, int):
            # if the number is an integer, change it to a floating point and assign it
            self.xaxis_min = float(val)

    # sets the maximum value for the xaxis limit
    """ method used to assign the xaxis maximum as double. """
    def set_xaxis_max (self, val = None):
        # check the type passed as the maximum
        if isinstance (val, float):
            # if the maximum value passed to the method is a float, assign it
            self.xaxis_max = val
        elif isinstance (val, int):
            # if the value passed to the method is an integer, change it to a float and assign it
            self.xaxis_max = float(val)

    # gets the minimum value for the xaxis limit
    """ method that returns the minimum value assigned to the xaxis limit. returns 'None' if unassigned. """
    def get_xaxis_min (self):
        return self.xaxis_min

    # gets the maximum value for the xaxis limit
    """ method that returns the maximum value assigned to the xaxis limit. returns 'None' if unassigned. """
    def get_xaxis_max (self):
        return self.xaxis_max

    # set the xaxis as either linear or logscale
    """ method sets the xaxis as either a linear or logscale. logscale base set the logscale base as default if not specified by user. """
    def set_xaxis_scale (self, linear = False, log = False, logscale_base = default_logscale_base):

        # check arguments passed to method
        if not linear and not log:
            # scale was not specified, exit method
            return
        if linear and log:
            # both scales were specified, exit method
            return

        # if linear scale was specified
        if linear:
            self.xaxis_scale = None # if no scale is specified, then it is linear

        # if the logscale was specified
        if log:
            # check the basis
            if logscale_base < 0:
                # if the logscale basis is negative, use the positive
                logscale_base = -logscale_base
            elif abs(logscale_base) < minimum_logscale_base:
                # if the logscale base is less than the minimum, assign the default
                logscale_base = default_logscale_base

            # assign the logscale base
            self.xaxis_scale = logscale_base # if the base has a value, it is log rather than linear

    # check if xaxis is logscale
    """ method returns boolean determining if the xaxis is logscale or not. """
    def xaxis_is_logscale(self):
        return (self.xaxis_scale is not None)

    # returns the scale assigned to xaxis
    """ method that returns scale used for xaxis as either 'log' or 'linear' """
    def get_xaxis_scale (self):
        if self.xaxis_scale is None:
            return "linear"
        else:
            return "log"

    # returns scale base assigned to xaxis
    """ method that returns log base assigned to xaxis. if linear scale, returns None. """
    def get_xaxis_scale_base (self):
        return self.xaxis_scale

    # method used to set the tick marks and tick labels for the major and minor xaxis
    """ method sets the tick marks used for the xaxis."""
    def set_xaxis_ticks(self, minval = None, maxval = None, nmajorticks = default_number_major_ticks, nminorticks = default_number_minor_ticks, format_string = default_major_format_string):

        # check arguments passed to method
        # minimum value used for x axis
        if minval is None and self.get_xaxis_min() is None:
            # if no minimum value was passed to the method
            # and a minimum value has not already been assigned
            # assign the minimum value in the xcol as the minimum value
            self.set_xaxis_min(self.df[self.xcol].min())
        elif minval is not None:
            # otherwise a minimum value was passed to the method
            # assign that value passed to the method as the minimum for the x axis
            # overwrite the old minimum (if there was one)
            self.set_xaxis_min(minval)

        # maximum value used for xaxis
        if maxval is None and self.get_xaxis_max() is None:
            # if no maximum value was passed to the method
            # and a minimum value has not already been assigned
            # assign the maximum value in the x column as the max value
            self.set_xaxis_max(self.df[self.xcol].max())
        elif maxval is not None:
            # other a maximum value was passed to the method
            # assign the value passed to the method as the maximum for the y axis
            # overwrite the old maximum, if one was assigned
            self.set_xaxis_max(maxval)

        if self.xaxis_is_logscale():
            # if the yaxis is logscale, determine the tick marks along a logscale
            pass
        else:
            # otherwise, determine the tick marks along a linearscale
            pass

    # assigns values to ticks used for

    ## YAXIS ##

    # initializes the yaxis label
    """ method for initializing Figure y-axis label. """
    def reset_yaxis_label (self, l = default_yaxis_label, s = None):
        self.yaxis_label = Label(l, s)

    # adjust yaxis label properties
    def set_yaxis_label (self, l = None, s = None):
        if l is not None:
            self.yaxis_label.set_label(l)
        if s is not None:
            self.yaxis_label.set_size(s)

    # gets the y-axis label as object
    """ returns the y-axis label as string. """
    def get_yaxis_label (self):
        return self.yaxis_label

    # gets the y-axis label as string
    """ returns the y-axis label as string. """
    def get_yaxis_label_str (self):
        return str(self.yaxis_label)

    # set the yaxis minimum and maximum values
    """ method that assigns minimum and maximum values to the yaxis limits. """
    def set_yaxis_limits (self, l = None, min_val = None, max_val = None):

        # check if a list was passed to the method
        if l is not None:
            # use the list to assign min and max values, if unassigned
            if min_val is None:
                min_val = min(l)
            if max_val is None:
                max_val = max(l)

        # assign the minimum and maximum values
        self.set_yaxis_min(min_val)
        self.set_yaxis_max(max_val)

    # reset the yaxis minimum and maximum limits
    """ method used to initialize and reset yaxis limits to 'None' type. """
    def reset_yaxis_limits(self):
        self.yaxis_min = None
        self.yaxis_max = None

    # set yaxis minimum limit
    """ method that assigns a double as the yaxis minimum limit. """
    def set_yaxis_min (self, val):
        # check the value type passed to the method
        if isinstance (val, float):
            # if the value is a float, assign it
            self.yaxis_min = val
        elif isinstance (val, int):
            # if the value is an integer, change to a float and assign it
            self.yaxis_min = float(val)

    # set yaxis maximum limit
    """ method that assigns a double as the yaxis maximum limit. """
    def set_yaxis_max (self, val):
        # check the value type passed to the method
        if isinstance (val, float):
            # if the value is a float, assign it
            self.yaxis_max = val
        elif isinstance (val, int):
            # if the value is an integer, change to a float and assign it
            self.yaxis_max = float(val)

    # get the yaxis minimum limit
    """ method that returns the yaxis minimum limit as double. returns 'None' if unassigned. """
    def get_yaxis_min (self):
        return self.yaxis_min

    # get the yaxis maximum limit
    """ method that returns the yaxis maximum limit as double. returns 'None' is unassigned. """
    def get_yaxis_max (self):
        return self.yaxis_max

    # set the yaxis as either linear or logscale
    """ method sets the yaxis as either a linear or logscale. logscale base set the logscale base as default if not specified by user. """
    def set_yaxis_scale (self, linear = False, log = False, logscale_base = default_logscale_base):

        # check arguments passed to method
        if not linear and not log:
            # scale was not specified, exit method
            return
        if linear and log:
            # both scales were specified, exit method
            return

        # if linear scale was specified
        if linear:
            self.yaxis_scale = None # if no scale is specified, then it is linear

        # if the logscale was specified
        if log:
            # check the basis
            if logscale_base < 0:
                # if the logscale basis is negative, use the positive
                logscale_base = -logscale_base
            elif abs(logscale_base) < minimum_logscale_base:
                # if the logscale base is less than the minimum, assign the default
                logscale_base = default_logscale_base

            # assign the logscale base
            self.yaxis_scale = logscale_base # if the base has a value, it is log rather than linear

    # check if yaxis is logscale
    """ method returns boolean determining if the yaxis is logscale or not. """
    def yaxis_is_logscale(self):
        return (self.yaxis_scale is not None)

    # returns the scale assigned to yaxis
    """ method that returns scale used for yaxis as either 'log' or 'linear' """
    def get_yaxis_scale (self):
        if self.yaxis_scale is None:
            return "linear"
        else:
            return "log"

    # returns scale base assigned to yaxis
    """ method that returns log base assigned to yaxis. if linear scale, returns None. """
    def get_yaxis_scale_base (self):
        return self.yaxis_scale

    # method used to set the tick marks and tick labels for the major and minor yaxis
    """ method sets the tick marks used for the yaxis."""
    def set_yaxis_ticks(self, minval = None, maxval = None, nmajorticks = default_number_major_ticks, nminorticks = default_number_minor_ticks, format_string = default_major_format_string):

        # check arguments passed to method
        # minimum value used for yaxis
        if minval is None and self.get_yaxis_min() is None:
            # if no minimum value was passed to the method
            # and a minimum value has not already been assigned
            # assign the value as the minimum for the y axis
            self.set_yaxis_min(self.df[self.ycol].min())
        elif minval is not None:
            # otherwise a minimum value was passed to the method
            # assign that value passed to the method as the minimum for the y axis
            # overwrite the old minimum (if there was one)
            self.set_yaxis_min(minval)

        # maximum value used for yaxis
        if maxval is None and self.get_yaxis_max() is None:
            # if no maximum value was passed to the method
            # and a minimum value has not already been assigned
            # assign the maximum value from the yaxis as the maximum value
            self.set_yaxis_max(self.df[self.ycol].max())
        elif maxval is not None:
            # other a maximum value was passed to the method
            # assign the value passed to the method as the maximum for the y axis
            # overwrite the old maximum, if one was assigned
            self.set_yaxis_max(maxval)

        if self.yaxis_is_logscale():
            # if the yaxis is logscale, determine the tick marks along a logscale
            pass
        else:
            # otherwise, determine the tick marks along a linearscale
            pass

    ## DPI ##

    # sets the figure dpi
    """ sets the dpi for the figure. """
    def set_dpi (self, d = None):
        # check that the value passed to the method is an integer greater than the minimum
        if d is None or not isinstance(d, int):
            self.dpi = default_dpi
        else:
            if d < minimum_dpi:
                self.dpi = minimum_dpi

    # gets the figure dpi
    """ return the figure dpi (dots per inch). """
    def get_dpi (self):
        return self.dpi

    ## LINEAR OR LOG SCALE ##

    # sets scale for all axis to be linear
    """ method sets the scale for x, y, and c axis to be linear. """
    def set_linscale (self):

        # set x scale to be linear
        self.set_xaxis_scale(linear = True)

        # set y scale to be linear
        self.set_yaxis_scale(linear = True)

        # set c scale to be linear
        # self.set_caxis_scale(linear = True)

    def set_logscale (self, base = default_logscale_base, logx = True, logy = True):

        if logx:
            # set the xaxis scale to be log
            self.set_xaxis_scale(log = True, logscale_base = base)
            self.set_xaxis_ticks()

        if logy:
            # set the yaxis scale to be log
            self.set_yaxis_scale(log = True, logscale_base = base)
            self.set_yaxis_ticks()

        # set the xaxis scale to be log
        # self.set_caxis_scale(log = True, logscale_base = base)

############
## SCRIPT ##
############

# none
