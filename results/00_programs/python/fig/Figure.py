## Matthew Dorsey
## @sunprancekid
## 19.02.2025
## class for containing information for figure properties which is consisten across all plots

##############
## PACKAGES ##
##############

import pandas as pd
import itertools # used for iterating over markers


################
## PARAMETERS ##
################

## constants, defaults for Label class
default_label_size = 12
minimum_label_size = 6

## constants, defaults for Figure class
default_title_label = None
default_subtitle_label = None
default_xaxis_label = None
default_yaxis_label = None
default_dpi = 200
minimum_dpi = 100
default_file_location = "./"
default_file_name = "figure.png"
marks = itertools.cycle(("D", "^", "v", "<", "o", "s", "p", "*")) # list of random markers


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
        if s is None or not isinstant(s, int):
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


## Figure class
class Figure (object):

    """ standard initialization routine for Figure object. """
    def __init__ (self):

        ## related to figure formatting, labelling
        self.set_title_label()
        self.set_subtitle_label()
        self.set_xaxis_label()
        self.reset_xaxis_limits()
        # self.set_xaxis_major_ticks()
        # self.set_xaxis_minor_ticks()
        self.set_yaxis_label()
        self.reset_yaxis_limits()
        # self.set_yaxis_major_ticks()
        # self.set_yaxis_minor_ticks()
        # self.set_cbar_label()
        # self.set_cscheme()
        # self.set_max_cbar()
        # self.set_min_cbar()
        # self.set_logscale()
        self.set_dpi()
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

    # initialize set of random set of markers that can be used for each unique ival in icol
    """ method generates a random set of markers than can be used with matplotlib. """
    def reset_markers(self):
        # create empty dictionary
        self.marker_dict = {}
        for i in self.get_unique_ivals():
            self.marker_dict.update({i: next(marks)}) # assign random marker to each ival

    # initialize list of labels that correspons to each unique ival in icol
    """ method initializes labels used to describe each unique ival in plots as that ival stored within that Figure dataframe. """
    def reset_labels(self):
        # create empty dictionary
        self.label_dict = {}
        for i in self.get_unique_ivals():
            self.label_dict.update({i: i}) # assign random marker to each ival

    # adjusts one marker in marker dictionary
    """ method changes one marker in the marker dictionary to a new marker type. the marker that is changed is the one that corresponds to the ival used as a key in the marker dictionary. """
    def set_marker(self, ival = None, marker = None):
        pass

    # return marker corresponding to ival
    """ method returns marker that correspons to ival in marker dictionary. """
    def get_marker (self, ival = None):
        return self.marker_dict[ival]

    # adjusts one label in label dictionary
    """ method changes one label in label dictionary to new string (not Label class). the label that is changed is the one that correspons to the ival used as a key in the label dictionary. """
    def set_label (self, ival = None, label = None):
        pass

    # returns one label in label dictionary
    """ method returns label that correspons to ival in label dictionary. """
    def get_label (self, ival = None):
        return self.label_dict[ival]

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

    # method that returns unique values for the isolation column
    """ returns list of all unique values contained within icol. """
    def get_unique_ivals (self):
        return self.df[self.icol].unique()

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

    ## TITLE ##

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

    # sets the xaxis label
    """ method for setting Figure x-axis label. """
    def set_xaxis_label(self, l = default_xaxis_label, s = None):
        self.xaxis_label = Label(l, s)

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


    ## YAXIS ##

    # sets the yaxis label
    """ method for setting Figure y-axis label. """
    def set_yaxis_label (self, l = default_yaxis_label, s = None):
        self.yaxis_label = Label(l, s)

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



############
## SCRIPT ##
############

# none
