## Matthew Dorsey
## @sunprancekid
## 14.03.2025
## methods used to handle problems relating to derivatives

## PACKAGES
# conda packages
import sys, os, math
import pandas as pd 
import numpy as np
import statistics
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy import interpolate # spline for curve fitting
from scipy import signal as sig # used for monotonic curve smoothing
# local packages
from analysis import monofit
from fig.Figure import Figure
from plot.plot import gen_plot
from plot.scatter_plot import gen_scatter

## PARAMETERS
# none

#####################################################
## METHODS USED FOR MONOTONIC SPLINE CURVE FITTING ##
#####################################################
#!/usr/bin/env python
"""https://stackoverflow.com/questions/56551114/fully-monotone-interpolation-in-python """
# see also
# https://en.wikipedia.org/wiki/Monotone-spline aka I-spline
# https://scikit-learn.org/stable/modules/isotonic.html
# denis 2 March 2020

def butter_filtfilt( x, Wn=0.5, axis=0 ):
    """ butter( 2, Wn ), filtfilt
        axis 0 each col, -1 each row
    """
    b, a = sig.butter( N=2, Wn=Wn )
    return sig.filtfilt( b, a, x, axis=axis, method="gust" )  # twice, forward backward

def ints( x ):
    return x.round().astype(int)

def minavmax( x ):
    return "min av max %.3g %.3g %.3g" % (
            x.min(), x.mean(), x.max() )

def pvec( x ):
    n = len(x) // 25 * 25
    return "%s \n%s \n" % (
        minavmax( x ),
        ints( x[ - n : ]) .reshape( -1, 25 ))

def monofit( y, Wn=0.1, verbose = False):
    """ monotone-increasing curve fit """
    y = np.asarray(y).squeeze()
    if verbose: print( "\n{ monofit: y %d %s  Wn %.3g " % (len(y), minavmax( y ), Wn ))
    ygrad = np.gradient( y )
    if verbose: print( "grad y:", pvec( ygrad ))

        # lowpass filter --
    gradsmooth = butter_filtfilt( ygrad, Wn=Wn )
    if verbose: print( "gradsmooth:", pvec( gradsmooth ))

    ge0 = np.fmax( gradsmooth, 0 )

    ymono = np.cumsum( ge0 )  # integrate, sensitive to first few
    ymono += (y - ymono).mean()

    err = y - ymono
    if verbose: print( "y - ymono:", pvec( err ))
    errstr = "average |y - monofit|: %.2g" % np.abs( err ).mean()
    if verbose: print( errstr )
    if verbose: print( "} \n" )

    return ymono, err, errstr

## METHODS
# equation used for linear fits
def linear_fit (x, m, b):
    return (x * m) + b

# equation used for linear fits with y-intercept at zero
def linear_fit_no_intercept (x, m):
    return linear_fit (x, m, 0.)

# converts linear scale to log scale
def lin2log(x, base = 10):
    return math.log(x) / math.log(base)

# converts log scale to linear scale
def log2lin (x, base = 10):
    return math.pow(base, x)

# method used to determine the slope of a line
# uses Savitzky-Golay filter to reduce noise in data
def filter_slope (x = None, y = None, log = False, order = 2, window = None, skip_int = 5):

    # convert data to logscale if requested
    x0 = [] # first order x-data
    y0 = [] # first order y-data
    if log:
        for i in range(len(x)):
            # if numbers are less than zero, convert to absolute scale
            if y[i] < 0.000000001:
                y[i] = abs(y[i])
            else:
                x0.append(lin2log(x[i], 10.))
                y0.append(lin2log(y[i], 10.))
    else:
        for i in range(len(x)):
            x0.append(x[i])
            y0.append(y[i])

    # based on the window, filter data sets
    dx_check = []
    x_plot = []
    y_plot = []
    y_plot_sg = []
    y_plot_data = []
    l_plot = []
    # initialize the sg filter averaging
    sg_x_cum = []
    sg_x_count = []
    sg_y_cum = []
    sg_y_count = []
    for i in range(len(x0) - skip_int * 2):
        sg_x_cum.append(0.)
        sg_x_count.append(0)
        sg_y_cum.append(0.)
        sg_y_count.append(0)

    for i in range(len(x0) - (window - 1)):
        # generate data set within window
        x0_fit = [] # zero order x-data corresponding to window
        y0_fit = [] # zero order y-data corresponding to window
        dx = 0. # used to determine average step size 
        # NOTE: SG filter assumes all data points are equally spaced ..
        #       so, dx should be the equal for all windows
        n = 0 
        x_avg = 0.

        # parse window data from set
        for j in range(window):
            x0_fit.append(x0[i + j])
            y0_fit.append(y0[i + j])
            x_avg += x0_fit[-1]
            if j > 0:
                # accumulate the window spacing
                dx = x0_fit[j] - x0_fit[j - 1]
                n += 1

        # apply SG filter
        x_avg = x_avg / n
        dx = dx / n
        dx_check.append(dx)
        y0_filter = sig.savgol_filter(y0_fit, window, order, delta = dx)

        # cumulate the filter for averaging
        for j in range(len(y0_filter)):
            if (j > skip_int - 1) and (j < (window - skip_int)):
                sg_x_cum[i + j - skip_int] += x0_fit[j]
                sg_x_count[i + j - skip_int] += 1
                sg_y_cum[i + j - skip_int] += y0_filter[j]
                sg_y_count[i + j - skip_int] += 1

    # convert everything back to linear scale
    for i in range(len(sg_y_cum)):
        sg_y_cum[i] = sg_y_cum[i] / sg_y_count[i] # average sg filter cumulation
        sg_x_cum[i] = sg_x_cum[i] / sg_x_count[i] 
        l_plot.append("SG filter")
        if log:
            x_plot.append(log2lin(sg_x_cum[i]))
            y_plot.append(log2lin(sg_y_cum[i]))
        else:
            x_plot.append(x0[i])
            y_plot.append(sg_y_cum[i])

    return x_plot, y_plot

    # for i in range(len(y0)):
    #     x_plot.append(x0[i])
    #     y_plot.append(log2lin(y0[i]))
    #     l_plot.append("Simulation Data")

    # # plot the fit against the real data, to double check
    # fig_check = Figure()
    # fig_check.load_data(d = pd.DataFrame.from_dict({'X': x_plot, 'Y': y_plot, 'L':l_plot}), xcol = 'X', ycol = 'Y', icol = 'L')
    # fig_check.reset_markers(['o', 'D'])
    # fig_check.set_logscale()
    # fig_check.set_yaxis_min(0.001)
    # fig_check.set_yaxis_max(1.1)
    # fig_check.set_xaxis_min(0.001)
    # fig_check.set_xaxis_max(10.)
    # fig_check.set_saveas('~/Desktop/', 'SG_filter')
    # fig_check.save_data()
    # gen_plot(fig = fig_check, markersize = 5, show = True, save = False)


# method used to determine the slope of a line
# options to smooth slope using a monotonic or spline function
def numerical_slope (x = None, y = None, log = False, monotonic = False, spline = False, monotonic_parameter = 0.25, average_int = 2):

    # check that the average int is greater than 2
    if average_int < 2:
        average_int = 2

    # convert data to a log scale
    x0 = []
    y0 = []
    if log:
        for i in range(len(x)):
            # remove numbers that are less than zero
            if y[i] < 0.000000001:
                continue
            else:
                x0.append(lin2log(x[i], 10.))
                y0.append(lin2log(y[i], 10.))
    else:
        for i in range(len(x)):
            x0.append(x[i])
            y0.append(y[i])

    # if smoothing was requested
    if monotonic:   
        # use monotonic function to calculate slope
        ymono, err, errstr = monofit (y0, Wn = monotonic_parameter)
    elif spline:
        pass

    # use smoothed function to calculate the slope
    xder = []
    yder = []
    for i in range(len(x0) - (average_int - 1)):
        x_fit = []
        y_fit = []
        for j in range(average_int):
            x_fit.append(x0[i + j])
            if monotonic:
                y_fit.append(ymono[i + j])
            else:
                y_fit.append(y0[i + j])
        # fit data to a line
        popt, pcov = curve_fit(f = linear_fit, xdata =  x_fit, ydata = y_fit)
        xder.append(statistics.mean(x_fit))
        yder.append(popt[0])

    # if log was called, return the xder to the linear scale
    if log:
        for i in range(len(xder)):
            xder[i] = log2lin(xder[i], 10.)

    return xder, yder

# method used to find the portion of an x-y relationship which corresponds to a certain slope
def find_regime (x = None, y = None, slope = None, log = False, tol = 0.1, spline = False, monotonic = False, average_int = 2, x_start = None, x_end = None, count_threshold = 10):

    ## CHECK ARGUMENTS
    # check that the slope has been specified as a float or an integer
    if slope is None or not isinstance(slope, float):
        if isinstance(slope, int):
            # if the slope is an integer, type cast it as a floating point number
            slope = float(slope)
        else:
            print(" smoothie :: find_regime :: ERROR :: 'slope' must be either integer or floating point number.")
            exit()

    ## CALCULATE THE NUMERICAL SLOPE
    # for i in range(len(x)):
    #     print(i, x[i], y[i])
    x0, dydx = numerical_slope(x = x, y = y, log = log, spline = spline, monotonic = monotonic, average_int = average_int)

    ## FIND THE REGION GRAPH WHICH CORRESPONDS TO THE TARGET SLOPE
    x_regime = [] # contains x data points corresponding to the desired regime
    m_regime = [] # contains slope data points corresponding to the desired regime

    # determine where to start the search
    if x_start is None or x_start > max(x):
        idx_start = 0
    else:
        for i in range(len(x0)):
            if x0[i] >= x_start:
                idx_start = i
                break

    # determing where to end the search
    if x_end is None or x_end > max(x):
        idx_end = len(x0)
    else:
        for i in range(idx_start, len(x0)):
            if x0[i] >= x_end:
                idx_end = i
                break

    # initialize other variables
    slope_exists = False
    lower_slope_lim = slope - tol
    upper_slope_lim = slope + tol
    while not slope_exists:
        # loop through each point, determine if the slope meets the criteria
        count = 0
        for i in range(idx_start, idx_end):
            # print(i, x0[i], dydx[i])
            if (dydx[i] > lower_slope_lim) and (dydx[i] < upper_slope_lim):
                # if the slope is within the range, cummulate the count
                count += 1
            else:
                # if the slope does not exist, reset the count
                count = 0
            # if the count surpases the threshold, the slope exists
            if count >= count_threshold:
                x_start = x0[i]
                x_start_idx = i - count_threshold + 1
                slope_exists = True
                break

        # if the algorithm gets to this point, it was unable to find a region corresponding to the slope above the threshold
        if not slope_exists:
            # reduce the threshold and restart
            if count_threshold > 2:
                count_threshold -= 1
            else:
                # the algorithm is unable to find any points within the threshold 
                # exit with none
                return None, None

    # once the start has been found, find the end of the desired regime
    count = count_threshold
    for i in range(x_start_idx, len(x0)):
        count += 1
        if (dydx[i] < lower_slope_lim) or (dydx[i] > upper_slope_lim):
            x_end = x0[i]
            x_end_index = i - count + 1
            break
        else:
            x_regime.append(x0[i])
            m_regime.append(dydx[i])
    
    return x_regime, m_regime


## ARGUMENTS
# none


## SCRIPT
# none