## Matthew Dorsey
## @sunprancekid
## 14.03.2025
## methods used to handle problems relating to derivatives

## PACKAGES
import sys, os, math
import pandas as pd 
from analysis import monofit

## PARAMETERS
# none

#############################################
## USED FOR MONOTONIC SPLINE CURVE FITTING ##
#############################################
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
# converts linear scale to log scale
def lin2log(x, base):
    return math.log(x) / math.log(base)

# converts log scale to linear scale
def log2lin (x, base):
    return math.pow(base, x)

# calculate slope using a monotonic fit
def monotonic_slope (x = None, y = None, log = False, monotonic_parameter = 0.25):

	# convert data to a log scale
	if log:
		# TODO :: convert data to logscale if requested
		pass

	# use monotonic function to calculate slope
	ymono, err, errstr = monofit (y, Wn = monotonic_parameter)

	# use smoothed function to calculate the slope
	xder = []
	yder = []
	for i in range(len(x) - 1):
		# calculate slope
		xder.append((x[i + 1] + x[i]) / 2.)
		yder.append((ymono[i + 1] - ymono[i]) / (x[i + 1] - x[i]))

	return xder, yder


## ARGUMENTS
# none


## SCRIPT
# none