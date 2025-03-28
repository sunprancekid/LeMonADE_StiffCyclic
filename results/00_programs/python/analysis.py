## PACKAGES
# from conda
import sys, os, math
import csv
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit
from scipy.stats import norm
from scipy import interpolate # spline for curve fitting
from scipy import signal as sig # used for monotonic curve smoothing
import itertools # used for iterating over markers
from decimal import Decimal # used for formating large numbers to scientific notation
# local
from fig.Figure import Figure
from plot.scatter_plot import gen_scatter

################
## PARAMETERS ##
################

# file that contains parameters for scaling simulations
scaling_parmcsv = '01_raw_data/scaling/scaling.csv'
# file that contains parameter for forceExtension simulations
forceExtension_parmcsv = '01_raw_data/forceExtension/forceExtension.csv'
# file that contains parameters for bendingPARM simulations
bendingPARM_parmcsv = '01_raw_data/bendingPARM/bendingPARM.csv'
# default image resolution
default_dpi = 200
# generic tolerance for curve fitting
TOL=.0001
# equilibrium bond vector distribution for real polymer chain without bending potential
equildict = {'vec' : [0, 1, 2, 3, 4, 5], 'height': [0.0632855, 0.206359, 0.215854, 0.0550903, 0.231382, 0.228029]}
# list of random markers
marks = itertools.cycle(("v", "<", "o", "s", "p", "*", "D"))
# minimum cut off for logscale
min_logscale = 0.00001
# minimum force for pincus regime
min_pinus_force = 0.1
# maxmimum force for pincus regime
max_pincus_force = 1.
# maximum amount of time that any simulation can equilibriate
max_equilibriate_mcs = 100000000
# maximum amount of time that any one simulation can run
max_simulation_mcs = 5000000000

#############
## METHODS ##
#############

# format number to scientific notation string
# n  :: floating number
# sd :: number of significant digits (default is 2)
def format_e(n, sd = 2):
    d = Decimal(str(n))
    a = '%E' % d
    return a.split('E')[0].rstrip('0').rstrip('.') + 'E' + a.split('E')[1]

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

#...............................................................................
def monofit( y, Wn=0.1, verbose = False ):
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


############################
## USED FOR EQUATION FITTING
############################

# converts linear scale to log scale
def lin2log(x, base):
    return math.log(x) / math.log(base)

# converts log scale to linear scale
def log2lin (x, base):
    return math.pow(base, x)

# model equation for power law fitting
def power_fit(x, A, B):
    return A * (x ** B)

# model equation for power law fitting in the hookean regime
def power_fit_hook(x, A):
    return power_fit(x, A, 1.)

# model equation for fitting data to a line
def linear_fit (x, m, b):
    return (x * m) + b

# model equation for fitting line with no decay
def no_decay_fit (x, A):
    return power_fit(x, A, 0.)

# model equation for exponential decay fitting
def exp_decay_fit(x, A, B):
    return A * (math.e ** (- x / B))

# model equation for logistic5 equation
def log5_fit (x, A, B, C, D, E, F):
    return A + ((B - A) / ((1. + (C / (x - F)) ** D) ** E))

# calculate normalized end to end distance according to the normalized persistence length
# X :: persistence length normalized by total polymer length
# Y :: end-to-end distance squared, divided by total polymer length squared
def WLC_lp2R (x, Lo = 1.):
    return Lo * math.sqrt((2 * (x / Lo)) - 2 * (math.pow((x / Lo), 2.)) * (1 - math.exp(-1. / (x / Lo))))

#################################################
## USED FOR PARSING SIMULATION RESULTS FROM FILES
#################################################

# clean file, remove any comments from file structure
def clean_file(filepath, commentchar = '#'):

    # read file into array containing each line
    with open(filepath, 'r', encoding='UTF-8') as file:
        lines = [line.rstrip() + "\n" for line in file]

    # remove any lines that contain the comment character,
    # or that are empty. Write everything else to the same file
    with open(filepath, "w") as file:
        for l in lines:
            if ((len(l) > 1) and ( not commentchar in l)):
                file.write(l)
            # else:
            #     print(l)

# parse data from end-to-end file, return average
# filepath :: path to file containing data
# avgcol :: column header that contains relevant data
def parse_data(filepath, avgcol = None, avgrow = None, header = None, tab = False):

    # check that the file exists
    file = Path(filepath)
    if not file.exists():
        print(filepath + " does not exist!")
        return 0

    # load the file as a data frame
    if tab:
        df = pd.read_csv(filepath, sep = '\t', header = header)
    else:
        df = pd.read_csv(filepath, header = header)

    if avgrow is not None:
        rowlist = []
        row = df.iloc[avgrow]
        for i in range(1,len(row)):
            rowlist.append(row[i])
        return rowlist
    elif avgcol is not None:
        return df[avgcol].to_list()
    else:
        return []

# method for getting simulation results from files
def parse_results(parms = None, dir = None, simfile = None, col = None, row = None, title = None, bootstrapping = False, M1 = False, M2 = False, var = False, tabsep = False, plot = False):

    # make sure required information has been provided
    if title is None:
        pass
    if simfile is None:
        pass
    if parms is None:
        pass
    if col is None:
        pass

    # initialize arrays
    M1_arr = []
    M2_arr = []
    var_arr = []

    # loop through all simulation directories
    for i, r in parms.iterrows():
        if not os.path.exists(dir + r['path'] + simfile) or os.stat(dir + r['path'] + simfile).st_size < 1:
            # if the path does not exists or the file is empty, add NaN to each array
            M1_arr.append(np.nan)
            M2_arr.append(np.nan)
            var_arr.append(np.nan)
            # skip to next entry
            continue
        clean_file(dir + r['path'] + simfile)
        if os.stat(dir + r['path'] + simfile).st_size < 1:
            # the file is empty after cleaning, add NaN to each array
            M1_arr.append(np.nan)
            M2_arr.append(np.nan)
            var_arr.append(np.nan)
            # skip to next entry
            continue
        print(dir + r['path'] + simfile)
        if col is not None:
            vals = parse_data(dir + r['path'] + simfile, avgcol = col, tab = tabsep)
        elif row is not None:
            vals = parse_data(dir + r['path'] + simfile, avgrow = row, tab = tabsep)
        else:
            print("parse_results :: must specify either column or row for " + dir + r['path'] + simfile)
            exit()
        if title == 'ROG':
            vals = [math.sqrt(i) for i in vals]
        # calculate the first moment if, specified
        avg = 0
        nvals = len(vals)
        if bootstrapping and (nvals > 50):
            # fit data to normal distribution, determine average
            avg, std = norm.fit(vals)
            # avg = mu
        else:
            # calculate using the formula
            for i in vals:
                avg += i
            avg = avg / nvals

        # append average if requested
        if M1:
            M1_arr.append(avg)

        # calculate the second moment of average, if requested
        if M2 and (nvals > 2):
            avg2 = 0.
            for i in vals:
                avg2 += pow(i,2)
            avg2 = avg2 / len(vals)
            M2_arr.append(math.sqrt(avg2))

        # calculate the variance, if requested
        if var and (nvals > 2):
            if not bootstrapping:
                # if boot strapping was not used calculate using formula
                std = 0
                for i in vals:
                    std += pow(i - avg, 2)
                std = std / (len(vals) - 1)
            var_arr.append(math.sqrt(std)) # variance is sqrt of std

        # plot distribution if requested
        if plot:
            n_bins = 50
            fig, ax = plt.subplots()
            ax.hist(vals, bins=n_bins, density=True, label = "Simulation Data")
            xmin, xmax = plt.xlim()
            if bootstrapping:
                x = np.linspace(xmin, xmax, 100)
                p = norm.pdf(x, avg, std)
                plt.plot(x, p, 'k', linewidth=2, label = "Normal Distribution")
            plt.plot([], [], ' ', label = "$\mu$ = %.3f,  $\sigma$ = %.3f" % (avg, std))
            plt.title("Distribution of " + title + " Values ($n_{{bins}}$ = {:d})".format(n_bins))
            plt.legend()
            plt.xlabel("Property Value")
            plt.ylabel("Distribution")
            plt.savefig(dir + r['path'] + "propdist_" + title + ".png", dpi = 100)
            plt.close()


    # add to parameters data frame and return to user
    if M1: parms[title + "_M1"] = M1_arr
    if M2: parms[title + "_M2"] = M2_arr
    if var: parms[title + "_var"] = var_arr
    return parms

# method that calculates bond vector distribution error, generates diagram
def check_bvd (parms = None, dir = None, dpi = None, show = False, plot = False):

    # check that mandatory information was passed to method
    if parms is None:
        print("ERROR :: check_bvd :: must provide dataframe to method.")
        exit()
    if dir is None:
        print("ERROR :: check_bvd :: must provide path to main directory contain simulation files to method.")
        exit()

    # specify defaults if not provided by user
    if dpi is None:
        dpi = default_dpi

    # loop through each simulation, check for bvd files
    for i, r in parms.iterrows():
        # if the file does not exist, skip to the next file
        simfile = dir + r['path'] + 'BVD.dat'
        if not os.path.exists(simfile):
            continue
        # else, the file exists so make the graph
        pad = 0.1
        bar_width = 0.4
        # check if the file already has a header
        df = pd.read_csv(simfile, names = ['vec', 'height'])
        equildf = pd.DataFrame.from_dict(equildict)
        # TODO surpress warnings
        vec_diff = []
        sim_bar_xloc = []
        equil_bar_xloc = []
        mse = 0.
        for i in df.index:
            # calculate the difference between the equilibrium and actual BVD distribution
            val = df.at[i, 'height'] - equildf.at[i, 'height']
            mse += math.pow(val, 2.)
            vec_diff.append(val)
            # set bar center for plots
            if plot:
                sim_bar_xloc.append(df.at[i, 'vec'] - bar_width / 2)
                equil_bar_xloc.append(equildf.at[i, 'vec'] + bar_width / 2)
        # the vector difference to the dataframe and rewrite to the
        df['me'] = vec_diff
        df.to_csv(dir + r['path'] + "BVD.csv", header = False)
        if plot:
            fig, ax = plt.subplots()
            ax.bar(sim_bar_xloc, df['height'], label = "Excluded Volume + Bending Pot.", color = "tab:orange", edgecolor = "black", width = bar_width)
            ax.bar(equil_bar_xloc, equildf['height'], label = "Excluded Volume", color = "c", edgecolor = "black", width = bar_width)
            ax.set_xticks(np.linspace(0., 5., 6))
            ax.set_xticklabels(("$|2 0 0|$", "$|1 2 0|$", "$|2 1 1|$", "$|3 0 0|$", "$|2 2 1|$", "$|3 1 0|$"))
            ax.set_xlabel("Bond Vector")
            ax.set_ylabel("Probability Distribution")
            ax.legend()
            ax.set_title(f"Bond Vector Distribution ({r['id']})")
            plt.savefig(dir + r['path'] + 'BVD_dist.png', dpi = dpi)
            if show:
                plt.show()
            plt.close()

# method that calculates persistence length from bond-bond correlation
def check_bbc (parms = None, dir = None, dpi = None, show = False, plot = False, lp_decay = False, lp_theta = False):

    # check that mandatory information was passed to method
    if parms is None:
        print("ERROR :: check_bbc :: must provide dataframe to method.")
        exit()
    if dir is None:
        print("ERROR :: check_bbc :: must provide path to main directory contain simulation files to method.")
        exit()

    if (not show) and (not plot) and (not lp_decay) and (not lp_theta):
        # exit the routine, nothing to do here
        return

    # specify defaults if not provided by user
    if dpi is None:
        dpi = default_dpi
    if lp_decay:
        lpd = []
    if lp_theta:
        lpt = []

    # loop through each simulation, check for bvd files
    for i, r in parms.iterrows():
        # if the file does not exist, skip to the next file
        simfile = dir + r['path'] + 'BBC.dat'
        if not os.path.exists(simfile):
            continue

        # parse results from local file
        df = pd.read_csv(simfile, names = ['s', 'corr'])

        # calculate the persistence length according to exponential decay
        x = []
        y = []
        for j in range(1, 6):
            val = df.iloc[j]['corr']
            x.append(j)
            y.append(val)
        popt, pcov = curve_fit(exp_decay_fit, x, y)
        d = popt[1]
        if lp_decay:
            lpd.append(d)

        # calculate the persistance length according to local angle
        t = - 1. / math.log(df.iloc[1]['corr'])
        if lp_theta:
            lpt.append(t)

        # plot or show the bbc against the bond distance of seperation
        if plot or show:
            x = df['s'].tolist()
            y = df['corr'].tolist()
            x_fit = [ ((max(x) - x[1]) / (100 - 1)) * i + x[1] for i in range(100)]
            y_fit = [ exp_decay_fit(i, *popt) for i in x_fit]
            fig, ax = plt.subplots()
            ax.scatter(x[1:], y[1:], label = "Bond-Bond Correlation $g(s)$", marker = next(marks), s = 20)
            ax.plot(x_fit, y_fit, 'k--', label = "Expotential Decay Fit $e^{\\ell_{{p, decay}}}$")
            ax.plot([], [], ' ', label = "$\\ell_{{p, decay}}$" + " = {:.2f}".format(d))
            ax.plot([], [], ' ', label = "$\\ell_{p, \\theta}$" + " = {:.2f}".format(t))
            fig.suptitle("Bond-Bond Correlation and Persistence Length Calculations")
            ax.set_title("({:s})".format(r['id']))
            ax.set_ylabel("Bond-Bond Correlation ($\\langle cos \\theta \\rangle$)")
            ax.set_xlabel("Bond-Bond Distance ($s$)")
            # format axis
            ax.legend()
            # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), fancybox=True, shadow=True, ncol=4)
            ax.set_xlim(0, max(x))
            ax.set_ylim(.01, 1.)
            ax.set_yscale('log', base = math.e)
            # adjust y-axis labels to expotential log scale
            def ticks(y, pos):
                n = "{:.0f}".format(np.log(y))
                return r'$e^{:s}$'.format(n)
            ax = plt.gca()
            ax.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))
            if plot:
                fig.savefig(dir + r['path'] + 'BBC.png', dpi = dpi, bbox_inches='tight')
            if show:
                plt.show()
            plt.close()

        # save persistance length calculations as dataframe to csv
        dfsave = pd.DataFrame.from_dict({'lp_t': [t], 'lp_d': [d]})
        dfsave.to_csv(dir + r['path'] + 'BBC.csv', header = False)

# method that calculates the slope of the scattering factor
def check_skq (parms = None, dir = None, dpi = None, show = False, plot = False, porod = False, real = False, monotonic = False, monotonic_parameter = 0.25, spline = False , tol_mk = 0.2, debug = False, plot_fit = False):

    # check that mandatory information was passed to method
    if parms is None:
        print("ERROR :: check_sqk :: must provide dataframe to method.")
        exit()
    if dir is None:
        print("ERROR :: check_sqk :: must provide path to main directory contain simulation files to method.")
        exit()

    # specify defaults if not provided by user
    if dpi is None:
        dpi = default_dpi

    # loop through each simulation, check for sqk files
    for i, r in parms.iterrows():
        # if the file does not exist, skip to the next file
        simfile = dir + r['path'] + 'SKQ.dat'
        if not os.path.exists(simfile):
            continue

        if debug: print(simfile)

        # parse results from local file
        clean_file(simfile)
        df = pd.read_csv(simfile, sep = ' ', names = ['q', 'S', 'n'])

        # determine slope
        x = df['q'].tolist()
        y = df['S'].tolist()
        x_fit = []
        y_fit = []
        for i in range(len(x)):
            if x[i] > 0.1 and x[i] < 1.:
                x_fit.append(x[i])
                y_fit.append(y[i])
        # fit data with in range to linear slope
        popt, pcov = curve_fit(power_fit, x_fit, y_fit)
        m = popt[1]
        y_fit = [ power_fit(i, *popt) for i in x_fit ]

        # plot
        if plot or show:
            fig, ax = plt.subplots()
            ax.plot(x, y, 'k', label = 'S(q)')
            ax.plot(x_fit, y_fit, 'r--', label = 'Power Law Fit $A \cdot x^B$')
            ax.plot([], [], ' ', label = "A = {:.2f}".format(popt[0]))
            ax.plot([], [], ' ', label = "B = {:.2f}".format(m))
            ax.set_yscale('log')
            ax.set_xscale('log')
            fig.suptitle("Scattering Function")
            ax.set_title("({:s})".format(r['id']))
            ax.set_ylabel("S(q)")
            ax.set_xlabel("|q|")
            ax.legend()
            if plot:
                fig.savefig(dir + r['path'] + 'SKQ.png', dpi = dpi, bbox_inches='tight')
            if show:
                plt.show()
            plt.close()

        # create porod-kratsky plot
        if porod:
            # boolean determining if the cross over regime exists
            xc_exist = True

            # create plot by scaling the y-axis by scaled wave number
            y_log = []
            x_log = []
            nu = 0.5
            if real:
                nu = 0.588
            for i in range(len(y)):
                y[i] = y[i] * (pow(x[i], 1. / nu))
                # convert data to logscale for slope calculation
                y_log.append(lin2log(y[i], 10.))
                x_log.append(lin2log(x[i], 10.))

            # fit function, if specified, calculate the slope of the curve
            xder = []
            yder = []
            if monotonic:
                # use monofit to smooth function according to a parameter (Wn)
                ymono, err, errstr = monofit(y_log, Wn = monotonic_parameter)
                # use the smoothed data to calculate the derivative
                for i in range(len(x_log) - 1):
                    # calculate, convert from log scale to linear scale (except the slope)
                    xder.append(log2lin((x_log[i + 1] + x_log[i]) / 2., 10.))
                    yder.append((ymono[i + 1] - ymono[i]) / (x_log[i + 1] - x_log[i]))
                    ymono[i] = log2lin(ymono[i], 10.)
            elif spline:
                # use spline to interpolate data
                tck = interpolate.splrep(x_log, y_log, s=0)
                xder = np.arange(min(x_log), max(x_log), (max(x_log) - min(x_log)) / 1000)
                yspl = interpolate.splev(xder, tck, der=0)
                yder = interpolate.splev(xder, tck, der=1)
                for i in range(len(xnew)):
                    # convert from log scale to linear scale (except the slope)
                    xder[i] = log2lin(xder[i], 10.)
                    yspl[i] = log2lin(yspl[i], 10.)
            else:
                # calclate the slope from the raw data
                for i in range(len(x_log) - 1):
                    yder.append((y_log[i + 1] - y_log[i]) / (x_log[i + 1] - x_log[i]))
                    xder.append(log2lin(((x_log[i + 1] + x_log[i]) / 2.), 10.))

            # identify regimes with 'm = 0'
            tol_m0 = tol_mk
            while xc_exist:
                # find the beginning of the regime with a slope of 0
                for i in range(len(xder)):
                    if yder[i] < (0. + tol_m0):
                        break
                idx_m0_start = i

                # next, find the end of the regime with slope of 0
                x_m0 = []
                y_m0 = []
                for i in range(idx_m0_start, len(xder)):
                    if yder[i] > (0. + tol_m0):
                        break
                    else:
                        # if the slope is still within the linear regime, add the data to the list
                        x_m0.append(x[i])
                        y_m0.append(y[i])
                idx_m0_end = i

                # if the length of the array is greater than 1, proceed with the line fitting
                if len(x_m0) > 1:
                    # fit a line to the data within the linear regime, determine slope
                    popt_m0, pcov = curve_fit(power_fit, x_m0, y_m0)
                    # create fit for the linear regime
                    nfit = 10
                    xfit_m0 = []
                    yfit_m0 = []
                    xfit_m0_start = x_m0[0] * (1. / 3.)
                    xfit_m0_end = x_m0[-1] * (3.)
                    for i in range(nfit):
                        xfit_m0.append(xfit_m0_start + (i * (xfit_m0_end - xfit_m0_start) / (nfit - 1)))
                        yfit_m0.append(power_fit(xfit_m0[i], popt_m0[0], popt_m0[1]))
                    break
                else:
                    # increment the tolerance and repeat
                    tol_m0 += .01
                    # if the tolerance for identifying the cross-over regime has
                    # surpassed the threshold, the cross over regime cannot be idenitified or does not exist
                    if tol_m0 > (1. - tol_mk):
                        xc_exist = False


            # identify regimes with 'm = 1'
            tol_m1 = tol_mk
            if xc_exist:
                # find the beginning of the regime with a slope of 1
                for i in range(idx_m0_end, len(xder)):
                    if yder[i] > (1. - tol_m1):
                        break
                idx_m1_start = i

                # next, find the end of the regime with slope of 0
                x_m1 = []
                y_m1 = []
                for i in range(idx_m1_start, len(xder)):
                    if yder[i] > (1. + tol_m1):
                        break
                    else:
                        # if the slope is still within the linear regime, add the data to the list
                        x_m1.append(x[i])
                        y_m1.append(y[i])
                idx_m1_end = i

                # if the length of the array is greater than 1, proceed with the line fitting
                if len(x_m1) > 1:
                    # fit a line to the data within the rod-like regime, determine slope
                    popt_m1, pcov = curve_fit(power_fit, x_m1, y_m1)
                    # create fit for the linear regime
                    xfit_m1 = []
                    yfit_m1 = []
                    xfit_m1_start = x_m1[0] * (1. / 3.)
                    xfit_m1_end = x_m1[-1] * (3.)
                    for i in range(nfit):
                        xfit_m1.append(xfit_m1_start + (i * (xfit_m1_end - xfit_m1_start) / (nfit - 1)))
                        yfit_m1.append(power_fit(xfit_m1[i], popt_m1[0], popt_m1[1]))
                else:
                    # otherwise increment the tolerance and repeat
                    tol_m1 += 0.01

            # if the cross over regime exists, identify it
            if xc_exist:
                # use the two lines to identify the cross over regime, where to two lines intersect
                b_m0 = power_fit(1., popt_m0[0], popt_m0[1])
                b_m1 = power_fit(1., popt_m1[0], popt_m1[1])
                x_x = math.pow(10., math.log(popt_m0[0] / popt_m1[0], 10.) / (popt_m1[1] - popt_m0[1]))
                y_x = power_fit(x_x, popt_m0[0], popt_m0[1])

            # plot
            if plot or show:
                # create figure axis
                fig, ax = plt.subplots()
                ax2 = ax.twinx()

                # if the cross over regime exists, plot it
                if xc_exist:
                    # plot scattering, scattering slope as well as fit if specified
                    ax.plot(x, y, 'k', label = '$I(q) \cdot q^{{\\frac{{1}}{{\\nu}}}}$')
                    if monotonic:
                        # plot monotonic fit
                        if plot_fit: ax.plot(x, ymono, 'r', label = "Monotonic Fit ($W_{{n}}$ = {:.02f})".format(monotonic_parameter))
                        ax2.plot(xder, yder, color = 'b', alpha = 0.8)
                        ax.plot([], [], color = 'b', alpha = 0.8, label = 'Slope of Fit')
                    elif spline:
                        if plot_fit: ax.plot(xder, yspl, 'r', label = 'Cubic Spline Fit')
                        ax2.plot(xder, yder, color = 'b', alpha = 0.8)
                        ax.plot([], [], color = 'b', alpha = 0.8, label = 'Slope of Fit')
                    else:
                        ax2.plot(xder, yder, color = 'b', alpha = 0.8)
                        ax.plot([], [], color = 'b', alpha = 0.8, label = 'Slope of Data')

                    # plot lines that identify the two regimes
                    # slope with 'm = 0'
                    ax2.plot([xder[idx_m0_start], max(x)], [popt_m0[1], popt_m0[1]], '--', color = 'tab:purple') # line corresponding to slope
                    ax2.plot([xder[idx_m0_start], xder[idx_m0_start]], [popt_m0[1], max(yder)], '--', color = 'tab:purple') # start of linear regime
                    ax2.plot([xder[idx_m0_end], xder[idx_m0_end]], [popt_m0[1], max(yder)], '--', color = 'tab:purple') # end of linear regime
                    ax.plot(xfit_m0, yfit_m0, color = 'tab:purple', label = "Linear Regime \n (m = {:.02f})".format(popt_m0[1])) # line fit to data within linear regime plus data label

                    # slope with 'm = 1'
                    ax2.plot([xder[idx_m1_start], max(x)], [popt_m1[1], popt_m1[1]], '--', color = 'tab:orange') # line corresponding to rod_like regime
                    ax2.plot([xder[idx_m1_start], xder[idx_m1_start]], [popt_m1[1], max(yder)], '--', color = 'tab:orange') # start of linear regime
                    ax2.plot([xder[idx_m1_end], xder[idx_m1_end]], [popt_m1[1], max(yder)], '--', color = 'tab:orange') # end of linear regime
                    ax.plot(xfit_m1, yfit_m1, color = 'tab:orange', label = "Rod-Like Regime \n (m = {:.02f})".format(popt_m1[1]))

                    # plot the cross over
                    ax.plot(x_x, y_x, 'ro', markerfacecolor='none') # plot location of intersection of two lines
                    ax.plot([x_x, x_x], [min(y) , max(y)], '--', color = 'r') # plot line corresponding to cross over

                # set scales, labels
                ax.set_yscale('log')
                ax.set_xscale('log')
                if min(yder) > 0.: ax2.set_ylim(0., None)
                fig.suptitle("Kratky-Porod Plot")
                ax.set_title("({:s})".format(r['id']))
                ax.set_ylabel("$I(q) \cdot q^{{\\frac{{1}}{{\\nu}}}}$")
                ax2.set_ylabel("Slope")
                ax.set_xlabel("|q|")
                ax.legend()
                if plot:
                    fig.savefig(dir + r['path'] + 'SKQ_KP.png', dpi = dpi, bbox_inches='tight')
                if show:
                    plt.show()
                plt.close()
            # add cross over regime to fit
            dfsave = pd.DataFrame.from_dict({'A': [popt[0]], 'B': [m], 'C': [1. / x_x]})
            dfsave.to_csv(dir + r['path'] + 'SKQ.csv', header = False)
        else:
            # save fit as dataframe to csv with out cross over regime
            dfsave = pd.DataFrame.from_dict({'A': [popt[0]], 'B': [m]})
            dfsave.to_csv(dir + r['path'] + 'SKQ.csv', header = False)

# calculate hysteresis (area in loop) for force amplitude correlation
def calc_hys (f, r, omega = 1., f_a = 1., norm = 1.):

    # initialize hysteresis value
    a = 0.

    # loop though each neighboring pair, calculate area under curve according to trapazoidal method
    for i in range(len(f)):
        # establish integers used to access arrays
        i_1 = i
        i_2 = i + 1
        if (i_2 >= len(f)): i_2 = 0
        # accumulate hysteresis integral
        a += norm * 0.5 * (r[i_2] + r[i_1]) * (f[i_2] - f[i_1])

    return -a

# method that collects the hysteresis area, plots hysteresis curve
def check_hys (parms = None, dir = None, dpi = None, show = False, plot = False, check = True, clear = False):

    # check that mandatory information was passed to method
    if parms is None:
        print("ERROR :: check_hys :: must provide dataframe to method.")
        exit()
    if dir is None:
        print("ERROR :: check_hys :: must provide path to main directory contain simulation files to method.")
        exit()

    # specify defaults if not provided by user
    if dpi is None:
        dpi = default_dpi

    # loop through each simulation, check for sqk files
    for i, r in parms.iterrows():
        # if the file does not exist, skip to the next file
        simfile = dir + r['path'] + 'HYS.dat'
        if not os.path.exists(simfile):
            # if the hysteresis file is not present, the simulatin could be for calculating the no force chain extension
            simfile = dir + r['path'] + 'RE2E.dat'
            if os.path.exists(simfile):
                # if the end to end file does exist, parse the end to end distance as the hysteresis
                clean_file(simfile)
                data_E2E = parse_data(simfile, avgcol = 4, tab = True)
                data_ROG = parse_data(simfile, avgcol = 4, tab = True)
                dfsave = pd.DataFrame.from_dict({'A': [np.mean(data_E2E)], 'B': [np.var(data_ROG)], 'C': [np.mean(data_ROG)]})
                dfsave.to_csv(dir + r['path'] + 'HYS.csv', header = False)
                continue
            else:
                # file doesn't exist, write NaN
                dfsave = pd.DataFrame.from_dict({'A': [np.NaN], 'B': [np.NaN], 'C': [np.NaN]})
                dfsave.to_csv(dir + r['path'] + 'HYS.csv', header = False)
                print(simfile + " does not exist.")
                continue
        # parse the hysteresis area (if nan, delete the file so that it can be rerun)
        df = pd.read_csv(simfile)
        n = len(df.index)
        f_e = [ ] # collects the expected force values associated with the hysteresis
        f_a = [ ] # colelcts the averaged force values associated with the hysteresis
        f_s = [ ] # collects the standard deviation of force values associated with hysteresis
        r_a = [ ] # collects the averaged end to end distance values associated with the hysteresis
        r_s = [ ] # collects the standard deviation of the end to end distance associated with hysteresis
        hasnan = False

        if 'f_exp' not in df:
            # old hys file that has not been updated yet
            continue
        print(r['path'])
        # collect the values from the dataframe
        for j in range(1, n):
            f_e.append(df.loc[j,'f_exp'])
            f_a.append(df.loc[j,'f_avg'])
            f_s.append(df.loc[j,'f_std'])
            r_a.append(df.loc[j,'E2E_avg'])
            r_s.append(df.loc[j,'E2E_std'])
            if math.isnan(f_a[-1]) or math.isnan(r_a[-1]):
                hasnan = True

        # if the file has nan, delete the file and continue
        if clear and hasnan:
            print(simfile + " has NaN.")
            os.remove(simfile)
            continue

        # use the average end-to-end distance and force to estimate the hysteresis
        a = calc_hys(f_a, r_a, omega = 2. * math.pi / r['T'], f_a = r['fA'], norm = 300)

        # write the hysteresis value and variance to a csv file
        dfsave = pd.DataFrame.from_dict({'A': [abs(df.loc[0,'E2E_avg'])], 'B': [df.loc[0,'E2E_std']], 'C': [a]})
        dfsave.to_csv(dir + r['path'] + 'HYS.csv', header = False)

        if check and (plot or show):
            # check the correlation between the expected and measured forces
            # plt.errorbar(f_e, f_a, yerr = f_s, fmt = '.', ms = 5, lolims = True, uplims = True)
            plt.scatter(f_e, f_a, s = 2, color = 'r', label = "Actual Correlation")
            plt.xlabel("Expected Force ($f$)")
            plt.ylabel("Measured Force ($f$)")
            if max(f_e) > max(f_a):
                maxy = max(f_e)
            else:
                maxy = max(f_a)
            if min(f_e) < min(f_a):
                miny = min(f_e)
            else:
                miny = min(f_a)
            # plt.plot([min(f_e), max(f_e)], [min(f_e), max(f_e)], '--', linewidth=2, color = 'k', label = "Expected Correlation")
            plt.ylim(miny*1.1, maxy*1.1)
            plt.legend()
            plt.suptitle("Correlaton between Measured and Expected Force")
            plt.title("({:s})".format(r['id']))
            if plot:
                plt.savefig(dir + r['path'] + 'HYS_check.png', dpi = dpi, bbox_inches='tight')
            if show:
                plt.show()
            plt.close()

        if plot or show:
            # plot end to end distance against averaged force
            # would be cool to circularly color code each data point
            # plt.errorbar(f_a, r_a, yerr = r_s, fmt = 'o', lolims = True, uplims = True)
            plt.plot(f_a, r_a)
            plt.xlabel("Force ($f$)")
            plt.ylabel("Extension ($R$)")
            plt.ylim(-1., 1.)
            plt.suptitle("Hysteresis Force Extension")
            plt.title("({:s})".format(r['id']))
            if plot:
                plt.savefig(dir + r['path'] + 'HYS2.png', dpi = dpi, bbox_inches='tight')
            if show:
                plt.show()
            plt.close()

# add column to data frame which converts hysteresis simulation period frequency
def period2freq (parms = None, period_col = None, freq_col = None, timescale = 1.):

    # check for mandatory parameters
    if parms is None:
        print("ERROR :: period2freq :: must pass parameters as dataframe to method.")
        exit()
    if period_col is None:
        print("ERROR :: period2freq :: must specify title of column containing period data as 'period_col'.")
        exit()

    # fill optional parameters
    if freq_col is None:
        freq_col = 'frequency'

    # loop through all parameters, get period, covert to frequency
    f = [] # empty list containing frequency
    for i, r in parms.iterrows():
        t = r[period_col]
        if t == 0.:
            f.append(0.)
        else:
            f.append(timescale * 2. * math.pi / t)

    parms[freq_col] = f
    return parms


#######################################
## USED FOR PLOTTING SIMULATION RESULTS
#######################################

# method for getting scaling simulation data (N vs. R)
def plot_scaling(parms, N_col = None, R_col = None, fit = False, Title = None, X_label = None, Y_label = None, data_label = None, dpi = None, logscale_x = False, logscale_y = False, x_min = None, y_min = None, x_max = None, y_max = None, saveas = None, error = False, color = None, plot_log5 = False, plot_slope = False, flip_axis = False):
    # make sure that the proper parameters were passed to the method
    if N_col is None:
        exit()
    if R_col is None:
        exit()
    # establish image parameters based on parameters passed to method
    if dpi is None:
        dpi = default_dpi
    # set up for figure
    fig, ax1 = plt.subplots()
    for c in range(len(R_col)):
        # get values, average duplicates
        x = []
        y = []
        v = []
        for u in parms[N_col].unique():
            x.append(u)
            y_mean = parms[parms[N_col] == u]
            if logscale_y:
                y.append(abs(y_mean[R_col[c] + '_M1'].mean()))
            else:
                y.append(y_mean[R_col[c] + '_M1'].mean())
            v.append(y_mean[R_col[c] + '_var'].mean() / 2)
        # plot values
        if error:
            # plot with error bars
            if flip_axis:
                plt.errorbar(y, x, xerr = v, label = data_label[c], fmt = 'o', color = color[c])
            else:
                plt.errorbar(x, y, yerr = v, label = data_label[c], fmt = 'o', color = color[c])
        else:
            # plot without error bars
            if flip_axis:
                ax1.plot(y, x, 'o', label = data_label[c], fillstyle = 'none', color = color[c])
            else:
                ax1.plot(x, y, 'o', label = data_label[c], fillstyle = 'none', color = color[c])
        # fit data to power low, if specified
        if fit:
            # fit data to power law equation, determine parameters
            popt, pcov = curve_fit(power_fit, x, y)
            x_fit = [ ((max(x) - min(x)) / (100 - 1)) * i + min(x) for i in range(100)]
            y_fit = [ power_fit(i, *popt) for i in x_fit]
            ax1.plot(x_fit, y_fit, '--', label = data_label[c] + " Fit", color = color[c])
            # plt.plot([], [], ' ', label="A = {:.3f}".format(popt[0]))
        if plot_log5 or plot_slope:
            # if plotting either the log5 curve or the slope of the log5 curve fit
            # convert that dataset to a logscale
            x_log = []
            y_log = []
            v_log = []
            for i in range(len(x)):
                if abs(y[i]) < TOL or math.isnan(y[i]):
                    continue
                x_log.append(math.log(x[i]) / math.log(10.))
                y_log.append(math.log(abs(y[i])) / math.log(10.))
                if v[i] > TOL:
                    v_log.append(v[i])
                else:
                    v_log.append(TOL)
            # generate initial parameters for log5 curve
            lb = [min(y_log), max(y_log) - (TOL / 100), -100, -100, -40, min(x_log) - 2 * TOL]
            ub = [max(y_log), max(y_log), 200, 100, 40, min(x_log) - TOL]
            parameterBounds = (lb, ub)
            init = [min(y_log), max(y_log), 4., 50., 0.01, min(x_log) - (3. * TOL / 2.)]
            # fit data to power log5 equation
            popt, pcov = curve_fit(f = log5_fit, xdata =  x_log, ydata = y_log, p0 = init, bounds = parameterBounds, maxfev = 100000000) # sigma = v_log,
            x_fit = [ ((max(x_log) - min(x_log)) / (100 - 1)) * i + min(x_log) for i in range(100)]
            y_fit = [ log5_fit(i, *popt) for i in x_fit ]
            # calculate the slope of the line
            if plot_slope:
                # calculate the of the line
                h = (max(x_fit) - min(x_fit)) / len(x_fit) # distance between points
                x_d1 = []
                y_d1 = []
                for i in range(len(y_fit) - 1):
                    x_d1.append(math.pow(10, (x_fit[i+1] + x_fit[i]) / 2.))
                    y_d1.append((y_fit[i+1] - y_fit[i]) / h)
            # convert everthing from the log scale, back the linear scale; plot
            if plot_log5:
                for i in range(len(x_fit)):
                    x_fit[i] = math.pow(10, x_fit[i])
                    y_fit[i] = math.pow(10, y_fit[i])
                ax1.plot(x_fit, y_fit, '--', label = "log 5 Curve", color = color[c], alpha = 0.8)
            if plot_slope:
                # plot slope already calculated
                ax2 = ax1.twinx()
                ax2.plot(x_d1, y_d1, '-', label = "First Derivative of log5", color = color[c], alpha = 0.8)
                ax1.plot([], [], '-', label = "First Derivative of log5", color = color[c], alpha = 0.8)
                ax2.set_ylabel("Change in " + Y_label)
                ax2.set_ylim(0, 3.)
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    if logscale_x:
        ax1.set_xscale('log')
    if logscale_y:
        ax1.set_yscale('log')
    ax1.legend(loc = 'upper left')
    if Title is not None:
        ax1.set_title(Title)
    if Y_label is not None:
        if flip_axis:
            ax1.set_xlabel(Y_label)
        else:
            ax1.set_ylabel(Y_label)
    if X_label is not None:
        if flip_axis:
            ax1.set_ylabel(X_label)
        else:
            ax1.set_xlabel(X_label)
    if saveas is not None:
        fig.savefig(saveas, dpi = dpi, bbox_inches='tight')
    plt.show()
    plt.close()

# method for plotting data from force extension simulations
def plot_force_extension(parms, X_col = None, Y_col = None, iso_col = None, isovals = None, isolabel = None, logscale_x = False, logscale_y = False, x_min = None, x_max = None, y_min = None, y_max = None, X_label = None, Y_label = None, Title = None, saveas = None, plot_log5 = False, plot_data = False, plot_slope = False, dpi = None, M2 = False, error = False, flip_axis = False, horizbar = False, vertbar = False, show = False, legend_loc = "best", norm = False, pincus = False, hookean = False, timescale = 1., hys = False, fa = None, diff = False):

    if X_col is None or Y_col is None or iso_col is None:
        print("plot_force_extension :: Must specify columns from data frame to plot!")
        exit()

    # establish image parameters based on parameters passed to method
    if dpi is None:
        dpi = default_dpi

    # booleans for documenting if the pincus regime and the hookean regime have been added to the figures
    added_pinus = False
    added_hookean = False
    added_hys_pos = False
    added_hys_neg = False
    A_max_neg = None
    A_max_pos = None
    y_min_check = None
    y_max_check = None

    # used for determing minimum and maximum axis
    y_min_find = None
    y_max_find = None
    x_max_find = None
    x_min_find = None

    if diff and hys:
        # initialize diffusion correlation for hysteresis
        x_diff = []
        y_diff = []

    # parse data from data frame, plot
    if isovals is None:
        isovals = parms[iso_col].unique()
    for iso in isovals:

        y_col_m = ""
        if M2:
            y_col_m = Y_col + '_M2'
        else:
            y_col_m = Y_col + '_M1'
        # for each unique isolated value
        iso_df = parms[parms[iso_col] == iso]
        # remove any data points that are negative, if log scale is being used
        rmv_row = [] # list of rows to remove from data frame
        for i, row in iso_df.iterrows():
            if ((logscale_x and (row[X_col] < min_logscale)) or (logscale_y and (row[y_col_m] < min_logscale))):
                rmv_row.append(i)
        if (len(rmv_row) > 0):
            iso_df.drop(rmv_row)

        x = iso_df[X_col].tolist()
        y = iso_df[y_col_m].tolist()
        # y_norm = iso_df['E2Etot_M1'].tolist()
        if norm and hys:
            # if the normalization is true
            # find the zero force extension value
            idx = x.index(0.)
            norm_val = y[idx] # this value corresponds to the equilibrium chain length
            # remove from the values from the list as they do not actually correspond to hystresis
            del x[idx]
            del y[idx]
            # identify the maximum frequency at which the hysteresis occurs
            timescale = 1. / x[y.index(max(y))]
            # use the maximum frequency to normalize the
            for i in range(len(x)):
                x[i] = x[i] * timescale
            if fa is not None:
                for i in range(len(y)):
                    y[i] = y[i] * pow(fa, 2.) / norm_val
            if diff:
                # if calculating diffusion determine bending strength and approximate scaling
                x_diff.append(iso)
                y_diff.append(timescale / pow(norm_val, 2.))
        # normalize the force and end-to-end distance if asked
        elif norm and not hys:
            # get the no-force chain extension
            idx = x.index(0.)
            norm_val = y[idx]
            # remove from the list
            del x[idx]
            del y[idx]
            # normalize each element by the no force chain extension
            for i in range(len(x)):
                x[i] = x[i] * norm_val
                y[i] = y[i] / norm_val
        if error:
            v = iso_df[Y_col + '_var'].tolist()
        if logscale_y is True:
            y = [abs(i) for i in y]
        # fit the log data to slope
        if plot_slope:
            # convert data to log scale
            x_log = []
            y_log = []
            v_log = []

            # for data points with regime
            for i in range(len(x)):
                if abs(y[i]) < TOL or math.isnan(y[i]):
                    continue
                # if withing the pincus regime
                if (x[i] > min_pinus_force) and (x[i] < max_pincus_force):
                    x_log.append(math.log(x[i]) / math.log(10.))
                    y_log.append(math.log(abs(y[i])) / math.log(10.))
                    if error:
                        if v[i] > TOL:
                            v_log.append(v[i])
                        else:
                            v_log.append(TOL)

            # fit to a line
            popt, pcov = curve_fit(f = linear_fit, xdata =  x_log, ydata = y_log)

            # fit data, convert back to logscale
            x_fit = [ (((max(x_log) + 0.5) - (min(x_log) - 1.)) / (100 - 1)) * i + (min(x_log) - 1.) for i in range(100)]
            y_fit = [ linear_fit(i, *popt) for i in x_fit]
            for i in range(len(x_fit)):
                x_fit[i] = math.pow(10, x_fit[i])
                y_fit[i] = math.pow(10, y_fit[i])

        # plot the data
        label = isolabel.format(iso)
        if norm and hys:
            label += " ($\\tau_{{R}} \\approx$ {:.01E})".format(Decimal(timescale))
        if (norm and not hys):
            label += " (X = {:.02f})".format(norm_val)
        if flip_axis:
            line = plt.plot(y, x, label = label, fillstyle = 'full', marker = next(marks), ls = ' ')
        else:
            line = plt.plot(x, y, label = label, fillstyle = 'full', marker = next(marks), ls = ' ')
            # check y range
            if y_max_check is None:
                y_max_check = max(y)
            elif y_max_check < max(y):
                y_max_check = max(y)

            if y_min_check is None:
                y_min_check = min(y)
            elif y_min_check > min(y):
                y_min_check = min(y)


        # plot the slope if requested
        if plot_slope:
            # plot the slope of the line fit to the simulation data
            plt.plot(x_fit, y_fit, '--', label = "m = {:.03f}".format(popt[0]), color = line[-1].get_color())

        # add pincus regime to graph
        if (norm and not hys) and pincus and (not added_pinus):
            # plot a line with a slope of 2/3
            # plot the line tangent to one point within the pincus regime
            # find one point within the pincus regime from the data
            for i in range(len(x)):
                if x[i] > 2.:
                    break
            # use the data point to identify the y-intercept of the line, on log scale
            x_log = lin2log(x[i], 10.)
            y_log = lin2log(y[i], 10.)
            A = math.pow(10., y_log - (2. / 3.) * x_log)
            # create that data points on a linear scale
            x_pincus = np.linspace(0.5, 1000, 20)
            y_pincus = [ power_fit(i, A, 2./3.) for i in x_pincus ]
            # plot
            # if norm:
            #     plt.plot(x_pincus, y_pincus, '--', color = 'k')
            # else:
            #     plt.plot(x_pincus, y_pincus, '--', color = line[-1].get_color())
            added_pinus = True

        if (norm and not hys) and hookean and (not added_hookean):
            # plot a line with a slope of 1
            # plot the line tangent to one point within the hookean regime
            # find one point within the hookea regime from the data
            for i in range(len(x)):
                if x[i] > 0.5:
                    break
            # use the data point to identify the y_intercept of the line, on log scale
            x_log = lin2log(x[i], 10.)
            y_log = lin2log(y[i], 10.)
            A = math.pow(10., y_log - (1.) * x_log)
            # create data points on a linear scale
            x_hookean = np.linspace(0.0001, 5., 20)
            y_hookean = [ power_fit_hook(i, A) for i in x_hookean ]
            # plot
            # if norm:
            #     plt.plot(x_hookean, y_hookean, ls = ':', color = 'k')
            # else:
            #     plt.plot(x_hookean, y_hookean, ls = ':', color = line[-1].get_color)
            added_hookean = True

        # add lines with slopes of one and negative point five
        if norm and hys:
            # find the y_intercept for a line with a slope of negative point five fit to data set
            # find data point to plot line
            for i in range(len(x)):
                if x[i] > (5.) and (max(x) > 5.):
                    break
                elif x[i] > 1. - TOL:
                    break
                elif math.isnan(y[i+1]):
                    break
            x_log = lin2log(x[i], 10.)
            y_log = lin2log(y[i], 10.)
            A = math.pow(10., y_log - (-0.5 * x_log))
            # assign the y-intercept for the line only if it is the greatest so far
            if A_max_neg is None:
                A_max_neg = A
            elif A_max_neg < A:
                A_max_neg = A
            # find the y_intercept for a line with a slope of one fit to data set
            # find data point to plot line
            for i in range(len(x)):
                if x[i] < (0.5) and (min(x) < 0.5):
                    break
                elif (x[i] <= 1. + TOL) and (not min(x) < 0.5):
                    break
                elif math.isnan(y[i+1]):
                    break
            print(x[i])
            x_log = lin2log(x[i], 10.)
            y_log = lin2log(y[i], 10.)
            A = math.pow(10., y_log - (1. * x_log))
            # assign y-intercept for the line only if it is the greatest so far
            if A_max_pos is None:
                A_max_pos = A
            elif A_max_pos < A:
                A_max_pos = A


    ## plot for hysterseis
    if norm and hys and not added_hys_neg and (A_max_neg is not None):
        A_max_neg = A_max_neg * 1.5
        x_hys_neg = np.linspace(1., x_max, 20)
        y_hys_neg = [power_fit(i, A_max_neg, -0.54) for i in x_hys_neg]
        plt.plot(x_hys_neg, y_hys_neg, ls = '--', color = 'k', label = "$A \\propto \\omega^{{-0.54}}$")
        added_hys_neg = True

    if norm and hys and not added_hys_pos and (A_max_pos is not None):
        A_max_pos = A_max_pos * 1.5
        x_hys_pos = np.linspace(x_min, 1., 20)
        y_hys_pos = [power_fit(i, A_max_pos, 1.) for i in x_hys_pos]
        plt.plot(x_hys_pos, y_hys_pos, ls = ':', color = 'k', label = "$A \\propto \\omega^{{1.}}$")

    ## plot for force extension
    # if the pincus regime has been add, include a label
    if added_pinus:
        plt.plot(x_pincus, y_pincus, ls = '--', color = 'k', label = "Pincus Regime ($R \\propto f^{{2/3}}$)")
    # add hookean regmie to graph
    if added_hookean:
        plt.plot(x_hookean, y_hookean, ls = ':', color = 'k', label = "Hookean Regime ($R \\propto f$)")

    if vertbar is True:
        if flip_axis:
            plt.plot([0., 0.], [x_max, x_min], color = "k", ls = "dashed")
        else:
            plt.plot([0., 0.], [y_max, y_min], color = "k", ls = "dashed")
    if horizbar is True:
        if flip_axis:
            plt.plot([y_min, y_max], [0., 0.], color = "k", ls = "dashed")
        else:
            plt.plot([x_min, x_max], [0., 0.], color = "k", ls = "dashed")

    # finish formatting the graph, once all of the plots have been added
    if flip_axis:
        plt.xlim(y_min, y_max)
        plt.ylim(x_min, x_max)
    else:
        plt.xlim(x_min, x_max)
        # determine min and max ylimits based on check values
        if y_max is None:
            # if a maximum value was not assigned, use the check value
            y_max = y_max_check * 5
        elif (y_max_check > y_max):
            # if the check value is greater than the assigned value, assign the check value
            y_max = y_max_check * 5
        if y_min is None:
            # if a minimum value was not assigned, assign the check value
            y_min = y_min_check * 0.5
        elif (y_min_check < y_min):
            # if the check value is less than the assigned value, assign the check value
            y_min = y_min_check * 0.5
        # set ylimits
        plt.ylim(y_min, y_max)

    if logscale_x:
        plt.xscale('log')
    if logscale_y:
        plt.yscale('log')
    plt.legend(loc = legend_loc)
    if Title is not None:
        plt.title(Title)
    if Y_label is not None:
        if plot_slope:
            Y_label = "Change in " + Y_label
        if flip_axis:
            plt.xlabel(Y_label)
        else:
            plt.ylabel(Y_label)
    if X_label is not None:
        if flip_axis:
            plt.ylabel(X_label)
        else:
            plt.xlabel(X_label)
    if saveas is not None:
        plt.savefig(saveas + ".png", dpi = dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()

    if diff and hys:
        # plot the diffusion correlation against the bending strength
        plt.scatter(x_diff, y_diff)
        plt.xlabel("Bending Strength ($\\kappa_{{\\theta}}$)")
        plt.ylabel("Diffusion ($\\tau_{{R}} / R_{{ROG}}$)")
        if saveas is not None:
            plt.savefig(saveas + "_diff.png", dpi = dpi, bbox_inches='tight')
        if show:
            plt.show()
        plt.close()

# method for plotting the bond bond correlation for bending potential simulation data against s, the bond distance
def plot_bendingpot_bbc (df = None, k = None, max_s = 10, plot_fit = False, saveas = None, logscale_y = True, Title = None, Y_label = None, X_label = None, dpi = None):

    if k is None:
        exit()
    if df is None:
        exit()
    if X_label is None:
        X_label = "Bond-Bond Distance ($s$)"
    if Y_label is None:
        Y_label = "Bond-Bond Correlation ($\\langle cos \\theta \\rangle$)"
    # if Title is None:
    #     Title = "Bond-Bond Correlation\nas a function Bond-Bond Seperation Distance"
    if dpi is None:
        dpi = 200


    fig, ax = plt.subplots()
    # for each k, average bond correlation associated with the bond distance
    for i in k:
        # get the df row associated with k
        kdf = df.loc[df['k'] == i]
        # get the data associated with k
        s = []
        bbc = []
        expect = []
        x = []
        y = []
        for j in range(1, max_s + 1):
            val = kdf.iloc[0]['l' + str(j) + '_M1']
            if val > 0.:
                s.append(j)
                # bbc.append(math.log(val))
                bbc.append(val)
                expect.append(math.exp(-j/i))
                if j < 5:
                    x.append(j)
                    y.append(val)
        ax.scatter(s, bbc, label = '$k_{\\theta}$ = ' + str(i), marker = next(marks), s = 20)
        if plot_fit:
            # fit first two data points to expotential decau curve, determine parameters
            popt, pcov = curve_fit(exp_decay_fit, x, y)
            x_fit = [ ((max(s) - min(s)) / (100 - 1)) * i + min(s) for i in range(100)]
            y_fit = [ exp_decay_fit(i, *popt) for i in x_fit]
            ax.plot(x_fit, y_fit, 'k--')
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width, box.height * 0.9])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), fancybox=True, shadow=True, ncol=4)
    ax.set_xlim(0, max_s)
    ax.set_ylim(.01, 1.)
    if logscale_y:
        ax.set_yscale('log', base = math.e)

    # adjust y-axis labels
    def ticks(y, pos):
        return r'$e^{:.0f}$'.format(np.log(y))
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))

    if Title is not None:
        ax.set_title(Title)
    if Y_label is not None:
        ax.set_ylabel(Y_label)
    if X_label is not None:
        ax.set_xlabel(X_label)
    if saveas is not None:
        fig.savefig(saveas, dpi = dpi, bbox_inches='tight')
    else:
        fig.show()
    plt.close()

# method for plotting persitance length from bending potential simulations
def plot_bending_lp (df = None, plot_expectation_CSA = False, plot_expectation_CA = False, saveas = None, Title = None, Y_label = None, X_label = None, dpi = None, logscale = False, fit_bbc = True, BVDMSE = False, show = False, xmax  = 100, ymax = 100):

    if df is None:
        exit()
    if X_label is None:
        X_label = "Bending Potential Strength ($k$)"
    if Y_label is None:
        Y_label = "Persistence Length ($l_{p} / \\langle l \\rangle$)"
    if dpi is None:
        dpi = 200

    k = df['k'].tolist()
    l1 = df['l1_M1'].tolist()
    y = []
    CS_expect = []
    CSA_expect = []
    fig, ax1 = plt.subplots()
    for i in range(len(k)):
        y.append(-1 / math.log(l1[i]))
    if plot_expectation_CA:
        n = [i for i in range(1,100)]
        n_CA = [i for i in n]
        ax1.plot(n, n_CA, 'r--', label = "$k$")
    if plot_expectation_CSA:
        n = [i for i in range(1, 100)]
        n_CSA = [ math.pow(math.pi * i, 0.5) for i in n]
        ax1.plot(n, n_CSA, 'r--', label = "$(\\pi k)^{{1 / 2}}$")
    ax1.scatter(k, y, s=20, label = "$\\ell_{p, \\theta}$")
    if fit_bbc:
        l_acc = []
        for i in k:
            kdf = df.loc[df['k'] == i]
            x = []
            y = []
            for j in range(1, 6):
                val = kdf.iloc[0]['l' + str(j) + '_M1']
                x.append(j)
                y.append(val)
            popt, pcov = curve_fit(exp_decay_fit, x, y)
            l_acc.append(popt[1])
        ax1.scatter(k, l_acc, s=20, label = "$\\ell_{p, ac}$")
    if BVDMSE:
        bdvmse = df['BVDMSE_M2'].tolist()
        # plot slope already calculated
        ax2 = ax1.twinx()
        ax2.plot(k, bdvmse, 'x', label = "BVD MSE", color = "tab:red")
        ax1.plot([], [], 'x', label = "BVD MSE", color = "tab:red")
        ax2.set_ylabel("Bond Vector Distribution Mean Square Error")
        ax2.set_ylim(0, 0.5)
    ax1.set_xlim(1, xmax)
    ax1.set_ylim(1, ymax)
    if logscale:
        ax1.set_xscale("log")
        ax1.set_yscale("log")
    ax1.legend(loc = "upper left")
    if X_label is not None:
        ax1.set_xlabel(X_label)
    if Y_label is not None:
        ax1.set_ylabel(Y_label)
    if Title is not None:
        ax1.set_title(Title)
    if saveas is not None:
        fig.savefig(saveas, dpi = dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()

# method for plotting simulated parameters R against various calulated persistence lengths
def plot_property_against_lp (df = None, R_col = None, Title = None, Y_label = None, X_label = None, dpi = None, logscale = False, lp_cos = False, lp_exp = False, xmin = None, ymin = None, xmax = None, ymax = None, saveas = None, show = False):

    # check that relevant parameters were passed to method
    if df is None:
        print("ERROR :: plot_property_against_lp :: must pass dataframe to method.")
        exit()
    if (not lp_cos) and (not lp_exp):
        print("ERROR :: plot_property_against_lp :: must specify at least one method for calculating persistence length.")
        exit()
    if R_col is None:
        print("ERROR :: plot_property_against_lp :: must specify property to calculate (R_col).")
        exit()

    # assign defaults
    if dpi is None:
        dpi = default_dpi
    if X_label is None:
        X_label = "Persistence Length ($\\ell_{p}$)"
    if Y_label is None:
        Y_label = "Property ($R$)"

    # parse data from data frame
    r = df[R_col].tolist()
    k = df['k'].tolist()
    if lp_cos:
        l1 = df['l1_M1'].tolist()
        # calculate the persistence length per
    if lp_exp:
        l_acc = []
        for i in k:
            kdf = df.loc[df['k'] == i]
            x = []
            y = []
            for j in range(1, 6):
                val = kdf.iloc[0]['l' + str(j) + '_M1']
                x.append(j)
                y.append(val)
            popt, pcov = curve_fit(exp_decay_fit, x, y)
            l_acc.append(popt[1])
        plt.plot(l_acc, r, 'o', fillstyle = 'full')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if logscale:
        plt.xscale('log')
        plt.yscale('log')
    if X_label is not None:
        plt.xlabel(X_label)
    if Y_label is not None:
        plt.ylabel(Y_label)
    if Title is not None:
        plt.title(Title)
    if saveas is not None:
        plt.savefig(saveas, dpi = dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()

# method for comparing the persistence length of two different bending potentials
def compare_bending_lp_pot (df = None, saveas = None, logscale = False, show = False, xmin = None, xmax = None,  ymin = None, ymax = None, lp_BBC = False, lp_theta = False, top = None, X_label = None, Y_label = None, Title = None, dpi = None):

    # check that the proper parameters have been passed to the method
    if df is None:
        exit()
    if top is None:
        exit()
    elif top == "CHAIN":
        top = False
    elif top == "RING":
        top = True
    else:
        exit()

    # default options
    if X_label is None:
        X_label = "Bending Constant ($k_{{\\theta}}$)"
    if Y_label is None:
        Y_label = "Persistence Length (${{\\ell_{p}}}$)"
    if dpi is None:
        dpi = default_dpi
    if xmin == None:
        xmin = 1.
    if ymin == None:
        ymin = 1.

    # create the data frame
    plotdf = df[df["R"] == top ]
    for p in df["pot"].unique():
        p_label = ""
        if p == "CSA":
            p_label = "Cosine Square Angle Potential"
        elif p == "CA":
            p_label = "Cosine Angle Potential"
        else:
            continue
        # generate the data fram corresponding to the potential and the topology
        potdf = plotdf[plotdf["pot"] == p]
        # calculate the persistence length and add to the plot
        bbc_lp = []
        bendcon = potdf["k"].tolist()
        if lp_BBC:
            for k in bendcon:
                kdf = potdf.loc[potdf['k'] == k]
                x = []
                y = []
                for j in range(1, 6):
                    val = kdf.iloc[0]['l' + str(j) + '_M1']
                    x.append(j)
                    y.append(val)
                popt, pcov = curve_fit(exp_decay_fit, x, y)
                bbc_lp.append(popt[1])
        if lp_theta:
            pass
        plt.scatter(bendcon, bbc_lp, label = p_label, marker = next(marks), s = 40)
    plt.legend(loc = "best")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if logscale:
        plt.xscale("log")
        plt.yscale("log")
    if X_label is not None:
        plt.xlabel(X_label)
    if Y_label is not None:
        plt.ylabel(Y_label)
    if Title is not None:
        plt.title(Title)
    if saveas is not None:
        plt.savefig(saveas, dpi = dpi, bbox_inches='tight')
    if show is not None:
        plt.show()
    plt.close()

def compare_propery_v_lp(df = None, saveas = None, R_col = None, logscale_x = False, logscale_y = False, show = False, xmax = None, ymax = None, ymin = None, xmin = None, lp_BBC = False, lp_theta = False, top = None, X_label = None, Y_label = None, Title = None, dpi = None, fit = False, WLC = False, Lo = 1., m = None, b = None):

    # check that the proper parameters have been passed to the method
    if df is None:
        print("ERROR :: compare_property_v_lp :: Must specify 'df', dataframe which contains simulations results.")
        exit()
    if top is None:
        print("ERROR :: compare_property_v_lp :: Must specify 'top' as either 'CHAIN' or 'RING'.")
        exit()
    if R_col is None:
        print("ERROR :: compare_property_v_lp :: Must specify 'R_col', coloumn which contains simulation property to plot against persistence length.")
        exit()
    elif top == "CHAIN":
        top = False
    elif top == "RING":
        top = True
    else:
        print("ERROR :: compare_property_v_lp :: Unknown topology '" + top + "'.")
        exit()

    # default options
    if X_label is None:
        X_label = "Persistence Length (${{\\ell_{p}}}$)"
    if Y_label is None:
        Y_label = "Polymer Property ($R$)"
    if dpi is None:
        dpi = default_dpi


    # add Worm-like Chain model persistence length scaling, if called
    if WLC:
        x_log = [ ((10. - (-1.)) / (100 - 1)) * i + (-1.) for i in range(100) ] # calculate samples along log scale
        x = [ math.pow(10., i) for i in x_log ] # convert to linear scale
        y = [ WLC_lp2R(i, Lo = Lo) for i in x ]
        plt.plot(x, y, label = "WLC Model", color = 'r')

    # create the data frame
    plotdf = df[df["R"] == top ]
    # plot for each bending potential
    for p in df["pot"].unique():
        p_label = ""
        if p == "CSA":
            p_label = "Cosine Square Angle Potential"
        elif p == "CA":
            p_label = "Cosine Angle Potential"
        else:
            continue
        # generate the data fram corresponding to the potential and the topology
        potdf = plotdf[plotdf["pot"] == p]
        # calculate the persistence length and add to the plot
        bbc_lp = potdf["lp_d_M1"].tolist()
        plt.scatter(bbc_lp, potdf[R_col].tolist(), label = p_label, marker = next(marks), s = 40)
        if fit and p == "CSA":
            # fit data to power law equation, determine parameters
            l = potdf[R_col].tolist()
            x_slope = []
            y_slope = []
            for i in range(len(bbc_lp)):
                if bbc_lp[i] > 1. and bbc_lp[i] < 3.:
                    x_slope.append(bbc_lp[i])
                    y_slope.append(l[i])
            popt, pcov = curve_fit(power_fit, x_slope, y_slope)
            # x_fit = []
            # for i in range(len(bbc_lp)):
            #     if bbc_lp[i] > 1. and bbc_lp[i] < 3.:
            #         x_fit.append(bbc_lp[i])
            x_fit = [ ((max(bbc_lp) - min(bbc_lp)) / (100 - 1)) * i + min(bbc_lp) for i in range(100)]
            y_fit = [ power_fit(i, *popt) for i in x_fit ]
            plt.plot(x_fit, y_fit, '--', label = "Power Law Fit", color = 'k')
            plt.plot([], [], ' ', label="Slope = {:.3f}".format(popt[1]))

    if m is not None and b is not None:
        # plot a line corresponding to y = mx + b
        # min x point
        if xmin is not None:
            x_start = xmin
        else:
            x_start = min(bbc_lp)
        # max x point
        if xmax is not None:
            x_end = xmax
        else:
            x_end = max(bbc_lp)
        y_start = x_start * m + b
        y_end = x_end * m + b
        plt.plot([x_start, x_end], [y_start, y_end], '--', color = 'k')

    # finish plot, show
    plt.legend(loc = "best")
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if logscale_x:
        plt.xscale("log")
    if logscale_y:
        plt.yscale("log")
    if X_label is not None:
        plt.xlabel(X_label)
    if Y_label is not None:
        plt.ylabel(Y_label)
    if Title is not None:
        plt.title(Title)
    if saveas is not None:
        plt.savefig(saveas, dpi = dpi, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()


#############
## CLASSES ##
#############

# none


###################
## CML ARGUMENTS ##
###################

# performing scaling analysis
scaling = ("scaling" in sys.argv)
# performe force extension analysis
forceExtension = ("forceExtension" in sys.argv)
# perform bending parameter analysis
bendingPARM = ("bendingPARM" in sys.argv)
# perform hysteresis analysis
hysteresis = ("hysteresis" in sys.argv)
# update the simulation results, even if the simulation results file already exists
update = ("update" in sys.argv)
# create all figures
figure = ("allfigs" in sys.argv)
# generate / re-generate figure 1
fig1 = ("fig1" in sys.argv)


############
## SCRIPT ##
############

# use pre-analyzed results to generate publication quality figures
if figure or fig1: ## generate figure one
    ## grab data from folder
    df = pd.read_csv("./03_figures/fig1/bendingPARM.csv")

    ## first plot compares bending parameter to persistence length for either beinding potential
    # load data and establish figure parameters
    fig1a = Figure()
    fig1a.load_data(d = df, xcol = 'k', ycol = 'lp_d_M1', ccol = 'k', icol = 'pot')
    fig1a.set_xaxis_label("Bending Constant ($\kappa_{{\\theta}}$)")
    fig1a.set_yaxis_label("Persistence Length ($\ell_{{p}}$)")
    fig1a.set_logscale() # assign all axis to be logscale
    fig1a.set_yaxis_ticks(minval = 1., maxval = 100., nmajorticks = 3, nminorticks = 2)
    fig1a.set_label(ival = "CSA", label = "Cosine Square Angle (CSA) Potential")
    fig1a.set_label(ival = "CA", label = "Cosine Angle (CA) Potential")
    fig1a.reset_markers(['D', 'o'])
    fig1a.set_saveas(savedir = './03_figures/fig1/', filename = 'fig1a')
    # TODO add publication and poster settings
    # TODO add ticks, tick size

    # regular plot
    gen_scatter(fig = fig1a, legendloc = 'upper left', markersize = 60, show = False, save = True)
    # publication plot
    fig1a.set_publication()
    gen_scatter(fig = fig1a, legendloc = 'upper left', markersize = 60, show = False, save = True)

    ## second plot compares the measured persistence length to the cross-over regime
    # gen_scatter(df)
    exit()

# compile and analyze results, as instructed
if scaling:
    # get simulation results and parameters
    scaling_parms = pd.read_csv(scaling_parmcsv)
    # parse results, generate results
    scaling_parms = parse_results(parms = scaling_parms, dir = '01_raw_data/scaling/', simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
    scaling_parms = parse_results(parms = scaling_parms, dir = '01_raw_data/scaling/', simfile = 'RE2E.dat', col = 1, title = 'E2Ex', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
    scaling_parms = parse_results(parms = scaling_parms, dir = '01_raw_data/scaling/', simfile = 'ROG.dat', col = 4, title = 'ROG', M1 = True, var = True, tabsep = True)
    # save the results
    if not os.path.exists("02_processed_data/scaling/"):
        os.mkdir("02_processed_data/scaling/")
    scaling_parms.to_csv("02_processed_data/scaling/scaling_results.csv")
    # perform analysis for each unqiue set of 'modules'
    for mod in scaling_parms['mod'].unique():
        # isolate the parameters corresponding to the module
        mod_parm = scaling_parms[scaling_parms['mod'] == mod]
        # establish filenames, etc.
        save_name = '02_processed_data/scaling/scaling_' + mod
        if mod == "linearChainReal":
            mod = "Real, Linear Polymer Chains"
        elif mod == "linearChainIdeal":
            mod = "Ideal, Linear Polymer Chains"
        # plot end-to-end distance
        plot_scaling (mod_parm, N_col = 'N', R_col = ['E2Etot', 'ROG'], logscale_x = True, logscale_y = True, x_min = 10, x_max = 1000, y_min = 1, y_max = 1500, X_label = "Number of Monomers ($N$)", Y_label = 'Equilibrium Property Value', data_label = ["End-to-End Distance", "Radius of Gyration"], Title = "Scaling for " + mod , saveas = save_name + "_scale.png", fit = True, error = True, color = ['tab:orange', 'tab:blue'])

if forceExtension:
    # lin = ['log', 'lin']
    lin = ['log']
    for l in lin:
        # force extension for linear data set
        job = 'forceExtension_' + l
        data_dir = "01_raw_data/" + job + "/"
        anal_dir = "02_processed_data/" + job + "/"
        parm_file = job + ".csv"
        if update or (not os.path.exists(anal_dir + parm_file)):
            # get simulation parameters
            FE_parms = pd.read_csv(data_dir + parm_file)
            # loop through parameters, relabel ring format
            top = []
            for i, r in FE_parms.iterrows():
                if r['R'] == 0:
                    top.append("CHAIN")
                elif r['R'] == 1:
                    top.append("RING")
                elif r['R'] == 2:
                    top.append("RINGx2")
                else:
                    print("forceExtension :: ERROR :: Unknown top code " + r['R'] + ".")
                    exit()
            FE_parms['top'] = top
            # average simulation properties
            FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
            # save results
            if not os.path.exists(anal_dir):
                os.mkdir(anal_dir)
            FE_parms.to_csv(anal_dir + parm_file)
        else:
            FE_parms = pd.read_csv(anal_dir + parm_file)

        ## MAKE GRAPHS
        # perform analysis for each unique set of bending parameters
        for k in FE_parms['K'].unique():
            for pot in FE_parms['pot'].unique():
                for N in FE_parms['N'].unique():
                    # get the unique df
                    k_df = FE_parms[(FE_parms['K'] == k) & (FE_parms['pot'] == pot) & (FE_parms['N'] == N)]
                    save_name = anal_dir + f"FE_k{k:0.2f}_{pot}_N{N}"
                    # master plot
                    plot_force_extension (k_df, Y_col = 'E2Etot', X_col = 'F', iso_col = 'top', isolabel = '{:s}', X_label = "Normalized External Force ($X \\cdot f$)", Y_label = "Normalized Chain Extension ($X^{{-1}} \\cdot R_{{E2E}}$)", saveas = save_name + '_norm_data', plot_data = True, y_max = 10., y_min = 0.01, x_min = 0.01, x_max = 100., show = True, logscale_x = True, logscale_y = True, plot_slope = False, norm = True, pincus = True, hookean = True) # saveas = save_name + '_norm_data.png',

        # perform analysis for each unqiue set of 'topologies'
        # loop through each topology, bending potential, bending potential strength, and chain length
        for top in FE_parms['top'].unique():
            for pot in FE_parms['pot'].unique():
                for N in FE_parms['N'].unique():

                    # establish data and name
                    top_df = FE_parms[(FE_parms['top'] == top) & (FE_parms['pot'] == pot) & (FE_parms['N'] == N)]
                    if top_df.empty:
                        # if the data frame is empty, skip
                        continue
                    save_name = anal_dir + f"FE_{top}_{pot}_N{N}"

                    # plot force extension data
                    if l == 'lin':
                        plot_force_extension (top_df, Y_col = 'E2Etot', X_col = 'F', iso_col = 'K', isolabel = '$k_{{\\theta}}$ = {:.02f}', X_label = "External Force", Y_label = "Chain Extension", saveas = save_name + '_data', plot_data = True, flip_axis = True, horizbar = True, vertbar = True, y_max = 300, y_min = -300, x_min = -2., x_max = 2., show = True)
                    elif l == 'log':
                        # regular plot
                        plot_force_extension (top_df, Y_col = 'E2Etot', X_col = 'F', iso_col = 'K', isolabel = '$k_{{\\theta}}$ = {:.02f}', X_label = "External Force ($f$)", Y_label = "Chain Extension ($R_{{E2E}}$)", saveas = save_name + '_data', plot_data = True, y_max = 500., y_min = 0.1, x_min = 0.0001, x_max = 5., show = True, logscale_x = True, logscale_y = True, plot_slope = False)

                        # master plot
                        plot_force_extension (top_df, Y_col = 'E2Etot', X_col = 'F', iso_col = 'K', isolabel = '$k_{{\\theta}}$ = {:.02f}', X_label = "Normalized External Force ($X \\cdot f$)", Y_label = "Normalized Chain Extension ($X^{{-1}} \\cdot R_{{E2E}}$)", saveas = save_name + '_norm_data', plot_data = True, y_max = 10., y_min = 0.01, x_min = 0.01, x_max = 100., show = True, logscale_x = True, logscale_y = True, plot_slope = False, norm = True, pincus = True, hookean = True)

if bendingPARM:
    tag = ['Ideal', 'Real']
    for t in tag:
        # force extension for linear data set
        job = 'bendingPARM_' + t
        data_dir = "01_raw_data/" + job + "/"
        anal_dir = "02_processed_data/" + job + "/"
        parm_file = job + ".csv"

        if update or (not os.path.exists(anal_dir + parm_file)):
            # load parameters, add property calculations from simulations
            bendingparms = pd.read_csv(data_dir + parm_file)

            ## PROCESS DATA
            check_skq(parms = bendingparms, dir = data_dir, show = False, plot = True, porod = True, real = (tag == 'REAL'), monotonic = True)
            # calculate the bvd error, plot difference if requested
            check_bvd(parms = bendingparms, dir = data_dir, plot = False)
            # parse persistence length from BBC, plot if requested
            check_bbc(parms = bendingparms, dir = data_dir, show = False, plot = True, lp_decay = True, lp_theta = True)
            # parse scattering factor slope from SKQ, plot if requested
            # check_skq(parms = bendingparms, dir = data_dir, show = False, plot = True, porod = True, real = (tag == 'REAL'), monotonic = False)

            ## PARSE DATA, ADD TO RESULTS DATAFRAME
            # get the end-to-end distance vector data
            bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
            # get radius of gyration data
            bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'ROG.dat', col = 4, title = 'ROG', M1 = True, bootstrapping = False, tabsep = True) # don't bootsrap ROG data, does not fit normal distribution
            # get the BVD error
            bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'BVD.csv', col = 3, title = 'BVDMSE', M1 = True, M2 = True)
            # get persistence length
            bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'BBC.csv', col = 2, title = 'lp_d', M1 = True)
            # get the scattering factor slope
            bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'SKQ.csv', col = 2, title = 'm_sq', M1 = True)
            # get scattering factor cross-over regime
            bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'SKQ.csv', col = 3, title = 'lp_xc', M1 = True)
            # collect the bond-bond correlation at different bond distances
            for i in range(30):
                bendingparms = parse_results(parms = bendingparms, dir = data_dir, simfile = 'BBC.dat', row = i, title = "l"+str(i), M1 = True)

            # save results
            if not os.path.exists(anal_dir):
                os.mkdir(anal_dir)
            bendingparms.to_csv(anal_dir + parm_file)
        else:
            bendingparms = pd.read_csv(anal_dir + parm_file)

        ## PLOT RESULTS
        if t == "Ideal":
            # plot the end to end distance against the measured persistence length
            compare_propery_v_lp(df = bendingparms,  R_col = 'E2Etot_M1', logscale_x = True, logscale_y = True, show = False, xmax = 1000, ymax = 1000, ymin = 30, xmin = .1, lp_BBC = True, top = "CHAIN", saveas = anal_dir + "CHAIN_RE2E_lp_BBC_WLC.png", Y_label = "End-to-End Distance ($R_{{E2E}}$)", fit = True, Lo = 600., WLC = True)
            compare_propery_v_lp(df = bendingparms, R_col = 'm_sq_M1', show = False, top = 'CHAIN', logscale_x = True, Y_label = "Scattering Factor", saveas = anal_dir + "CHAIN_SKQ_lp.png", xmin = .1, xmax = 200.)
            compare_propery_v_lp(df = bendingparms, R_col = 'lp_xc_M1', show = False, top = 'CHAIN', logscale_x = True, logscale_y = True, Y_label = "Estimated Cross Over Regime", saveas = anal_dir + "CHAIN_XC_lp.png", xmin = .1, xmax = 200., ymin = 0.1, m = 1., b = 0.)
        elif t == "Real":
            pass
        else:
            print("ERROR :: bendingPARM :: Unknown tag " + t + ".")

    exit()
    ## PLOT RESULTS
    # plot the rod-like transition for chains with either CA or CSA potentials
    compare_propery_v_lp(df = bendingparms, R_col = 'm_sq_M1', show = False, top = 'CHAIN', logscale_x = True, Y_label = "Scattering Factor Slope", saveas = "02_processed_data/bendingPARM/" + "CHAIN_SKQ_lp.png", xmin = .1, xmax = 80.)
    # compare the end-to-end distance scaling against the persistence length
    compare_propery_v_lp(df = bendingparms,  R_col = 'E2Etot_M1', logscale_x = True, logscale_y = True, show = False, xmax = 100, ymax = 300, ymin = 10, xmin = 1, lp_BBC = True, top = "CHAIN", saveas = "02_processed_data/bendingPARM/" + "CHAIN_RE2E_lp_BBC.png", Y_label = "End-to-End Distance ($R_{{E2E}}$)")
    compare_propery_v_lp(df = bendingparms,  R_col = 'E2Etot_M1', logscale_x = True, logscale_y = True, show = False, xmax = 100, ymax = 300, ymin = 10, xmin = 1, lp_BBC = True, top = "CHAIN", saveas = "02_processed_data/bendingPARM/" + "CHAIN_RE2E_lp_BBC_fit.png", Y_label = "End-to-End Distance ($R_{{E2E}}$)", fit = True)
    compare_propery_v_lp(df = bendingparms,  R_col = 'E2Etot_M1', logscale_x = True, logscale_y = True, show = False, xmax = 100, ymax = 300, ymin = 10, xmin = 1, lp_BBC = True, top = "CHAIN", saveas = "02_processed_data/bendingPARM/" + "CHAIN_RE2E_lp_BBC_WLC.png", Y_label = "End-to-End Distance ($R_{{E2E}}$)", fit = True, Lo = 300., WLC = True)

    # compare the persistence length for CSA and CA potentials
    compare_bending_lp_pot(df = bendingparms,logscale = True, show = True, xmax = 100, ymax = 200, lp_BBC = True, top = "CHAIN", saveas = "02_processed_data/bendingPARM/" + "CHAIN_lp_BBC.png", dpi = 300, ymin = 0.1, xmin = 0.1)
    compare_propery_v_lp(df = bendingparms,  R_col = 'BVDMSE_M2', logscale_x = True, logscale_y = True, show = True, xmax = 100, ymax = 1., ymin = .001, xmin = 1, lp_BBC = True, top = "CHAIN", saveas = "02_processed_data/bendingPARM/" + "CHAIN_BVD_lp_BBC.png", Y_label = "Bond Vector Distribution Mean Square Difference")

    # loop through each potential, type of structure (ring or chain)
    for pot in bendingparms['pot'].unique():
        for ring in bendingparms['R'].unique():
            # estblish title and save name depending on conditions
            if ring == True:
                chain = "Ring"
                N_string = "200"
            else:
                chain = "Chain"
                N_string = "100"
            plotdf = bendingparms.loc[(bendingparms['pot'] == pot) & (bendingparms['R'] == ring)]
            # plot the bending parameter scaling against s
            title_bbc = f"Bond-Bond Correlation for {chain}s with {pot} Potential (N = {N_string})"
            saveas_bbc = "BBC_" + pot + "_" + chain + ".png"
            # plot_bendingpot_bbc(df = plotdf, k = [1, 3, 5, 7, 10, 13, 20], max_s = 29, plot_fit = True, saveas = "02_processed_data/bendingPARM/" + saveas_bbc)
            # plot_bendingpot_bbc(df = plotdf, k = [1, 3, 5, 7, 10, 13, 20], max_s = 29, plot_fit = False, saveas = "02_processed_data/bendingPARM/" + "BBCnf_" + pot + "_" + chain + ".png")
            # plot the bending parameter scaling against k
            title_lp = f"Persistence Length for {chain}s with {pot} Potential (N = {N_string})"
            saveas_lp = "lp_" + pot + "_" + chain + ".png"
            plot_bending_lp(df = plotdf, plot_expectation_CSA = (pot == "CSA"), plot_expectation_CA = (pot == "CA"), Title = title_lp, saveas = "02_processed_data/bendingPARM/" + saveas_lp, logscale = True, BVDMSE = False, show = True, ymax = 200)
            # TODO plot the chain length against the persistence length
            plot_property_against_lp(df = plotdf, R_col = 'E2Etot_M1', lp_exp = True, logscale = True, Y_label = "End-to-End Distance", xmin = 1., xmax = 100., ymin = 10., ymax = 300., Title = f"End-to-End Scaling for {chain}s with {pot} Potential (N = {N_string})", saveas = "02_processed_data/bendingPARM/" + "Rvlp_" + pot + "_" + chain + ".png", show = True)
            # TODO plot the BVD MSE against the persistence length
            plot_property_against_lp(df = plotdf, R_col = 'BVDMSE_M2', lp_exp = True, logscale = True, Y_label = "Bond Vector Distribution Mean Square Difference", xmin = 1., xmax = 100., ymin = 0.001, ymax = .5, Title = f"Bond Vector Difference Scaling for {chain}s with {pot} Potential (N = {N_string})", show = True, saveas = "02_processed_data/bendingPARM/" + "BVDMSEvlp_" + pot + "_" + chain + ".png")

if hysteresis:

    # establish the names of the simulation results files
    job = "hysteresis"
    data_dir = "01_raw_data/" + job + "/"
    anal_dir = "02_processed_data/" + job + "/"
    parm_file = job + ".csv"

    # collect the simulaiton results
    if update or (not os.path.exists(anal_dir + parm_file)):

        # get simulation parameters
        hys_parms = pd.read_csv(data_dir + parm_file)
        hys_parms = period2freq (parms = hys_parms, period_col = 'T', timescale = 1.)
        # loop through parameters, relabel ring format
        top = []
        for i, r in hys_parms.iterrows():
            if r['R'] == 0:
                top.append("CHAIN")
            elif r['R'] == 1:
                top.append("RING")
            elif r['R'] == 2:
                top.append("RINGx2")
            else:
                print("forceExtension :: ERROR :: Unknown top code " + r['R'] + ".")
                exit()
        hys_parms['top'] = top

        # collect hysteresis results
        check_hys(parms = hys_parms, dir = data_dir, show = False, plot = True, check = True, clear = True)

        # add the hysteresis results to the bending parameter data frame
        hys_parms = parse_results(parms = hys_parms, dir = data_dir, simfile = 'HYS.csv', col = 1, title = 'A_avg', M1 = True)
        hys_parms = parse_results(parms = hys_parms, dir = data_dir, simfile = 'HYS.csv', col = 2, title = 'A_std', M1 = True)
        hys_parms = parse_results(parms = hys_parms, dir = data_dir, simfile = 'HYS.csv', col = 3, title = 'A_sl', M1 = True)

        # save results
        if not os.path.exists(anal_dir):
            os.mkdir(anal_dir)
        hys_parms.to_csv(anal_dir + parm_file)
    else:
        hys_parms = pd.read_csv(anal_dir + parm_file)


    # for each topology, plot the hysteresis for each of the bending potentials
    for fa in hys_parms['fA'].unique():
        for ring in hys_parms['R'].unique():
            plotdf = hys_parms.loc[(hys_parms['R'] == ring) & (hys_parms['fA'] == fa)]
            # file for plotting hysteresis against actual frequency
            plot_force_extension(plotdf, hys = True, Y_col = 'A_sl', X_col = 'frequency', iso_col = 'k', isolabel = '$k_{{\\theta}}$ = {:.02f}', logscale_x = True, logscale_y = True, show = True, y_min = 100., y_max = 500., x_min = .00000001, x_max = 0.001, saveas = anal_dir + f"HYS_R{ring}_fa{fa}", X_label = "Frequency ($\\omega$)", Y_label = "Hysteresis ($|A|$)")
            # file for plotting hysteresis against normalized frequency
            plot_force_extension(plotdf, norm = True, hys = True, Y_col = 'A_sl', X_col = 'frequency', iso_col = 'k', isolabel = '$k_{{\\theta}}$ = {:.02f}', logscale_x = True, logscale_y = True, show = True, y_min = 0.5, y_max = 10., x_min = .001, x_max = 1000, saveas = anal_dir + f"HYS_R{ring}_fa{fa}_norm", X_label = "Normalized Frequency ($\\omega \\cdot \\tau_{R}$)", Y_label = "Normalized Hysteresis ($|A| \\cdot f_{{A}}^{{2}} \\cdot R^{{-1}}_{{0}}$)", fa = fa, diff = True) #, diff = (ring == 0)

    # for each bending strength, compare the different topologies
    for fa in hys_parms['fA'].unique():
        for k in hys_parms['k'].unique():
            plotdf = hys_parms.loc[(hys_parms['k'] == k) & (hys_parms['fA'] == fa)]
            save_name = anal_dir + f"HYS_k{k:0.2f}"
            # file for plotting hysteresis against actual frequency
            plot_force_extension(plotdf, hys = True, Y_col = 'A_sl', X_col = 'frequency', iso_col = 'top', isolabel = '{:s}', logscale_x = True, logscale_y = True, show = True, y_min = 1., y_max = 500., x_min = .00000001, x_max = 0.001, saveas = save_name, X_label = "Frequency ($\\omega$)", Y_label = "Hysteresis ($|A|$)")
            # file for plotting hysteresis against normalized frequency
            plot_force_extension(plotdf, norm = True, hys = True, Y_col = 'A_sl', X_col = 'frequency', iso_col = 'top', isolabel = '{:s}', logscale_x = True, logscale_y = True, show = True, y_min = 0.005, y_max = 10., x_min = .001, x_max = 1000, saveas = save_name + "_norm", X_label = "Normalized Frequency ($\\omega \\cdot \\tau_{R}$)", Y_label = "Normalized Hysteresis ($|A| \\cdot f_{{A}}^{{2}} \\cdot R^{{-1}}_{{0}}$)", fa = fa)
