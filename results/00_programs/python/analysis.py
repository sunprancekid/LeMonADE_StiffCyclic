## PACKAGES
import sys, os, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm
# TODO :: add recombinator library for boot strapping

## PARAMETERS
# file that contains parameters for scaling simulations
scaling_parmcsv = '01_raw_data/scaling/scaling.csv'
# file that contains parameter for forceExtension simulations
forceExtension_parmcsv = '01_raw_data/forceExtension/forceExtension.csv'
# default image resolution
default_dpi = 200


## METHODS
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
def parse_data(filepath, avgcol = None, header = None):

    # check that the file exists
    file = Path(filepath)
    if not file.exists():
        print(filepath + " does not exist!")
        return 0

    # load the file as a data frame
    df = pd.read_csv(filepath, sep = '\t', header = header)
    return df[avgcol].to_list()

# method for getting simulation results from files
def parse_results(parms = None, dir = None, simfile = None, col = None, title = None, bootstrapping = False, M1 = False, M2 = False, var = False):

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
        if not os.path.exists(dir + r['path'] + simfile):
            # if the path does not exists, add NaN to each array
            M1_arr.append(np.nan)
            M2_arr.append(np.nan)
            var_arr.append(np.nan)
            # skip to next entry
            continue
        clean_file(dir + r['path'] + simfile)
        vals = parse_data(dir + r['path'] + simfile, avgcol = col)
        if title == 'ROG':
            vals = [math.sqrt(i) for i in vals]
        # calculate the first moment if, specified
        avg = 0
        if bootstrapping:
            # fit data to normal distribution, determine average
            n_bins = 50
            mu, std = norm.fit(vals)
            avg = mu
            fig, ax = plt.subplots()
            ax.hist(vals, bins=n_bins, density=True, label = "Simulation Data")
            xmin, xmax = plt.xlim()
            x = np.linspace(xmin, xmax, 100)
            p = norm.pdf(x, mu, std)
            plt.plot(x, p, 'k', linewidth=2, label = "Normal Distribution")
            plt.plot([], [], ' ', label = "$\mu$ = %.3f,  $\sigma$ = %.3f" % (mu, std))
            plt.title("Distribution of " + title + " Values ($n_{{bins}}$ = {:d})".format(n_bins))
            plt.legend()
            plt.xlabel("Property Value")
            plt.ylabel("Distribution")
            plt.savefig(dir + r['path'] + "propdist_" + title + ".png", dpi = 200)
            plt.close()
        else:
            # calculate using the formula
            for i in vals:
                avg += i
            avg = avg / len(vals)
        if M1:
            # if specified, append to array
            M1_arr.append(avg)

        # calculate the second moment, if specified
        if M2:
            avg2 = 0
            for i in vals:
                avg2 += pow(i,2)
            avg2 = avg2 / len(vals)
            M2_arr.append(math.sqrt(avg2))

        # calculate the variance, if specified
        if var:
            if not bootstrapping:
                # if boot strapping was not used
                std = 0
                for i in vals:
                    std += pow(i - avg, 2)
                std = std / (len(vals) - 1)
            var_arr.append(std)

    # add to parameters data frame and return to user
    if M1: parms[title + "_M1"] = M1_arr
    if M2: parms[title + "_M2"] = M2_arr
    if var: parms[title + "_var"] = var_arr
    return parms

# method for getting scaling simulation data (N vs. R)
def plot_scaling(parms, N_col = None, R_col = None, fit = False, Title = None, X_label = None, Y_label = None, data_label = None, dpi = None, logscale_x = False, logscale_y = False, x_min = None, y_min = None, x_max = None, y_max = None, saveas = None, error = False, color = None, fit_hook = False):
    # make sure that the proper parameters were passed to the method
    if N_col is None:
        exit()
    if R_col is None:
        exit()
    # establish image parameters based on parameters passed to method
    if dpi is None:
        dpi = default_dpi
    for c in range(len(R_col)):
        # get values, average duplicates
        x = []
        y = []
        v = []
        for u in parms[N_col].unique():
            x.append(u)
            y_mean = parms[parms[N_col] == u]
            y.append(y_mean[R_col[c] + '_M1'].mean())
            if error:
                v.append(y_mean[R_col[c] + '_var'].mean() / 2)
        # plot values
        if error:
            # plot with error bars
            plt.errorbar(x, y, yerr = v, label = data_label[c], fmt = 'o', color = color[c])
        else:
            plt.plot(x, y, 'o', label = data_label[c], fillstyle = 'none', color = color[c])
        # fit data to power low, if specified
        if fit:
            # fit data to power law equation, determine parameters
            popt, pcov = curve_fit(power_fit, x, y)
            x_fit = [ ((max(x) - min(x)) / (100 - 1)) * i + min(x) for i in range(100)]
            y_fit = [ power_fit(i, *popt) for i in x_fit]
            plt.plot(x_fit, y_fit, '--', label = data_label[c] + " Fit", color = color[c])
            # plt.plot([], [], ' ', label="A = {:.3f}".format(popt[0]))
        if fit_hook:
            # fit data to power law equation, determine parameters
            A = y[15] / x[15]
            x_fit = [ ((max(x) - min(x)) / (100 - 1)) * i + min(x) for i in range(100)]
            y_fit = [ power_fit_hook(i, A) for i in x_fit]
            plt.plot(x_fit, y_fit, '--', label = "Hookean Regime", color = "black", alpha = 0.5)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    if logscale_x:
        plt.xscale('log')
    if logscale_y:
        plt.yscale('log')
    plt.legend(loc = 'upper left')
    if Title is not None:
        plt.title(Title)
    if Y_label is not None:
        plt.ylabel(Y_label)
    if X_label is not None:
        plt.xlabel(X_label)
    if saveas is not None:
        plt.savefig(saveas, dpi = dpi, bbox_inches='tight')
    plt.show()


# model equation for power law fitting
def power_fit(x, A, B):
    return A * (x ** B)

# model equation for power law fitting in the hookean regime
def power_fit_hook(x, A):
    return A * (x ** 1.)

## CLASSES
# none


## ARGUMENTS
# performing scaling analysis
scaling = ("scaling" in sys.argv)
# performe force extension analysis
forceExtension = ("forceExtension" in sys.argv)

## SCRIPT
# load defaults / settings from yaml

# parse arguments from command line

# generate directories, as instructed

# compile and analze results, as instructed
if scaling:
    # get simulation results and parameters
    scaling_parms = pd.read_csv(scaling_parmcsv)
    # TODO :: add boot strapping
    scaling_parms = parse_results(parms = scaling_parms, dir = '01_raw_data/scaling/', simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True)
    scaling_parms = parse_results(parms = scaling_parms, dir = '01_raw_data/scaling/', simfile = 'RE2E.dat', col = 1, title = 'E2Ex', M1 = True, M2 = True, var = True, bootstrapping = True)
    scaling_parms = parse_results(parms = scaling_parms, dir = '01_raw_data/scaling/', simfile = 'ROG.dat', col = 4, title = 'ROG', M1 = True, var = True)
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
        # plot radius of gyration for real chains
        # plot_scaling (mod_parm, N_col = 'N', R_col = 'ROG_M1', logscale = True, x_min = 10, x_max = 1000, y_min = 10, y_max = 1500, X_label = "Number of Monomers ($N$)", Y_label = "Radius of Gyration", Title = "Radius of Gyration Scaling for " + mod , saveas = save_name + "_rog.png", fit = True)

        # combined plots


if forceExtension:
    # get simulation parameters
    FE_parms = pd.read_csv(forceExtension_parmcsv)
    # average simulation properties
    FE_parms = parse_results(parms = FE_parms, dir = '01_raw_data/forceExtension/', simfile = 'RE2E.dat', col = 1, title = 'E2Ex', M1 = True, M2 = True, var = True, bootstrapping = True)
    # save results
    if not os.path.exists("02_processed_data/forceExtension/"):
        os.mkdir("02_processed_data/forceExtension/")
    FE_parms.to_csv("02_processed_data/forceExtension/forceExtension.csv")
    # perform analysis for each unqiue set of 'modules'
    for mod in FE_parms['mod'].unique():
        # isolate the parameters corresponding to the module
        mod_parm = FE_parms[FE_parms['mod'] == mod]
        # establish filenames, etc.
        save_name = '02_processed_data/forceExtension/forceExtension_' + mod
        if mod == "FElinearChainReal":
            mod = "Real, Linear Polymer Chains"
        elif mod == "FElinearChainIdeal":
            mod = "Ideal, Linear Polymer Chains"
        # plot x_distance
        plot_scaling(mod_parm, N_col = 'F', R_col = ['E2Ex'], logscale_x = True, logscale_y = True, x_min = .00008, x_max = 10., y_min = 0.01, y_max = 500, X_label = "External Force ($f_{x}$)", Y_label = "Distance", data_label = ["Chain Extension in X-direction"], Title = "Force Extension for " + mod + " (N = 100)", saveas = save_name + '.png', error = False, color = ['tab:purple'], fit_hook = True)

