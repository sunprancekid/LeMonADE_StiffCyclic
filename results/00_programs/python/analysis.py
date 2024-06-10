## PACKAGES
import sys, os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# TODO :: add recombinator library for boot strapping

## PARAMETERS
# file that contains parameters for scaling simulations
scaling_parmcsv = '01_raw_data/scaling/scaling.csv'
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
    return df[avgcol].mean()

# method for getting simulation results from files
def parse_scaling_results(parms, simfile, col, title):
    # initialize empty arrays
    vals = []
    # loop through all simulation directories
    for i, r in parms.iterrows():
        clean_file("01_raw_data/scaling/" + r['path'] + simfile)
        vals.append(parse_data("01_raw_data/scaling/" + r['path'] + simfile, avgcol = col))

    # add to parameters data frame and return to user
    parms[title] = vals
    return parms

# method for getting scaling simulation data (N vs. R)
def plot_scaling(parms, N_col = None, R_col = None, fit = False, Title = None, X_label = None, Y_label = None, dpi = None, logscale = False, x_min = None, y_min = None, x_max = None, y_max = None, saveas = None):
    # make sure that the proper parameters were passed to the method
    if N_col is None:
        exit()
    if R_col is None:
        exit()
    # establish image parameters based on parameters passed to method
    if dpi is None:
        dpi = default_dpi
    # get values, average duplicates
    x = []
    y = []
    for u in parms[N_col].unique():
        x.append(u)
        y_mean = parms[parms[N_col] == u]
        y.append(y_mean[R_col].mean())
    # plot values
    plt.plot(x, y, 'o', label = "Simulation Data", fillstyle = 'none')
    # fit data to power low, if specified
    if fit:
        # fit data to power law equation, determine parameters
        popt, pcov = curve_fit(power_fit, x, y)
        x_fit = [ ((max(x) - min(x)) / (100 - 1)) * i + min(x) for i in range(100)]
        y_fit = [ power_fit(i, *popt) for i in x_fit]
        plt.plot(x_fit, y_fit, '--', label = "Power Fit ($y = A \cdot x^{{B}}$)")
        plt.plot([], [], ' ', label="A = {:.3f}".format(popt[0]))
        plt.plot([], [], ' ', label="B = {:.3f}".format(popt[1]))
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xscale('log')
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

## CLASSES
# none


## ARGUMENTS
# performing scaling analysis
scaling = ("scaling" in sys.argv)

## SCRIPT
# load defaults / settings from yaml

# parse arguments from command line

# generate directories, as instructed

# compile and analze results, as instructed
if scaling:
    # get simulation results and parameters
    scaling_parms = pd.read_csv(scaling_parmcsv)
    # TODO :: add boot strapping
    scaling_parms = parse_scaling_results(scaling_parms, 'RE2E.dat', 4, 'E2E')
    scaling_parms = parse_scaling_results(scaling_parms, 'ROG.dat', 4, 'ROG')
    # perform analysis for each unqiue set of 'modules'
    for mod in scaling_parms['mod'].unique():
        # isolate the parameters corresponding to the module
        mod_parm = scaling_parms[scaling_parms['mod'] == mod]
        # establish filenames, etc.
        save_name = '01_raw_data/scaling/scaling_' + mod
        if mod == "linearChainReal":
            mod = "Real, Linear Polymer Chains"
        elif mod == "linearChainIdeal":
            mod = "Ideal, Linear Polymer Chains"
        # plot end-to-end distance
        plot_scaling (mod_parm, N_col = 'N', R_col = 'E2E', logscale = True, x_min = 10, x_max = 1000, y_min = 10, y_max = 100, X_label = "Number of Monomers ($N$)", Y_label = "End-to-End Distance", Title = "End-to-End Distance Scaling for " + mod , saveas = save_name + "_e2e.png", fit = True)
        # plot radius of gyration for real chains
        plot_scaling (mod_parm, N_col = 'N', R_col = 'ROG', logscale = True, x_min = 10, x_max = 1000, y_min = 10, y_max = 1500, X_label = "Number of Monomers ($N$)", Y_label = "Radius of Gyration", Title = "Radius of Gyration Scaling for " + mod , saveas = save_name + "rog.png", fit = True)

# analyze results, as instructed

# generate graphs, as instructed
