## PACKAGES
import sys, os, math
import csv
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit
from scipy.stats import norm

## PARAMETERS
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
# equildict = {'vec' : ["$|2 0 0|$", "$|1 2 0|$", "$|2 1 1|$", "$|3 0 0|$", "$|2 2 1|$", "$|3 1 0|$"], 'height': [0.0632855, 0.206359, 0.215854, 0.0550903, 0.231382, 0.228029]}
equildict = {'vec' : [0, 1, 2, 3, 4, 5], 'height': [0.0632855, 0.206359, 0.215854, 0.0550903, 0.231382, 0.228029]}


## METHODS

## USED FOR EQUATION FITTING
# model equation for power law fitting
def power_fit(x, A, B):
    return A * (x ** B)

# model equation for power law fitting in the hookean regime
def power_fit_hook(x, A):
    return A * (x ** 1.)

# model equation for exponential decay fitting
def exp_decay_fit(x, A, B):
    return A * (math.e ** (- x / B))

# model equation for logistic5 equation
def log5_fit (x, A, B, C, D, E, F):
    return A + ((B - A) / ((1. + (C / (x - F)) ** D) ** E))

## USED FOR PARSING SIMULATION RESULTS FROM FILES
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
def parse_results(parms = None, dir = None, simfile = None, col = None, row = None, title = None, bootstrapping = False, M1 = False, M2 = False, var = False, tabsep = False):

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
            avg = avg / nvals
        if M1:
            # if specified, append to array
            M1_arr.append(avg)

        # calculate the second moment, if specified
        if M2 and (nvals > 2):
            avg2 = 0
            for i in vals:
                avg2 += pow(i,2)
            avg2 = avg2 / len(vals)
            M2_arr.append(math.sqrt(avg2))

        # calculate the variance, if specified
        if var and (nvals > 2):
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

# method that generates bond vector diagrams for simulations which have generated that file
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
            ax.bar(sim_bar_xloc, df['height'], label = "Simulation", color = "tab:orange", edgecolor = "black", width = bar_width)
            ax.bar(equil_bar_xloc, equildf['height'], label = "Real Polymer Chain", color = "c", edgecolor = "black", width = bar_width)
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

## USED FOR PLOTTING SIMULATION RESULTS
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
            y_fit = [ log5_fit(i, *popt) for i in x_fit]
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
def plot_force_extension(parms, X_col = None, Y_col = None, iso_col = None, isolabel = None, logscale_x = False, logscale_y = False, x_min = None, x_max = None, y_min = None, y_max = None, X_label = None, Y_label = None, Title = None, saveas = None, plot_log5 = False, plot_data = False, plot_slope = False, dpi = None, M2 = False, error = False, flip_axis = False):

    if X_col is None or Y_col is None or iso_col is None:
        print("plot_force_extension :: Must specift columns from data frame to plot!")
        exit()

    # establish image parameters based on parameters passed to method
    if dpi is None:
        dpi = default_dpi

    # parse data from data frame, plot
    iso_vals = parms[iso_col].unique()
    for iso in iso_vals:
        # for each unique isolated value
        iso_df = parms[parms[iso_col] == iso]
        x = iso_df[X_col].tolist()
        if M2:
            y = iso_df[Y_col + '_M2'].tolist()
        else:
            y = iso_df[Y_col + '_M1'].tolist()
        if error:
            v = iso_df[Y_col + '_var'].tolist()
        if logscale_y is True:
            y = [abs(i) for i in y]
        # fit the data to a log5 curve
        if plot_log5 or plot_slope:
            # convert data to log scale
            x_log = []
            y_log = []
            v_log = []
            for i in range(len(x)):
                if abs(y[i]) < TOL or math.isnan(y[i]):
                    continue
                x_log.append(math.log(x[i]) / math.log(10.))
                y_log.append(math.log(abs(y[i])) / math.log(10.))
                if error:
                    if v[i] > TOL:
                        v_log.append(v[i])
                    else:
                        v_log.append(TOL)
            # initialize paramaters for curve fitting
            lb = [min(y_log), max(y_log) - (TOL / 100), -100, -100, -40, min(x_log) - 2 * TOL]
            ub = [max(y_log), max(y_log), 200, 100, 40, min(x_log) - TOL]
            parameterBounds = (lb, ub)
            init = [min(y_log), max(y_log), 4., 50., 0.01, min(x_log) - (3. * TOL / 2.)]
            # fit data to log5 equation
            popt, pcov = curve_fit(f = log5_fit, xdata =  x_log, ydata = y_log, p0 = init, bounds = parameterBounds, maxfev = 100000000) # sigma = v_log,
            x_fit = [ ((max(x_log) - min(x_log)) / (100 - 1)) * i + min(x_log) for i in range(100)]
            y_fit = [ log5_fit(i, *popt) for i in x_fit]
            # calculate the slope of the line
            if plot_slope:
                h = (max(x_fit) - min(x_fit)) / len(x_fit) # distance between points
                x_d1 = []
                y_d1 = []
                for i in range(len(y_fit) - 1):
                    x_d1.append(math.pow(10, (x_fit[i+1] + x_fit[i]) / 2.))
                    y_d1.append((y_fit[i+1] - y_fit[i]) / h)
            # convert the log5 curve back to linear scale
            for i in range(len(x_fit)):
                x_fit[i] = math.pow(10, x_fit[i])
                y_fit[i] = math.pow(10, y_fit[i])
        # plot either the slope or the simulation data
        if plot_slope:
            # plot the slope of the line fit to the simulation data
            plt.plot(x_d1, y_d1, '-', label = isolabel.format(iso))
        else:
            # plot the simulation data (and the log5 fit, if called)
            if plot_data and plot_log5:
                line = plt.plot(x_fit, y_fit, '--')
                if flip_axis:
                    plt.plot(y, x, 'o', label = isolabel.format(iso), fillstyle = 'none', color = line[-1].get_color())
                else:
                    plt.plot(x, y, 'o', label = isolabel.format(iso), fillstyle = 'none', color = line[-1].get_color())
            elif plot_log5:
                line = plt.plot(x_fit, y_fit, '--', label = isolabel.format(iso))
            elif plot_data:
                if flip_axis:
                    plt.plot(y, x, 'o', label = isolabel.format(iso), fillstyle = 'full')
                else:
                    plt.plot(x, y, 'o', label = isolabel.format(iso), fillstyle = 'full')
    # finish formatting the graph, once all of the plots have been added
    if flip_axis:
        plt.xlim(y_min, y_max)
        plt.ylim(x_min, x_max)
    else:
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
    if logscale_x:
        plt.xscale('log')
    if logscale_y:
        plt.yscale('log')
    if plot_slope:
        plt.legend(loc = 'upper left')
    else:
        plt.legend(loc = 'lower right')
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
        plt.savefig(saveas, dpi = dpi, bbox_inches='tight')
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
    if Title is None:
        Title = "Bond-Bond Correlation\nas a function Bond-Bond Seperation Distance"
    if dpi is None:
        dpi = 200

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
        plt.scatter(s, bbc, label = '$k_{\\theta}$ = ' + str(i), marker = 'o', s = 20)
        if plot_fit:
            # fit first two data points to expotential decau curve, determine parameters
            popt, pcov = curve_fit(exp_decay_fit, x, y)
            x_fit = [ ((max(s) - min(s)) / (100 - 1)) * i + min(s) for i in range(100)]
            y_fit = [ exp_decay_fit(i, *popt) for i in x_fit]
            plt.plot(x_fit, y_fit, 'k--')
    plt.legend()
    plt.xlim(0, max_s)
    plt.ylim(.01, 1.)
    if logscale_y:
        plt.yscale('log', base = math.e)

    # adjust y-axis labels
    def ticks(y, pos):
        return r'$e^{:.0f}$'.format(np.log(y))
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(ticks))

    if Title is not None:
        plt.title(Title)
    if Y_label is not None:
        plt.ylabel(Y_label)
    if X_label is not None:
        plt.xlabel(X_label)
    if saveas is not None:
        plt.savefig(saveas, dpi = dpi, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

# method for plotting presitance length from bending potential simulations
def plot_bending_lp (df = None, plot_expectation_CSA = False, plot_expectation_CA = False, saveas = None, Title = None, Y_label = None, X_label = None, dpi = None, logscale = False, fit_bbc = True, BVDMSE = False, show = False):

    if df is None:
        exit()
    if X_label is None:
        X_label = "Bending Potential Strength ($k$)"
    if Y_label is None:
        Y_label = "Presistance Length ($l_{p} / \\langle l \\rangle$)"
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
    ax1.set_xlim(1, 100)
    ax1.set_ylim(1, 100)
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
def plot_property_against_lp (df = None, R_col = None, Title = None, Y_label = None, X_label = None, dpi = None, logscale = False, lp_cos = False, lp_exp = False, xmin = None, ymin = None, xmax = None, ymax = None):

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
    plt.show()



## CLASSES
# none


## ARGUMENTS
# performing scaling analysis
scaling = ("scaling" in sys.argv)
# performe force extension analysis
forceExtension = ("forceExtension" in sys.argv)
# perform bending parameter analysis
bendingPARM = ("bendingPARM" in sys.argv)
# update the simulation results, even if the simulation results file already exists
update = ("update" in sys.argv)

## SCRIPT
# load defaults / settings from yaml

# parse arguments from command line

# generate directories, as instructed

# compile and analze results, as instructed
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
    if update or (not os.path.exists("02_processed_data/forceExtension/forceExtension.csv")):
        # get simulation parameters
        FE_parms = pd.read_csv(forceExtension_parmcsv)
        # calcuate bond vector difference
        check_bvd(parms = FE_parms, dir = '01_raw_data/forceExtension/', plot = True)
        # average simulation properties
        FE_parms = parse_results(parms = FE_parms, dir = '01_raw_data/forceExtension/', simfile = 'BVD.csv', col = 3, title = 'BVDMSE', M1 = True, M2 = True)
        FE_parms = parse_results(parms = FE_parms, dir = '01_raw_data/forceExtension/', simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
        FE_parms = parse_results(parms = FE_parms, dir = '01_raw_data/forceExtension/', simfile = 'ROG.dat', col = 4, title = 'ROGtot', M1 = True, M2 = True, var = True, bootstrapping = False, tabsep = True)
        # save results
        if not os.path.exists("02_processed_data/forceExtension/"):
            os.mkdir("02_processed_data/forceExtension/")
        FE_parms.to_csv("02_processed_data/forceExtension/forceExtension.csv")
    else:
        FE_parms = pd.read_csv("02_processed_data/forceExtension/forceExtension.csv")
    # perform analysis for each unqiue set of 'modules'
    # loop through each topology, bending potential, bending potential strength, and chain length
    for top in FE_parms['R'].unique():
        if top == 0:
            top_NAME = "CHAIN"
            # continue
        elif top == 1:
            top_NAME = "RING"
        else:
            print(f"forceExtension :: unknown TOP code {top}")
            exit()
        for pot in FE_parms['pot'].unique():
            for N in FE_parms['N'].unique():
                # create plot with force extension of bending constants merged
                top_df = FE_parms[(FE_parms['R'] == top) & (FE_parms['pot'] == pot) & (FE_parms['N'] == N)]
                if top_df.empty:
                    continue
                save_name = f"02_processed_data/forceExtension/FE_{top_NAME}_{pot}_N{N}"
                # plot force extension data for all bending constants, as well as log5 fit
                plot_force_extension(top_df, Y_col = 'E2Etot', X_col = 'F', iso_col = 'K', isolabel = '$k_{{\\theta}}$ = {:2}', X_label = "External Force ($f_{x}$)", Y_label = "Chain Extension (X-Direction)", Title = f"Force Extension Data for {top_NAME} with {pot} potential (N = {N})", saveas = save_name + '_data.png', plot_data = True, flip_axis = True)
                # plot error in bond vector distribution against force for all bending constants
                plot_force_extension(top_df, Y_col = 'BVDMSE', X_col = 'F', iso_col = 'K', isolabel = '$k_{{\\theta}}$ = {:2}', X_label = "External Force ($f_{x}$)", Y_label = "Bond Vector Distribution Mean Square Error", Title = f"Force Extension Data for {top_NAME} with {pot} potential (N = {N})", saveas = save_name + '_BVDMSE.png', plot_data = True, M2 = True)
                # create unique force extension plots for each unique bending constant
                for k in FE_parms['K'].unique():
                    # establish file names
                    save_name = f"02_processed_data/forceExtension/FE_{top_NAME}_{pot}_N{N}_k{k}"
                    # isolate the parameters corresponding to the set
                    isoK_df = FE_parms[(FE_parms['R'] == top) & (FE_parms['pot'] == pot) & (FE_parms['N'] == N) & (FE_parms['K'] == k)]
                    if isoK_df.empty:
                        print(save_name + " does not exist!")
                        continue
                    plot_scaling(isoK_df, N_col = 'F', R_col = ['E2Etot'], X_label = "External Force ($f_{x}$)", Y_label = "Chain Extension (X-Direction)", data_label = ["Simulation Data"], Title = f"Force Extension for {top_NAME} with {pot} potential ($k_{{\\theta}}$ = {k}, N = {N})", saveas = save_name + 'M2.png', error = True, color = ['tab:purple'], flip_axis = True)

if bendingPARM:
    if update or (not os.path.exists("02_processed_data/bendingPARM/bendingPARM.csv")):
        # load parameters, add property calculations from simulations
        bendingparms = pd.read_csv(bendingPARM_parmcsv)
        check_bvd(parms = bendingparms, dir = '01_raw_data/bendingPARM/')
        bendingparms = parse_results(parms = bendingparms, dir = '01_raw_data/bendingPARM/', simfile = 'RE2E.dat', col = 1, title = 'E2Ex', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
        bendingparms = parse_results(parms = bendingparms, dir = '01_raw_data/bendingPARM/', simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True)
        bendingparms = parse_results(parms = bendingparms, dir = '01_raw_data/bendingPARM/', simfile = 'BVD.csv', col = 3, title = 'BVDMSE', M1 = True, M2 = True)
        for i in range(30):
            bendingparms = parse_results(parms = bendingparms, dir = '01_raw_data/bendingPARM/', simfile = 'BBC.dat', row = i, title = "l"+str(i), M1 = True)
        # save results
        if not os.path.exists("02_processed_data/bendingPARM/"):
            os.mkdir("02_processed_data/bendingPARM/")
        bendingparms.to_csv("02_processed_data/bendingPARM/bendingPARM.csv")
    else:
        bendingparms = pd.read_csv("02_processed_data/bendingPARM/bendingPARM.csv")
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
            plot_bendingpot_bbc(df = plotdf, k = [1, 3, 5, 7, 10, 13], max_s = 29, plot_fit = True, Title = title_bbc, saveas = "02_processed_data/bendingPARM/" + saveas_bbc)
            # plot the bending parameter scaling against k
            title_lp = f"Presitance Length for {chain}s with {pot} Potential (N = {N_string})"
            saveas_lp = "lp_" + pot + "_" + chain + ".png"
            plot_bending_lp(df = plotdf, plot_expectation_CSA = (pot == "CSA"), plot_expectation_CA = (pot == "CA"), Title = title_lp, saveas = "02_processed_data/bendingPARM/" + saveas_lp, logscale = True, BVDMSE = True, show = True)
            # TODO plot the chain length against the persistence length
            plot_property_against_lp(df = plotdf, R_col = 'E2Etot_M1', lp_exp = True, logscale = True, Y_label = "End-to-End Distance", xmin = 1., xmax = 100., ymin = 10., ymax = 300.)
            # TODO plot the BVD MSE against the persistence length
            plot_property_against_lp(df = plotdf, R_col = 'BVDMSE_M2', lp_exp = True, logscale = True, Y_label = "Bond Vector Distribution Mean Square Difference", xmin = 1., xmax = 100., ymin = 0.001, ymax = .5)

