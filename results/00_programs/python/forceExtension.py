## PACKAGES
# conda packages
import sys, os, math
import csv
from pathlib import Path
import numpy as np
import pandas as pd
# local packages
from fig.Figure import Figure
from fig.fit import Line, bfm_langevin_fit, fit_line_no_slope, fit_line_no_intercept, fit_line
from plot.scatter_plot import gen_scatter
from plot.plot import gen_plot
from analysis import parse_results, plot_force_extension
from util.smoothie import numerical_slope, find_regime


## PARAMETERS
# default name used for force extension jobs
default_jobdir = "forceExtension"

## METHODS
# calculate parallel and perpendicular non-linear elasticity constants from bond order parameters
def plot_elasticity_bondop (path = "./", label = None, df = None, F_col = None, R_col = None, K_parallel_col = None, K_perp_col = None, iso_col = None, logx = True, logy = True, show = True, save = False, norm = None, linear = False, nonlinear = False):
	
	# check that the correct parameters were passed to the method
	if df is None:
		# must pass dataframe to the method
		print("analysis :: forceExtension :: plot_elasticity_bondop :: ERROR :: must pass dataframe to method as 'df'.")

	if F_col is None and R_col is None:
		# must specifcy column containing force
		print("analysis :: forceExtension :: plot_elasticity_bondop :: ERROR :: must specifiy column in 'df' containing force values as 'F_col'.")
	
	if K_parallel_col is None and K_perp_col is None:
		# must specify at least the parallel or perpendicular columns containing the bond order parameter
		print("analysis :: forceExtension :: plot_elasticity_bondop :: ERROR :: must specify either column containing parallel bond order parameter (as 'K_parallel_col') or column containing perpendicular bond order parameter (as 'K_parallel_col'), or both.")

	normalized = False
	if norm is not None:
		normalized = True
	if (not isinstance(norm, float)) and (not isinstance(norm, int)):
		# if norm is not a float or an int, reassign 1
		norm = 1.
		normalized = False

	# loop through the each line in the df, process the data
	f = [] # force associated with chain extension
	r = [] # chain extension
	k = [] # linear elastic constant from bond order paramemer fluctuation
	l = [] # label determines if elastic constant is either parallel or perpendicular
	for i, row in df.iterrows():

		if (F_col is not None and row[F_col] < 0.000001) or (R_col is not None and row[R_col] < 0.000001):
			## TODO :: opportunity to pull normalizing values here
			# skip if the minimum force or chain extension has not been surpased
			# (this is challenging for log log plots when there are negative values, etc.)
			continue

		# determine the parallel linear elastic constant
		if K_parallel_col is not None:
			if F_col is not None:
				f.append(row[F_col] * norm)
			if R_col is not None:
				r.append(row[R_col] / norm)
			k.append(1. / row[K_parallel_col + "_var"])
			l.append('parallel')

		# determine the perpendicular elastic constant
		if K_perp_col is not None:
			if F_col is not None:
				f.append(row[F_col] * norm)
			if R_col is not None:
				r.append(row[R_col] / norm)
			k.append(1. / row[K_perp_col + "_var"])
			l.append('perp')

	# plot the elasticity constants against the force
	if F_col:
		# if linear is true, identify the linear constant
		fit_linear_regime = None
		if linear:
			# accumulate the data for the parallel constant
			x = []
			y = []
			for i in range(len(f)):
				if l[i] == 'parallel':
					x.append(f[i])
					y.append(k[i])
			# find the regime corresponding to a slope equal to zero
			x_linear, m_linear = find_regime(x = x, y = y, slope = 0., tol = 0.15, log = True, monotonic = True, average_int = 10)
			# fit that data to a line without slope
			fit_linear_regime = fit_line_no_slope(x = x, y = y, x_start = x_linear[0], x_end = x_linear[-1])
			parms = fit_linear_regime.get_parameters()
			fit_linear_regime.set_linecolor("k")
			fit_linear_regime.set_linestyle(":")
			K_linear = parms[1]
			fit_linear_regime.set_label(f"Linear Regime\n(K = {K_linear:.04f})")

		fit_nonlinear_regime = None
		if nonlinear:
			# accumulate the data for the parallel constant
			x = []
			y = []
			for i in range(len(f)):
				if l[i] == 'parallel':
					x.append(f[i])
					y.append(k[i])
			# find the regime corresponding to a slope equal to two
			x_nonlinear, m_nonlinear = find_regime(x = x, y = y, slope = 2.1, tol = 0.5, log = True, monotonic = True, average_int = 10, x_start = x_linear[-1])
			for i in range(len(x_nonlinear)):
				print(i, x_nonlinear[i], m_nonlinear[i])
			# fit that data to a line without slope
			fit_nonlinear_regime = fit_line(x = x, y = y, x_start = x_nonlinear[0], x_end = x_nonlinear[-1], log = True)
			parms = fit_nonlinear_regime.get_parameters()
			fit_nonlinear_regime.set_linecolor("k")
			fit_nonlinear_regime.set_linestyle("--")
			K_nonlinear = parms[0]
			fit_nonlinear_regime.set_label(f"Non-Linear Regime\n($K \\propto F^{{{K_nonlinear:.02f}}}$)")

		df = pd.DataFrame.from_dict({'x': f, 'y': k, 'i': l})
		fig = Figure()
		fig.load_data(d = df, xcol = 'x', ycol = 'y', icol = 'i')
		if normalized:
			fig.set_xaxis_label(f"Normalized Force ($F \\cdot X$ where (X = {norm:.2f}))")
			fig.set_title_label("Non-Linear Elasticity against Normalized Force")
		else:
			fig.set_xaxis_label("Force ($F$)")
			fig.set_title_label("Non-Linear Elasicity against Force")
		fig.set_yaxis_label("Non-Linear Elastic Constant ($K$)")
		fig.set_subtitle_label(label)
		fig.set_label(ival = 'parallel', label = "Parallel ($K_{\\parallel}$)")
		fig.set_label(ival = 'perp', label = "Perpendicular ($K_{\\perp}$)")
		fig.set_yaxis_min(1.)
		if logx is True or logy is True:
			fig.set_logscale(logx = logx, logy = logy)
		if save:
			# set save location and name
			fig.set_saveas(savedir = path, filename = label)
			# save data
			fig.save_data()

		# plot
		gen_plot(fig = fig, show = show, save = save, fit = [fit_linear_regime, fit_nonlinear_regime])

	# plot the elasticity constants against the chain extension
	if R_col:
		df = pd.DataFrame.from_dict({'x': r, 'y': k, 'i': l})
		fig = Figure()
		fig.load_data(d = df, xcol = 'x', ycol = 'y', icol = 'i')
		if normalized:
			fig.set_xaxis_label(f"Normalized Chain Extension ($R_{{E2E}} / X$ where (X = {norm:.2f}))")
			fig.set_title_label("Non-Linear Elasticity against Normalized Chain Extension")
		else:
			fig.set_xaxis_label("Chain Extension ($R_{E2E}$)")
			fig.set_title_label("Non-Linear Elasticity againt Chain Extension")
		fig.set_yaxis_label("Linear Elastic Constant ($K$)")
		fig.set_subtitle_label(label)
		fig.set_label(ival = 'parallel', label = "Parallel ($K_{\\parallel}$)")
		fig.set_label(ival = 'perp', label = "Perpendicular ($K_{\\perp}$)")
		fig.set_yaxis_min(1.)
		if logx is True or logy is True:
			fig.set_logscale(logx = logx, logy = logy)
		if save:
			# set save location and name
			fig.set_saveas(savedir = path, filename = label)
			# save data
			fig.save_data()

		# plot
		gen_plot(fig = fig, show = show, save = save)


# calculate elasticity numerically from force extension curve
def calc_elasticity_numerically (df = None, fcol = None, rcol = None, icol = None, ilabel = None, plot = False, show = False, norm = None, monotonic = False, spline = False, average_int = 5, savedir = None, real = False, ideal = False):

	## CHECK ARGUMENTS
	# check that the data frame was provided
	if df is None:
		# must pass dataframe to the method
		print("analysis :: forceExtension :: calc_elasticity_numerically :: ERROR :: must pass dataframe to method as 'df'.")
		return

	# check that the force column was specified
	if fcol is None:
		# must provide the label for the force column
		print("analysis :: forceExtension :: calc_elasticity_numerically :: ERROR :: must specify the column in 'df' which contains the force data as 'fcol'.")
		return

	# check that the chain extension column was specified
	if rcol is None:
		# must provide the label for the chain extension column
		print("analysis :: forceExtension :: calc_elasticity_numerically :: ERROR :: must specify the column in 'df' which contains the chain extensino data as 'rcol'.")
		return 

	# if smoothing factors were specified, check that one spline or monotonic is called
	if spline and monotonic:
		spline = False
	
	# if normalization factor was specified
	normalized = False
	if norm is not None:
		normalized = True
	if (not isinstance(norm, float)) and (not isinstance(norm, int)):
		# if norm is not a float or an int, reassign 1
		norm = 1.
		normalized = False

	## PARSE AND ANALYZE DATA
	# parse data
	f = df[fcol].to_list()
	r = df[rcol].to_list()

	# normalize if requested
	if normalized:
		f = f * norm
		r = r / norm

	# calculate the slope, using smoothing if specified
	r0, dfdr = numerical_slope(x = r, y = f, log = False, average_int = average_int, monotonic = monotonic, spline = spline)
	f0, drdf = numerical_slope(x = f, y = r, log = False, average_int = average_int, monotonic = monotonic, spline = spline)

	# remove and negative numbers
	for i in range(len(r0) - 1, -1, -1):
		if dfdr[i] < 0.:
			dfdr.pop(i)
			r0.pop(i)
			drdf.pop(i)
			f0.pop(i)

	## DETERMING THE LINEAR REGIME AND LINEAR ELASTICITY CONSTANT
	x_kvf_linear, m_kvf_linear = find_regime (x = f0, y = dfdr, slope = 0., tol = 0.3, log = True, monotonic = True, average_int = 10)
	# create fits for the linear and nonlinear regimes
	fit_kvf_linear = None
	if x_kvf_linear is not None:
		for i in range(len(x_kvf_linear)):
			print(i, x_kvf_linear[i], m_kvf_linear[i])
		fit_kvf_linear = fit_line_no_slope(x = f0, y = dfdr, x_start = x_kvf_linear[0], x_end = x_kvf_linear[-1])
		parms = fit_kvf_linear.get_parameters()
		fit_kvf_linear.set_linecolor("k")
		fit_kvf_linear.set_linestyle(":")
		K_linear = parms[1]
		fit_kvf_linear.set_label(f"Linear Regime\n($K = {K_linear:.04f}$)")
	else:
		K_linear = None

	# identify the pincus regime for real chains
	fit_kvf_nonlinear = None
	fit_kvr_nonlinear = None
	if real and K_linear is not None:
		# the elastic constant vs. the force should scale with approximately 1/3 in the pinucs regime
		x_kvf_nonlinear, __ = find_regime (x = f0, y = dfdr, slope = (2./3.), tol = 0.3, log = True, monotonic = True, average_int = 10, count_threshold = 5, x_start = x_kvf_linear[-1])

		# fit the nonlinear data to a line against the force
		if x_kvf_nonlinear is not None:
			fit_kvf_nonlinear = fit_line(x = f0, y = dfdr, x_start = x_kvf_nonlinear[0], x_end = x_kvf_nonlinear[-1], log = True)
			parms = fit_kvf_nonlinear.get_parameters()
			fit_kvf_nonlinear.set_linecolor("k")
			fit_kvf_nonlinear.set_linestyle("--")
			kvf_nonlinear = parms[0]
			fit_kvf_nonlinear.set_label(f"Non-Linear Regime\n($K \\propto F^{{{kvf_nonlinear:.04f}}}$)")

		# the elastic constant vs. the chain length should scale with approximately 1/2 in the pincus regime
		x_kvr_nonlinear, __ = find_regime (x = r0, y = dfdr, slope = (1.), tol = 0.3, log = True, monotonic = True, average_int = 10, count_threshold = 5, x_start = x_kvf_linear[-1])

		# fit the nonlinear data to a line against the chain length
		if x_kvr_nonlinear is not None:
			fit_kvr_nonlinear = fit_line(x = r0, y = dfdr, x_start = x_kvr_nonlinear[0], x_end = x_kvr_nonlinear[-1], log = True)
			parms = fit_kvr_nonlinear.get_parameters()
			fit_kvr_nonlinear.set_linecolor("k")
			fit_kvr_nonlinear.set_linestyle("--")
			kvr_nonlinear = parms[0]
			fit_kvr_nonlinear.set_label(f"Non-Linear Regime\n($K \\propto R^{{{kvr_nonlinear:.04f}}}$)")

	# for ideal chains, plot the derivative of the bond-fluctuation langevin function
	label = ["Simulation Data"] * len(f0)
	if ideal:
		# create the bfm fit
		fit = bfm_langevin_fit()
		bfm_f = fit.get_xval_list(lims = (min(f0), max(f0)), n = 100, log = True)
		bfm_r = fit.get_yval_list(lims = (min(f0), max(f0)), n = 100, log = True)
		bfm_r0, bfm_dfdr = numerical_slope(x = bfm_r, y = bfm_f, average_int = average_int, monotonic = monotonic, spline = spline, log = False)
		bfm_f0, __ = numerical_slope(x = bfm_f, y = bfm_r, average_int = average_int, monotonic = monotonic, spline = spline, log = False)
		# append the values to the list
		for i in range(len(bfm_f0)):
			f0.append(bfm_f0[i])
			r0.append(bfm_r0[i])
			dfdr.append(bfm_dfdr[i])
			label.append("BFM Langevin Function")

	## PLOT EvF 
	fig = Figure()
	fig.load_data(d = pd.DataFrame.from_dict({'F':f0, 'dfdr': dfdr, 'l': label}), xcol = 'F', ycol = 'dfdr', icol = 'l')
	if normalized:
		fig.set_title_label("Normalized Elastic Constant against Normalized Force")
		fig.set_xaxis_label("Normalized Force ($F \\cdot X$)")
		fig.set_yaxis_label("Normalized Elastic Constant ($K = d(F \\cdot X) / d(R \\cdot X^{{-1}})$)")
	else:
		fig.set_title_label("Elastic Constant against Force")
		fig.set_xaxis_label("Force ($F$)")
		fig.set_yaxis_label("Elastic Constant ($K = R_{{max}} \\cdot d(F) / d(R)$)")
	fig.set_logscale()
	if ideal:
		fig.set_marker("BFM Langevin Function", "x")
	if savedir is not None:
		fig.set_saveas(savedir = savedir, filename = "EvF")
		fig.save_data()
	gen_plot(fig = fig, markersize = 2, linewidth = 1, show = False, save = (savedir is not None), fit = [fit_kvf_linear, fit_kvf_nonlinear], legendloc = 'upper left')

	# PLOT EvR
	fig = Figure()
	fig.load_data(d = pd.DataFrame.from_dict({'R':r0, 'dfdr': dfdr, 'l': label}), xcol = 'R', ycol = 'dfdr', icol = 'l')
	if normalized:
		fig.set_title_label("Normalized Elastic Constant against Normalized Chain Extension")
		fig.set_xaxis_label("Normalized Chain Extension ($R / X$)")
		fig.set_yaxis_label("Normalized Elastic Constant ($K = d(F \\cdot X) / d(R \\cdot X^{{-1}})$)")
	else:
		fig.set_title_label("Elastic Constant against Chain Extension")
		fig.set_xaxis_label("Normalized Chain Extension ($R / R_{{max}})$)")
		fig.set_yaxis_label("Elastic Constant ($K = R_{{max}} \\cdot d(F) / d(R)$)")
	fig.set_logscale()
	if savedir is not None:
		fig.set_saveas(savedir = savedir, filename = "EvR")
		fig.save_data()
	gen_plot(fig = fig, markersize = 2, linewidth = 1, show = False, save = (savedir is not None), fit = [fit_kvf_linear, fit_kvr_nonlinear], legendloc = 'upper left')

	return K_linear


def plot_force_extension (df = None, fcol = None, rcol = None, icol = None, ilabel = None, show = False, savedir = None, real = False, ideal = False, rmax = None, scale = None):

	# initialize fits
	fit_langevin = None
	fit_linear = None
	fit_pincus = None

	# loop through each icol, parse data, determine regimes
	f = []
	r = []
	l = []
	## TODO :: think there is a better way to run this without so much looping, but I am not entirely sure how yet
	if icol is not None:
		for i in df[icol].unique():
			# get the portion of the data frame corresponding to the unique value
			udf = df[df[icol] == i]
			# get the no chain extension data, drop from udf and reset
			r_nf = udf[udf[fcol] <= 0.]
			x_nf = r_nf.reset_index().loc[0, rcol] # no force chain extension, used for scaling
			for i in r_nf.index:
				udf = udf.drop(i).dropna()

			# accumulate the chain extension data
			f_temp = udf[fcol].to_list()
			r_temp = udf[rcol].to_list()
			l_temp = udf[icol].to_list()
			for i in range(len(f_temp)):
				# append the force
				if scale:
					f_temp[i] = f_temp[i] * x_nf
					f.append(f_temp[i])
				else:
					f.append(f_temp[i])
				# append the chain extension
				if scale:
					r_temp[i] = r_temp[i] / x_nf
					r.append(r_temp[i])
				elif rmax is not None:
					r_temp[i] = r_temp[i] / rmax
					r.append(r_temp[i])
				else:
					r.append(r_temp[i])
				# append the label 
				scale_label = " "
				if scale:
					scale_label += "\n(X = {:.02f})".format(x_nf)
				if ilabel is None or not isinstance(ilabel, str):
					l.append(l_temp[i] + scale_label)
				else:
					label = ilabel.format(l_temp[i])
					l.append(label + scale_label)

			# fit langevin data
			if fit_langevin is None:
				if ideal:
					fit_langevin = bfm_langevin_fit()

			# fit the linear regime
			if fit_linear is None:
				if real:
					x_linear, m_linear = find_regime (x = f_temp, y = r_temp, slope = 1., tol = 0.05, log = True, monotonic = True, average_int = 10)
					if x_linear is not None: # if the regime was found
						fit_linear = fit_line_no_intercept(x = f_temp, y = r_temp, x_start = x_linear[0], x_end = x_linear[-1])
						fit_linear.set_linecolor("k")
						fit_linear.set_linestyle("--")
						fit_linear.set_label(f"Linear Regime\n($R \\propto F$)")

			# fit the pincus regime
			if fit_pincus is None:
				if real:
					# for real chains, the pincus regime has slope of ~2/3
					x_pincus, m_pincus = find_regime (x = f_temp, y = r_temp, slope = (2. / 3.), tol = 0.05, log = True, monotonic = True, average_int = 10)
					if x_pincus is not None:
						fit_pincus = fit_line(x = f_temp, y = r_temp, x_start = x_pincus[0], x_end = x_pincus[-1], log = True)
						p = fit_pincus.get_parameters()
						fit_pincus.set_linestyle(":")
						fit_pincus.set_linecolor("k")
						fit_pincus.set_label(f"Pincus Regime \n ($R \\propto F^{{{p[0]:.03f}}}$)")
	else:
		print(" forceExtension :: plot_force_extension :: ERROR :: must specify 'icol'.")

	# plot the force extension on a regular scale
	fig = Figure()
	fig.load_data(d = pd.DataFrame.from_dict({'F': f, 'R': r, 'L': l}), xcol = 'F', ycol = 'R', icol = 'L')
	# set y axis label
	if scale:
		fig.set_yaxis_label(f"Chain Extension ($R_{{E2E}} \\cdot X^{{-1}}$)")
	elif rmax is not None:
		fig.set_yaxis_label(f"Normalized Chain Extension ($R_{{E2E}} / R_{{max}}$)")
	else:
		fig.set_yaxis_label(f"Chain Extension ($R_{{E2E}}$)")
	# set x axis label
	if scale:
		fig.set_xaxis_label(f"Normalized Chain Extension ($F \\cdot X$)")
	else:
		fig.set_xaxis_label(f"Force ($F$)")
	# set title
	fig.set_title_label(f"Force Extension Curve")
	# fig.set_subtitle_label(d)
	fig.set_xaxis_min(0.0)

	if savedir is not None:
		if scale:
			fig.set_saveas(savedir = savedir, filename = "RvF_scale")
		else:
			fig.set_saveas(savedir = savedir, filename = "RvF")
		fig.save_data()
	gen_plot(fig = fig, show = show, fit = fit_langevin, legendloc = 'lower right')

	# plot force extension on a log-log scale
	fig.set_logscale()
	fig.set_yaxis_min(min(r))
	fig.set_xaxis_min(min(f))
	if savedir is not None:
		if scale:
			fig.set_saveas(savedir = savedir, filename = "RvF_scale_log")
		else:
			fig.set_saveas(savedir = savedir, filename = "RvF_log")
	gen_plot(fig = fig, show = show, fit = [fit_linear, fit_pincus], legendloc = 'lower right')



## ARGUMENTS
# update the simulation results, even if the simulation results file already exists
update = ("update" in sys.argv)
# create all figures
figure = ("allfigs" in sys.argv)
# generate / re-generate figure 1
fig1 = ("fig1" in sys.argv)
# show figures on generation 
show = ("show" in sys.argv)
# analysis of real polymer chains
real = ("real" in sys.argv)
# analysis of ideal polymer chains
ideal = ("ideal" in sys.argv)
# parse the directory passed as an argument, if specified
if ("job" in sys.argv):
	jobdir = sys.argv[sys.argv.index("job") + 1]
else:
	jobdir = default_jobdir

## SCRIPT
# collect raw simulation data
lin = ['log']
for l in lin:
	# force extension for linear data set
	job = jobdir + '_' + l
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
		# FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True, plot = False)
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 5, title = 'E2Eproj', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True, plot = False)
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 8, title = 'cos_theta_parallel', M1 = True, M2 = True, var = True, bootstrapping = False, tabsep = True, plot = False)
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 9, title = 'cos_theta_perp', M1 = True, M2 = True, var = True, bootstrapping = False, tabsep = True, plot = False)
		# save results
		if not os.path.exists(anal_dir):
			os.makedirs(anal_dir)
		FE_parms.to_csv(anal_dir + parm_file, index = False)
	else:
		FE_parms = pd.read_csv(anal_dir + parm_file, index_col = False)

	## FORCE-EXTENSION, ELASTICITY ANALYSIS
	for n in FE_parms['N'].unique(): # chain length
		for top in FE_parms['top'].unique(): # unique topology

			## COMPARE EACH FORCE EXTENSION CURVE FOR EACH TOPOLOGY
			# get the unique dataframe, establish save name and directory
			save = f"{top}_N{n:03d}"
			savedir = f"./02_processed_data/{job:s}/{save:s}/"
			df = FE_parms[(FE_parms['top'] == top) & (FE_parms['N'] == n)]

			## PLOT FORCE-EXTENSION
			# regular
			plot_force_extension(df = df, fcol = 'F', rcol = 'E2Eproj_M1', icol = 'K', ilabel = "$\\kappa_{{\\theta}}$ = {:.02f}", real = real, ideal = ideal, show = False, savedir = savedir, rmax = (3 * n))
			# add scaling
			plot_force_extension(df = df, fcol = 'F', rcol = 'E2Eproj_M1', icol = 'K', ilabel = "$\\kappa_{{\\theta}}$ = {:.02f}", real = real, ideal = ideal, show = False, savedir = savedir, scale = True)
			exit()


			for k in FE_parms['K'].unique(): # strength of bending potential
				## ANALYZE EACH UNIQUE FORCE-EXTENSION SIMULATION
				# for each, create a unique directory
				d = f"{top}_N{n:03d}_k{k:0.2f}"
				df = FE_parms[(FE_parms['top'] == top) & (FE_parms['K'] == k) & (FE_parms['N'] == n)]
				df = df.reset_index()
				print(d)

				# get the no force change extension and then remove it from the dataframe
				r_nf = df.loc[0, 'E2Eproj_M1']
				df = df.reset_index().drop(0).dropna()
				# get the max chain extension
				r_max = df.loc[(len(df.index - 1)), 'E2Eproj_M1']

				## PLOT THE FORCE EXTENSION CURVE (RvF), compare to Langevin model
				df_norm = df.copy(deep = True)
				df_norm['F'] = df['F'] #* 2.69
				df_norm['E2Eproj_M1'] = df_norm['E2Eproj_M1'] / r_max
				df_norm['label'] = ["Simulation Data"] * len(df['F'].to_list())
				plot_force_extension (df = df_norm, rcol = 'E2Eproj_M1', fcol = 'F', real = real, ideal = ideal, show = False, savedir = f"./02_processed_data/{job:s}/{d:s}/")

				## CALCULATE THE ELASTICITY NUMERICALLY, determine linear constants
				k_lin = calc_elasticity_numerically(df = df_norm, F_col = 'F', R_col = 'E2Eproj_M1', monotonic = False, average_int = 5, savedir = f"./02_processed_data/{job:s}/{d:s}/", real = real, ideal = ideal, icol = 'label')

				## PARALLEL AND PERPENDICULAR ELASTICITY
				# plot the parallel and perpendicular elasticity against the applied force
				# determine the slope of the linear and non-linear regimes
				plot_elasticity_bondop(df = df_norm, F_col = 'F', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_force", linear = True, nonlinear = True)
				plot_elasticity_bondop(df = df_norm, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension", linear = True, nonlinear = True, logx = False)
				continue
				# same plot with normalization
				plot_elasticity_bondop(df = df, norm = float(r_nf), F_col = 'F', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_force_norm")
				## TODO get the transition force

				# plot the parallel and perpendicular elasticity against the chain extension
				plot_elasticity_bondop(df = df, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension")
				# same plot but with normalization
				plot_elasticity_bondop(df = df, norm = r_nf, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension_norm")
				# plot the nonlinear elasticity on a linear scale against the unitary chain extension
				plot_elasticity_bondop(df = df, norm = r_max, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension_max", logx = False)
		

	## ADD HEAT MAP OF KLIN PER TOPOLOGY AND BENDING PARM
	## ADD SCATTER PLOT OF TRANSITION LENGTH OF FORCE PER BENDING PARM AND TOPOLOGY

	# force extension curve comparing topologies with constant elasticity
	if not os.path.exists(anal_dir + "FE_top/"):
		os.makedirs(anal_dir + "FE_top/")
	for k in FE_parms['K'].unique():
		for pot in FE_parms['pot'].unique():
			for N in FE_parms['N'].unique():
				# get the unique df
				k_df = FE_parms[(FE_parms['K'] == k) & (FE_parms['pot'] == pot) & (FE_parms['N'] == N)]
				save_name = anal_dir + "FE_top/" + f"FE_k{k:0.2f}_{pot}_N{N}"
				# master force extension plot
				plot_force_extension (k_df, Y_col = 'E2Eproj', X_col = 'F', iso_col = 'top', isolabel = '{:s}', X_label = "Normalized External Force ($X \\cdot f$)", Y_label = "Normalized Chain Extension ($X^{{-1}} \\cdot R_{{E2E}}$)", saveas = save_name + '_norm_data', plot_data = True, y_max = 10., y_min = 0.001, x_min = 0.01, x_max = 100., show = show, logscale_x = False, logscale_y = False, plot_slope = False, norm = True, pincus = True, hookean = True) # saveas = save_name + '_norm_data.png',
				# plot without normalization
				plot_force_extension (k_df, Y_col = 'E2Eproj', X_col = 'F', iso_col = 'top', isolabel = '{:s}', X_label = "External Force ($f$)", Y_label = "Chain Extension ($R_{{E2E}}$)", saveas = save_name + '_data', plot_data = True, y_max = 500., y_min = 0.1, x_min = 0.0001, x_max = 5., show = show, logscale_x = True, logscale_y = True, plot_slope = False)

	# force extension curve comparing elasticity with constant topology
	if not os.path.exists(anal_dir + "FE_bend/"):
		os.makedirs(anal_dir + "FE_bend/")
	for top in FE_parms['top'].unique():
		for pot in FE_parms['pot'].unique():
			for N in FE_parms['N'].unique():
				top_df = FE_parms[(FE_parms['top'] == top) & (FE_parms['pot'] == pot) & (FE_parms['N'] == N)]
				if top_df.empty:
					# if the data frame is empty, skip
					continue
				save_name = anal_dir + "FE_bend/" + f"FE_{top}_{pot}_N{N}"
				# normalized master force extension plot
				plot_force_extension (top_df, Y_col = 'E2Eproj', X_col = 'F', iso_col = 'K', isolabel = '$\\kappa_{{\\theta}}$ = {:.02f}', X_label = "Normalized External Force ($X \\cdot f$)", Y_label = "Normalized Chain Extension ($X^{{-1}} \\cdot R_{{E2E}}$)", saveas = save_name + '_norm_data', plot_data = True, y_max = 10., y_min = 0.001, x_min = 0.01, x_max = 100., show = show, logscale_x = True, logscale_y = True, plot_slope = False, norm = True, pincus = True, hookean = True) # saveas = save_name + '_norm_data.png',
				# plot without normalization
				plot_force_extension (top_df, Y_col = 'E2Eproj', X_col = 'F', iso_col = 'K', isolabel = '$k_{{\\theta}}$ = {:.02f}', X_label = "External Force ($f$)", Y_label = "Chain Extension ($R_{{E2E}}$)", saveas = save_name + '_data', plot_data = True, y_max = 500., y_min = 0.1, x_min = 0.0001, x_max = 5., show = show, logscale_x = True, logscale_y = True, plot_slope = False)



	# parallel and perpendicular elasticity
	# elasticity from force-extension curve
