## PACKAGES
# conda packages
import sys, os, math
import csv
from pathlib import Path
import numpy as np
import pandas as pd
# local packages
from fig.Figure import Figure
from fig.scatter_plot import gen_scatter
from analysis import parse_results, plot_force_extension
from util.smoothie import monotonic_slope
from fig.scatter_plot import gen_scatter


## PARAMETERS
# none

## METHODS
# calculate parallel and perpendicular non-linear elasticity constants from bond order parameters
def plot_elasticity_bondop (path = "./", label = None, df = None, F_col = None, R_col = None, K_parallel_col = None, K_perp_col = None, iso_col = None, log = True, show = True, save = False):
	
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
				f.append(row[F_col])
			if R_col is not None:
				r.append(row[R_col])
			k.append(1. / row[K_parallel_col + "_var"])
			l.append('parallel')

		# determine the perpendicular elastic constant
		if K_perp_col is not None:
			if F_col is not None:
				f.append(row[F_col])
			if R_col is not None:
				r.append(row[R_col])
			k.append(1. / row[K_perp_col + "_var"])
			l.append('perp')

	# plot the elasticity constants against the force
	if F_col:
		df = pd.DataFrame.from_dict({'x': f, 'y': k, 'i': l})
		fig = Figure()
		fig.load_data(d = df, xcol = 'x', ycol = 'y', icol = 'i')
		fig.set_xaxis_label("Force ($F$)")
		fig.set_yaxis_label("Linear Elastic Constant ($K$)")
		fig.set_title_label("Non-Linear Elasicity against Force")
		fig.set_subtitle_label(label)
		fig.set_label(ival = 'parallel', label = "Parallel ($K_{\\parallel}$)")
		fig.set_label(ival = 'perp', label = "Perpendicular ($K_{\\parallel}$)")
		if log is True:
			fig.set_logscale()
		if save:
			if label is None:
				fig.set_saveas(savedir = ".", filename = f"elasticity_v_force")
			else:
				fig.set_saveas(savedir = f"{path:s}{label:s}/", filename = f"elasticity_v_force")
			# save data
			fig.save_data()

		# plot
		gen_scatter(fig = fig, show = show, save = save, markersize = 36)

	# plot the elasticity constants against the chain extension
	if R_col:
		df = pd.DataFrame.from_dict({'x': r, 'y': k, 'i': l})
		fig = Figure()
		fig.load_data(d = df, xcol = 'x', ycol = 'y', icol = 'i')
		fig.set_xaxis_label("Chain Extension ($R_{E2E}$)")
		fig.set_yaxis_label("Linear Elastic Constat ($K$)")
		fig.set_title_label("Non-Linear Elasticity againt Chain Extension")
		fig.set_subtitle_label(label)
		fig.set_label(ival = 'parallel', label = "Parallel ($K_{\\parallel}$)")
		fig.set_label(ival = 'perp', label = "Perpendicular ($K_{\\parallel}$)")
		if log is True:
			fig.set_logscale()
		if save:
			if label is None:
				fig.set_saveas(savedir = ".", filename = f"elasticity_v_extension")
			else:
				fig.set_saveas(savedir = f"{path:s}{label:s}/", filename = f"elasticity_v_extension")
			# save data
			fig.save_data()

		# plot
		gen_scatter(fig = fig, show = show, save = save, markersize = 36)

# calculate the linear elastic constant from the force extension curve
def calc_linear_elastic_constant(df = None, F_col = None, R_col = None, plot = False, monotonic = True):

	## check that the correct arguments were passed to the method
	# check that the data frame was provided
	if df is None:
		# must pass dataframe to the method
		print("analysis :: forceExtension :: calc_linear_elastic_constant :: ERROR :: must pass dataframe to method as 'df'.")
		return

	# check that the force column was specified
	if F_col is None:
		# must provide the label for the force column
		print("analysis :: forceExtension :: calc_linear_elastic_constant :: ERROR :: must specify the column in 'df' which contains the force data as 'F_col'.")
		return

	if R_col is None:
		# must provide the label for the chain extension column
		print("analysis :: forceExtension :: calc_linear_elastic_constant :: ERROR :: must specify the column in 'df' which contains the chain extensino data as 'R_col'.")
		return 

	## determine the region where the relationship between the chain extension and the force is linear
	# parse data
	y = df[F_col]
	x = df[R_col]
	# calculate the slope
	if monotonic:
		# calculate the slope using a monotonic (smoothed function)
		xder, yder = monotonic_slope(x = x, y = y)
		for i in range(len(xder)):
			print(x[i], y[i], xder[i], yder[i])
		exit()
	# calculate the slope using a spline function
	# calculate the slope using a standard derivative
	# determine the region where the slope is constant (within the variance)

	## plot the force against the chain extension with the linear region identified

# calculate elasticity numerically from force extension curve
def plot_elasticity_numerically ():
	pass

## ARGUMENTS
# update the simulation results, even if the simulation results file already exists
update = ("update" in sys.argv)
# create all figures
figure = ("allfigs" in sys.argv)
# generate / re-generate figure 1
fig1 = ("fig1" in sys.argv)
# show figures on generation 
show = ("show" in sys.argv)

## SCRIPT
# collect raw simulation data
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
		# FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True, plot = False)
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 5, title = 'E2Eproj', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True, plot = False)
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 8, title = 'cos_theta_parallel', M1 = True, M2 = True, var = True, bootstrapping = False, tabsep = True, plot = False)
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 9, title = 'cos_theta_perp', M1 = True, M2 = True, var = True, bootstrapping = False, tabsep = True, plot = False)
		# save results
		if not os.path.exists(anal_dir):
			os.makedirs(anal_dir)
		FE_parms.to_csv(anal_dir + parm_file)
	else:
		FE_parms = pd.read_csv(anal_dir + parm_file)

	# process data 
	# estimate the linear elasticity constant
	for top in FE_parms['top'].unique(): # unique topology
		for k in FE_parms['K'].unique(): # strength of bending potential
			# for each, create a unique directory
			d = f"{top}_k{k:0.2f}"
			## PARALLEL AND PERPENDICULAR ELASTICITY
			# plot the parallel and perpendicular elasticity against the applied force
			# plot the parallel and perpendicular elasticity against the chain extension
			plot_elasticity_bondop(df = FE_parms[(FE_parms['top'] == top) & (FE_parms['K'] == k)], F_col = 'F', R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/forceExtension_{l:s}/", label = d)
			# plot the force extension curve along with the numerically calcualted elasticity
			# determine the linear elastic constant
			k_lin = calc_linear_elastic_constant(df= FE_parms[(FE_parms['top'] == top) & (FE_parms['K'] == k)], F_col = 'F', R_col = 'E2Eproj_M1')

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
