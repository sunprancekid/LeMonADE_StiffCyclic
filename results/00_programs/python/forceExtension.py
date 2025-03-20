## PACKAGES
# conda packages
import sys, os, math
import csv
from pathlib import Path
import numpy as np
import pandas as pd
# local packages
from fig.Figure import Figure
from fig.fit import Line
from plot.scatter_plot import gen_scatter
from plot.plot import gen_plot
from analysis import parse_results, plot_force_extension
from util.smoothie import numerical_slope


## PARAMETERS
# default name used for force extension jobs
default_jobdir = "forceExtension"

## METHODS
# calculate parallel and perpendicular non-linear elasticity constants from bond order parameters
def plot_elasticity_bondop (path = "./", label = None, df = None, F_col = None, R_col = None, K_parallel_col = None, K_perp_col = None, iso_col = None, logx = True, logy = True, show = True, save = False, norm = None):
	
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
		df = pd.DataFrame.from_dict({'x': f, 'y': k, 'i': l})
		fig = Figure()
		fig.load_data(d = df, xcol = 'x', ycol = 'y', icol = 'i')
		if normalized:
			fig.set_xaxis_label(f"Normalized Force ($F \\cdot X$ where (X = {norm:.2f}))")
			fig.set_title_label("Non-Linear Elasticity against Normalized Force")
		else:
			fig.set_xaxis_label("Force ($F$)")
			fig.set_title_label("Non-Linear Elasicity against Force")
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
		gen_scatter(fig = fig, show = show, save = save, markersize = 36)

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
		gen_scatter(fig = fig, show = show, save = save, markersize = 36)

# calculate the linear elastic constant from the force extension curve
def calc_linear_elastic_constant(df = None, F_col = None, R_col = None, plot = False, monotonic = True, smooth_int = 10):

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
	y = df[F_col].to_list()
	x = df[R_col].to_list()
	# calculate the slope
	if monotonic:
		# calculate the slope using a monotonic (smoothed function)
		xder, yder = numerical_slope(x = x, y = y, log = True, average_int = smooth_int)
		for i in range(len(xder)):
			print(f"{x[i]:.04f}\t{y[i]:.04f}\t{xder[i]:.05f}\t{yder[i]:.05f}")
		# exit()
	# calculate the slope using a spline function
	# calculate the slope using a standard derivative
	# determine the region where the slope is constant (within the variance)

	## plot the force against the chain extension with the linear region identified

# calculate elasticity numerically from force extension curve
def calc_elasticity_numerically (df = None, F_col = None, R_col = None, plot = False, show = False, norm = None, monotonic = False, spline = False):

	## CHECK ARGUMENTS
	# check that the data frame was provided
	if df is None:
		# must pass dataframe to the method
		print("analysis :: forceExtension :: calc_elasticity_numerically :: ERROR :: must pass dataframe to method as 'df'.")
		return

	# check that the force column was specified
	if F_col is None:
		# must provide the label for the force column
		print("analysis :: forceExtension :: calc_elasticity_numerically :: ERROR :: must specify the column in 'df' which contains the force data as 'F_col'.")
		return

	# check that the chain extension column was specified
	if R_col is None:
		# must provide the label for the chain extension column
		print("analysis :: forceExtension :: calc_elasticity_numerically :: ERROR :: must specify the column in 'df' which contains the chain extensino data as 'R_col'.")
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
	f = df[F_col].to_list()
	r = df[R_col].to_list()

	# normalize if requested
	if normalized:
		f = f * norm
		r = r / norm

	# calculate the slope, using smoothing if specified
	r0, dfdr = numerical_slope(x = r, y = f, log = False, average_int = 5, monotonic = monotonic, spline = spline)
	f0, drdf = numerical_slope(x = f, y = r, log = False, average_int = 5, monotonic = monotonic, spline = spline)

	# determine the linear region and the linear elastic constant
	# the slope with the line regime is constant and zero
	# find the regine with slope close to zero
	# calculat the slope of the slope
	f1, fdrdf = numerical_slope(x = f0, y = dfdr, log = False, average_int = 5, monotonic = True)
	for i in range(len(f1)):
		print(f"{f1[i]:.04f}\t{fdrdf[i]:.04f}")
	exit()

	# TODO determine the pincus regime and the non-linear elastic pincus constant
	# plot the elastic constant
	fig = Figure()
	fig.load_data(d = pd.DataFrame.from_dict({'F':f0, 'dfdr': dfdr}), xcol = 'F', ycol = 'dfdr')
	if normalized:
		fig.set_title_label("Normalized Elastic Constant against Normalized Force")
		fig.set_xaxis_label("Normalized Force ($F \\cdot X$)")
		fig.set_yaxis_label("Normalized Elastic Constant ($K = d(F \\cdot X) / d(R \\cdot X^{{-1}})$)")
	else:
		fig.set_title_label("Elastic Constant against Force")
		fig.set_xaxis_label("Force ($F$)")
		fig.set_yaxis_label("Elastic Constant ($K = d(F) / d(R)$)")
	fig.set_logscale()
	fig.set_yaxis_ticks(minval = 0.9, maxval = 5000., nmajorticks = 3, nminorticks = 2)
	gen_plot(fig = fig, markersize = 2, linewidth = 1, show = True, save = False)
	exit()



## ARGUMENTS
# update the simulation results, even if the simulation results file already exists
update = ("update" in sys.argv)
# create all figures
figure = ("allfigs" in sys.argv)
# generate / re-generate figure 1
fig1 = ("fig1" in sys.argv)
# show figures on generation 
show = ("show" in sys.argv)
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
		FE_parms.to_csv(anal_dir + parm_file)
	else:
		FE_parms = pd.read_csv(anal_dir + parm_file)

	# process data 
	## IDEAL CHAIN

	## REAL CHAIN
	# estimate the linear elasticity constant
	i = 0
	for top in FE_parms['top'].unique(): # unique topology
		for k in FE_parms['K'].unique(): # strength of bending potential
			# for each, create a unique directory
			d = f"{top}_k{k:0.2f}"
			df = FE_parms[(FE_parms['top'] == top) & (FE_parms['K'] == k)]
			print(d)
			# get the no force change extension and then remove it from the dataframe
			r_nf = df.loc[0, 'E2Eproj_M1']
			df = df.reset_index().drop(0).dropna()
			# get the max chain extension
			r_max = df.loc[(len(df.index - 1)), 'E2Eproj_M1']

			## PLOT THE FORCE EXTENSION CURVE RvF and FvR
			df_norm = df.copy(deep = True)
			df_norm['F'] = df['F'] * r_nf
			df_norm['E2Eproj_M1'] = df_norm['E2Eproj_M1'] / r_nf
			fig = Figure()
			fig.load_data(d = df_norm, xcol = 'E2Eproj_M1', ycol = 'F')
			fig.set_xaxis_label(f"Normalized Chain Extension ($R_{{E2E}} / X$)")
			fig.set_yaxis_label(f"Normalized Force ($F \\cdot X$)")
			fig.set_title_label(f"Force Extension Curve ($X$ = {r_nf:.2f})")
			fig.set_subtitle_label(d)
			fig.set_yaxis_ticks(minval = 0.01, maxval = 500., nmajorticks = 3, nminorticks = 2)
			fig.set_logscale()
			fig.set_saveas(savedir = f"./02_processed_data/{job:s}/{d:s}/", filename = "FvR")
			gen_scatter(fig = fig, show = False, markersize = 20)

			fig = Figure()
			fig.load_data(d = df_norm, xcol = 'F', ycol = 'E2Eproj_M1')
			fig.set_yaxis_label(f"Normalized Chain Extension ($R_{{E2E}} / X$)")
			fig.set_xaxis_label(f"Normalized Force ($F \\cdot X$)")
			fig.set_title_label(f"Force Extension Curve ($X$ = {r_nf:.2f})")
			fig.set_subtitle_label(d)
			fig.set_xaxis_ticks(minval = 0.01, maxval = 500., nmajorticks = 3, nminorticks = 2)
			fig.set_logscale()
			fig.set_saveas(savedir = f"./02_processed_data/{job:s}/{d:s}/", filename = "RvF")
			gen_scatter(fig = fig, show = False, markersize = 20)


			## PARALLEL AND PERPENDICULAR ELASTICITY
			# plot the parallel and perpendicular elasticity against the applied force
			plot_elasticity_bondop(df = df, F_col = 'F', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_force")
			# same plot with normalization
			plot_elasticity_bondop(df = df, norm = float(r_nf), F_col = 'F', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_force_norm")
			## TODO get the transition force

			# plot the parallel and perpendicular elasticity against the chain extension
			plot_elasticity_bondop(df = df, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension")
			# same plot but with normalization
			plot_elasticity_bondop(df = df, norm = r_nf, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension_norm")
			# plot the nonlinear elasticity on a linear scale against the unitary chain extension
			plot_elasticity_bondop(df = df, norm = r_max, R_col = 'E2Eproj_M1', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp', show = show, save = True, path = f"./02_processed_data/{job:s}/{d:s}/", label = "elasticity_v_extension_max", logx = False)
			## TODO get the transition chain length
		
			## LINEAR ELASTICITY
			calc_elasticity_numerically(df = df_norm, F_col = 'F', R_col = 'E2Eproj_M1', monotonic = True)
			exit()
			# linear elastic constant is the y(x = 1) of line fit to data on log scale
			# use the transition chain extension to determine the linear elastic constant
			k_lin = calc_linear_elastic_constant(df = df, F_col = 'F', R_col = 'E2Eproj_M1')
			# i = i + 1
			# if i >= 4:
			#   # exit()

			## NUMERICAL ELASTIC CONSTANT
		

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
