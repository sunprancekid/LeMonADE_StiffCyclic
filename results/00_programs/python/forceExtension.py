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


## PARAMETERS
# none

## METHODS
# calculate parallel and perpendicular non-linear elasticity constants from bond order parameters
def plot_elasticity_bondop (df = None, F_col = None, K_parallel_col = None, K_perp_col = None, iso_col = None):
	
	# check that the correct parameters were passed to the method
	if df is None:
		# must pass dataframe to the method
		print("analysis :: forceExtension :: plot_elasticity_bondop :: ERROR :: must pass dataframe to method as 'df'.")

	if F_col is None:
		# must specifcy column containing force
		print("analysis :: forceExtension :: plot_elasticity_bondop :: ERROR :: must specifiy column in 'df' containing force values as 'F_col'.")
	
	if K_parallel_col is None and K_perp_col is None:
		# must specify at least the parallel or perpendicular columns containing the bond order parameter
		print("analysis :: forceExtension :: plot_elasticity_bondop :: ERROR :: must specify either column containing parallel bond order parameter (as 'K_parallel_col') or column containing perpendicular bond order parameter (as 'K_parallel_col'), or both.")

	# loop through the each line in the df, process the data
	f = [] 			# force associated with chain extension
	k_parallel = [] # linear elastic constant in parallel direction
	k_perp = [] 	# linear elastic constant in perpendicular direction
	i = [] 			# from isloation column
	for i, r in df.iterrows():
		f.append(r[F_col])
		# determine the parallel linear elastic constant
		if K_parallel_col is not None:
			k_parallel.append(1. / r[K_parallel_col + "_var"])
		# determine the perpendicular elastic constant
		if K_perp_col is not None:
			k_perp.append(1. / r[K_perp_col + "_var"])
		# pull  the isolation tag
		if iso_col is not None:
			i.append(r[iso_col])

	# plot the data

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
		FE_parms = parse_results(parms = FE_parms, dir = data_dir, simfile = 'RE2E.dat', col = 4, title = 'E2Etot', M1 = True, M2 = True, var = True, bootstrapping = True, tabsep = True, plot = False)
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
				plot_force_extension (k_df, Y_col = 'E2Eproj', X_col = 'F', iso_col = 'top', isolabel = '{:s}', X_label = "Normalized External Force ($X \\cdot f$)", Y_label = "Normalized Chain Extension ($X^{{-1}} \\cdot R_{{E2E}}$)", saveas = save_name + '_norm_data', plot_data = True, y_max = 10., y_min = 0.001, x_min = 0.01, x_max = 100., show = show, logscale_x = True, logscale_y = True, plot_slope = False, norm = True, pincus = True, hookean = True) # saveas = save_name + '_norm_data.png',
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

	# estimate the linear elasticity constant
	for top in FE_parms['top'].unique(): # unique topology
		for k in FE_parms['K'].unique(): # strength of bending potential
			# for each, create a unique directory
			d = f"{t}_k{k:0.2f}"
			## PARALLEL AND PERPENDICULAR ELASTICITY
			# plot the parallel and perpendicular elasticity
			plot_elasticity_bondop(df = FE_parms[(FE_parms['top'] == top) & (FE_parms['K'] == k)], F_col = 'F', K_parallel_col = 'cos_theta_parallel', K_perp_col = 'cos_theta_perp')
			# plot the force extension curve along with the numerically calcualted elasticity
			exit()

	# parallel and perpendicular elasticity
	# elasticity from force-extension curve
