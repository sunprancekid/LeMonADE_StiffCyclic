#!/bin/bash
set -e

## Matthew Dorsey
## 2024.06.11
##  script for generating constant force extension simulations for real and ideal polymer chains



## PARAMETERS
## PARAMETERS -- BOOLEAN
# boolean determining if the script should be executed verbosely
declare -i BOOL_VERB=0
# boolean determining if the script should generate the simulation parameters
declare -i BOOL_GEN=0
# boolean determining if the script should submit simulations to SLURM
declare -i BOOL_SUB=0
# boolean determining if the script should upload a directory to the linux server
declare -i BOOL_UPP=0
# boolean determining if the script should download a directory from the could
declare -i BOOL_DWN=0
## PARAMETERS -- JOB
# number of unique forces to test
declare -i N_FORCE_VAL=30
# maxmium force to test
MAX_FORCE_VAL="5"
# minimum force to test
MIN_FORCE_VAL=".0001"
# array containing N to test
PARM_N=( 100 )
# array containing the different executables to test
PARM_EX=( "FElinearChainReal" "FElinearChainIdeal" )
# number of times to replicate each unique set of simualation conditions
declare -i PARM_R=1
# default directory for upload / download, generating parameters
MAINDIR="01_raw_data"
# default job name
JOB="forceExtension"
# linux server
LINUXSERV="gandalf"
# path to LeMonADE executables
EXECDIR="00_programs/build/bin/"
## PARAMETERS -- SIMULATION
# number of MCSs between each property calculation
declare -i N_MCS=100000
# number of time properties are calculation
declare -i N_runs=10000
# number of MCSs before equilibrium properties are calculated
declare -i t_equilibrium=100000000


## FUNCTIONS
# display script options
help () {
	echo -e "\nScript for generating polymer simulations with BFM model on linux clusters.\nUSAGE: ./gen_scaling simulations.sh << FLAGS >>\n"
	echo -e "\n"
	echo -e " ## SCRIPT PROTOCOL ##"
	echo -e " -h           | display script options, exit 0."
	echo -e " -v           | execute script verbosely."
	echo -e " -g           | generate simulation parameters / directories ."
	echo -e " -s           | submit job to SLURM."
	echo -e " -u           | upload local default directory to linux cluster."
	echo -e " -d           | sync local default directory with linux cluster."
	echo -e "\n"
	echo -e " ## SCRIPT PARAMETERS ##"
	echo -e " -j << ARG >> | specify job title (default is ${JOB})."
	echo -e " -l << ARG >> | specify linux cluster (default is ${LINUXSERV})."
	echo -e " -p << ARG >> | specify the local directory (default is ${MAINDIR})."
	echo -e " -e << ARG >> | specify path to LeMonADE executables (default is ${EXECDIR})."
	echo -e "\n"
}

# generate simulation parameters
generate_simulation_parameters() {
    ## TODO implement simulations parameters
}

## OPTIONS
while getopts "hvgsudj:l:p:" option; do
	case $option in
		h) # print script parameters to CLT
			help
			exit 0 ;;
        v) # exectue script verbosely
            declare -i BOOL_VERB=1 ;;
        g) # generate simulation parameters
            declare -i BOOL_GEN=1 ;;
		s) # submit simulations to SLURM
			declare -i BOOL_SUB=1 ;;
		u) # upload local directory to linux cluster
			declare -i BOOL_UPP=1 ;;
		d) # sync local directory with linux cluster
			declare -i BOOL_DWN=1 ;;
		j) # update job name
			JOB="${OPTARG}" ;;
		l) # update linux cluster
			LINUXSERV="${OPTARG}" ;;
		p) # update local directory path
			MAINDIR="${OPTARG}" ;;
        \?) # sonstiges
            help
            exit $NONZERO_EXITCODE
    esac
done


## ARGUMENTS
# none


## SCRIPT
# generate parameters save to file
if [ $BOOL_GEN -eq 1 ]; then
	gen_simparm
fi
