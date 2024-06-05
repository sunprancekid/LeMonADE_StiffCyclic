#!/bin/bash
set -e

## Matthew Dorsey
## 05.06.2024
## Programs for generating simulation parameters, directory hirarchy, and execution instructions


## PARAMETERS
## PARAMETERS -- BOOLEAN
# boolean determining if the script should be executed verbosely
declare -i BOOL_VERBOSE=0
# boolean determining if the script should generate the simulation parameters
declare -i BOOL_GEN=0
# boolean determining if the script should upload a directory to the linux server
declare -i BOOL_UP=0
# boolean determining if the script should download a directory from the could
declare -i BOOL_DWN=0
## PARAMETERS -- JOB
# array containing N to test
declare -A PARM_N=( 32 48 64 78 96 128 256 )
# array containing the different executables to test
declare -A PARM_EX=( "linearChainReal" "linearChainIdeal" )
# number of times to replicate each unique set of simualation conditions
declare -i PARM_R=3
# default directory for upload / download, generating parameters
MAINDIR="01_raw_data"
# default job name
JOB="scaling"
# executable
# ?
# linux server
LINUXSERV="gandalf"


## FUNCTIONS
# display script options
help () {

	# TODO :: write options so that if a variable is specified, its value is held constant
	# if the job already exists, the script skips, unless and overwrite flag is specified

	echo -e "\nScript for generating conH jobs on CHTC systems.\nUSAGE: ./conH.sh << FLAGS >>\n"
	echo -e " -h           | display script options, exit 0."
	echo -e " -v           | execute script verbosely."
	echo -e " -g           | boolean determing if simulation parameters / directories should be generated."
	echo -e " -s           | submit job to CHTC based on current status."
	echo -e " -o           | boolean determining if existing jobs should be overwritten."
	echo -e "\n"
	echo -e " ## JOB PARAMETERS ##"
	echo -e " -j << ARG >> | specify job title (default is ${JOB})"
	echo -e " -c << ARG >> | specify cell size (default is ${CELL})"
	echo -e " -e << ARG >> | specify events per (default is ${EVENTS})"
	echo -e " -f << ARG >> | specify annealing fraction (default is ${FRAC})"
	echo -e "\n"
	echo -e " ## SIMULATION PARAMETERS ##"
	echo -e " -r << ARG >> | integer representing the number of replicates to perform (default is 1)."
	echo -e "\n"
	echo -e " ## CHTC SUBMIT INSTRUCTIONS ##"
	echo -e " -t           | \"touch\" simulation directries, update files."
	echo -e " -r           | rerun anneal simulations that have already been performed."
}


## OPTIONS
while getopts "vg:" option; do
	case $option in
        v) # exectue script verbosely
            declare -i BOOL_VERBOSE=1
            ;;
        g) # generate simulation parameters
            declare -i BOOL_GEN=1
            ;;
        \?) # sonstiges
            help
            exit $NONZERO_EXITCODE

    esac
done


## ARGUMENTS
# none


## SCRIPT
# generate parameters save to file

