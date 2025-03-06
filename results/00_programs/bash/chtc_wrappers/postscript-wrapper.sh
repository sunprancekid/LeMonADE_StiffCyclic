#!/bin/bash 
set -e 

# Matthew Dorsey
# 2023-04-20
# wrapper for post-scripts of DAGMAN annealing jobs. 
# used to write status of DAGMAN jobs, standard out and error of chtc scripts.
# script moves files to their proper storage directory if certain criteria are met.

## PARAMETERS
# exit code for indicating to DAGMAN workflow that there is an error
declare -i NONZERO_EXITCODE=1
# boolean for checking the simulation length in MCS
declare -i BOOL_MCS=0
# passed as argument to check mcs time, minimum number of mcs simulation needs to surpass
declare -i MIN_MCS=0
# boolean for determining the variance of a columns values
declare -i BOOL_VAR=0
# passed as argument to check variance, maximum variance allowed within a column
declare -i MAX_VAR=0

## OPTIONS
# determine exit criteria
# flags are used for specifying exit criteria 
while getopts "t:v:" option; do
    case $option in
    	t) # flag for checking the length of the simulation
			
			# boolean determining if the script checks the simulation length
			declare -i BOOL_MCS=1 
			# value passed to 
			declare -i MIN_MCS=${OPTARG};;
		v) # flag for checking the variance of a property value

			# boolean determining if the script checks the variance of a property value
			declare -i BOOL_VAR=1
			# maximum variance that property value muss be less than
			declare -i MAX_VAR=${OPTARG};;
		\?) # default
			echo "must declare arguments" ;;
   esac
done
shift $((OPTIND-1))


## ARGUMENTS
# first argument: id associated with the annealing simulation
SIMID=$1
# second argument: exit code of the script
EXIT_VAL=$2
# third argument: number of node retries
RETRY_VAL=$3


## FUCNTIONS
# none


# SCRIPT
# establish the name used to write standard out and error from
# the pre and post processing scripts
OUTNAME="postscript_stdout.txt"

# determine the operation to perform
if [[ $BOOL_MCS -eq 1 ]]; then 
	./check_value.sh -f RE2E.dat -c 1 -m ${MIN_MCS} $SIMID $EXIT_VAL $RETRY_VAL >> ${OUTNAME} 2>&1
fi

exit 0
