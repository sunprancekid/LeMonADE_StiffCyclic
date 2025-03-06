# !/bin/bash
set -e

# Matthew Dorsey
# @sunprancekid
# checks that certain property values meet certain criteria.
# if they do the script exists successfully, otherwise the script exits with an error code.

## PARAMETERS
# boolean determining if a file name has been passed to the script
declare -i BOOL_FILE=0
# file name passed to method
CHECKFILE=""
# boolean determining of column corresponding to the property that could be checked
declare -i BOOL_COL=0
# column number corresponding to property value that should be checked
declare -i CHECK_COL=0
# boolean determining if a minimum value should be checked
declare -i BOOL_CHECK_MIN=0
# minimum value that propertu must surpace
declare -i MIN_VAL=0


## OPTIONS
# determine check criteria
while getopts "f:c:m:" option; do
    case $option in
    	f) # flag parsing file name
			
			# boolean determining that the file name has been passed
			declare -i BOOL_FILE=1 
			# file name
			CHECKFILE=${OPTARG};;

		c) # flag for parsing the check column

			# boolean determining that the column has been passed
			declare -i BOOL_COL=1
			# maximum variance that property value muss be less than
			declare -i CHECK_COL=${OPTARG};;

		m) # flag for parsing minimum value

			# boolean determining that the minimum value has been parse
			declare -i BOOL_CHECK_MIN=1
			# minimum value that the property must surpass
			declare -i MIN_VAL=${OPTARG};;
		\?) # default
			echo "must declare arguments." ;;
   esac
done
shift $((OPTIND-1))


## ARGUMENTS
# first argument: id associated with the annealing simulation
SIMID=$1
# second argument: exit code of the node being processed
RETURN_VAL=$2
# third argument: number of times the node has iterated
RETRY_VAL=$3 


## FUNCTIONS
# TODO add help, listing options and arguments
# TODO add check, that all of the correct parameters have been provided

# function that parses the number of lines in a file
# provided the name of the file
get_lines() {

	## PARAMETERS
	# none

	## OPTIONS
	# none

	## ARGUMENTS
	# name of file to parse lines from
	local filename=$1

	## SCRIPT
	# parse the number of lines from the file
	declare -i N_LINES=$(wc -l < $filename)

	# return the numbe of lines in the csv to the user
	echo $N_LINES
}


# function that gets a specific line from a csv file
get_line() {

	## PARAMETERS
	# none

	## OPTIONS
	# none

	## ARGUMENTS
	# name of file to parse lines from
	local filename=$1
	# line number of parse from file
	declare -i line_no=$2

	## SCRIPT 
	# parse the line from the file
	line=$( head -n ${line_no} ${filename} | tail -n 1 )

	# return the line
	echo $line
}


## SCRIPT
# TODO add check
# check

# parse the final line from the column
declare -i N=$(get_lines ${CHECKFILE})
get_lines ${CHECKFILE} ${N}