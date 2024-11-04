#!/bin/bash
set -e

## Matthew Dorsey
## 28.10.2024
## generates instructions for executing LeMonADE simulation

## PARAMETERS
# boolean determining if the simulation simulation executable is for real or ideal chains
declare -i BOOL_IDEAL=0 # assume false

## METHODS
# output instructions for script execution
help () {

    ## PARAMETERS
    # name of script
    local script_name="gen_LeMonADE_exec"

    ## OPTIONS
    # none

    ## ARGUMENTS
    # none

    ## SCRIPT
	echo -e "\nScript for generating LeMonADE execution instructions.\nUSAGE: ./${script_name}.sh << FLAGS >>"
	echo -e "\n"
	echo -e " ## SCRIPT PROTOCOL ##"
	echo -e " TODO :: add flags corresponding to simulation parameters."
	echo -e "\n"
}

# build script constants based on parameters passed to script
build () {

    ## PARAMETERS
    # none

    ## OPTIONS
    # none

    ## ARGUMENTS
    # none

    ## SCRIPT
    # check ideal vs. simulation
    if [[ BOOL_IDEAL -eq 1 ]]; then
        # title of ideal chain simulation script
        SIMULATION_EXECUTABLE="simulateIdealPolymerBFM"
    else
        # title of real chain simulation script
        SIMULATION_EXECUTABLE="simulateRealPolymerBFM"
    fi
}

# generates execution script
gen_exec () {

    ## PARAMETERS
    # name of file containing simulation execution instructions
    local exec_file="${SIMPATH}${SIMID}.sh"

    ## OPTIONS
    # none

    ## ARGUMENTS
    # none

    ## SCRIPT
    # write files directory, generate simulation executables
    echo "#!/bin/bash" > ${exec_file}
    echo "set -e" >> ${exec_file}
    echo "PATH=\"./\"" >> ${exec_file}
    echo "while getopts \"p:\" option; do" >> ${exec_file}
    echo -e "\tcase \$option in" >> ${exec_file}
    echo -e "\t\tp)" >> ${exec_file}
    echo -e "\t\t\tPATH=\${OPTARG}" >> ${exec_file}
    echo -e "\tesac" >> ${exec_file}
    echo "done" >> ${exec_file}
    echo "\${PATH}generatePolymerBFM ${GENFLAGS}" >> ${exec_file}
    echo "\${PATH}${SIMULATION_EXECUTABLE} ${SIMFLAGS}" >> ${exec_file}
    chmod u+x ${exec_file}
}

# OPTIONS
# parse options
while getopts "hi" option; do
	case $option in
		h) # print script parameters to CLT
			help
			exit 0 ;;
        i) # simulation for ideal chain rather than real chain
            declare -i BOOL_IDEAL=1
            ;;
        \?) # sonstiges
            help
            exit $NONZERO_EXITCODE
    esac
done


## ARGUMENTS
# first argument: simulation path
SIMPATH=$1
# second argument: simulation id
SIMID=$2
# third argument: flags for simulation generation
GENFLAGS="$3"
# fourth argument: flags for simulation simulation
SIMFLAGS="$4"

## SCRIPT
# build constants based on options passed to script
build

# generate the executable instructions
gen_exec
