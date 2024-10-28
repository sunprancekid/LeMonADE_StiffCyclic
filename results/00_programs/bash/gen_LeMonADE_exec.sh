!/bin/bash
set -e

## Matthew Dorsey
## 28.10.2024
## generates instructions for executing LeMonADE simulation

## PARAMETERS
# none

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

# OPTIONS
while getopts "h" option; do
	case $option in
		h) # print script parameters to CLT
			help
			exit 0 ;;
        \?) # sonstiges
            help
            exit $NONZERO_EXITCODE
    esac
done

## ARGUMENTS
# none

## SCRIPT
# none
