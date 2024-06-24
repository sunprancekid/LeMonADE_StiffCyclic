## Matthew A. Dorsey
## Intitute for Polymer Resarch - Institute for Polymer Theory
## North Carolina State University
## 2024.06.24

## script that submits a set of parameterized jobs for execution.
## the script on accepts a file that submits the simulations for execution.
## the script can either run on a local cluster, or submit on HPC cluster at IPF / ITP.

## PARAMETERS
# boolean that determines if the script should execute verbosely
declare -i BOOL_VERB=0
# boolean that determines if a parameterized file has been passed to the method
declare -i BOOL_FILE=0
# boolean that determines if a path for the simulation directory has been passed to the method
declare -i BOOL_PATH=0
# boolean that determines if a job title has been passed to the method
declare -i BOOL_JOB=0
# boolean that determines if the job should be submitted on a remote linux cluster
declare -i BOOL_SUBMIT_SLURM=0
# default linux server used for remote jobs
LINUXSERV="gandalf"
# default path to parameterized linux executables for job
EXECDIR="00_programs/build/bin/"

## FUNCTIONS
# report script used to command line terminal
help() {

    ## PARAMETERS
    # none

    ## ARGUMENTS
    # first argument: exit code
    local exitcode=$1

    ## SCRIPT
    # report usage
	echo -e "\nScript for executing parameterized jobs on linux clusters.\nUSAGE: ./submit_batch_job.sh << FLAGS >>\n"
	echo -e "\n"
	echo -e " ## SCRIPT PROTOCOL ##"
	echo -e " -h           | display script options, exit 0."
	echo -e " -v           | execute script verbosely."
	echo -e " -s           | submit job to SLURM via remote linux cluster (otherwise run locally in serial)."
# 	echo -e " -u           | upload local default directory to linux cluster."
# 	echo -e " -d           | sync local default directory with linux cluster."
	echo -e "\n"
	echo -e " ## SCRIPT PARAMETERS ##"
	echo -e " -j << ARG >> | MANDATORY: specify job title."
	echo -e " -p << ARG >> | specify path to simulation job directory hirearchy (default is 01_raw_data/\${JOB})."
	echo -e " -f << ARG >> | file which contains simulation parameters (default is 01_raw_data/\${JOB}/\${JOB}.csv)"
	echo -e " -l << ARG >> | specify remote linux cluster (default is ${LINUXSERV})."
	echo -e " -e << ARG >> | specify path to LeMonADE executables (default is ${EXECDIR})."
	echo -e "\n"

    # exit with exit code
    exit $exitcode
}

# checks that correct options have been passed to method before script execution
# initlaized variables for unspecified options
check () {

    ## PARAMETERS
    # none

    ## ARGUMENTS
    # none

    ## SCRIPT
    # check option booolean
    if [[ $BOOL_JOB -eq 0 ]]; then
        # if the job name was not specified, report to user and exit the program
        echo "./submit_batch_job.sh::ERROR:: must specify job name. See option arguments.\n"
        help $NONZERO_EXITCODE
    fi

    # initialize the path to the simulation directory hirearchy
    if [[ $BOOL_PATH -eq 0 ]]; then
        # the path has not been specified, assign the default
        PATH="01_raw_data/${JOB}/"
    else
        # the path has been specified, assign whatever was passed through the option
        PATH="${OPTPATH}"
    fi
    # check if the path exists. if it does not, make it
    if [ ! -d $PATH ]; then
        mkdir -p $PATH
    fi

    # initialize the file that contains the simulation parameters
    if [[ BOOL_FILE -eq 0 ]]; then
        # the file was not specified by the user, assign the default
        PARMFILE="01_raw_data/${JOB}/${JOB}.csv"
    else
        # the file was specified by the user, assign that value
        PARMFILE="${OPTFILE}"
    fi
    # check that the file exists
    if [ ! -f $PARMFILE ]; then
        # if the file does not exist, report to user and exit program
        echo "./submit_batch_job.sh::ERROR:: simulation parameter file ${FILE} does not exist.\n"
        help $NONZERO_EXITCODE
    fi

    # report instruvtions to user, if verbose execution
    if [[ $VERBOSE -eq 1 ]]; then
        ## TODO add verbose execution
    fi
}

## OPTIONS
while getopts "hvsj:p:f:l:e:" option; do
	case $option in
		h) # print script parameters to CLT
			help 0 ;;
        v) # exectue script verbosely
            declare -i BOOL_VERB=1 ;;
		s) # submit simulations to SLURM
			declare -i BOOL_SUBMIT_SLURM=1 ;;
		j) # update job name
            declare -i BOOL_JOB=1
			JOB="${OPTARG}" ;;
        p) # update path
            declare -i BOOL_PATH=1
            OPTPATH="${OPTARG}";;
        f) # update parameter file
            declare -i BOOL_FILE=1
            OPTFILE="${OPTARG}";;
		l) # update linux cluster
			LINUXSERV="${OPTARG}" ;;
        e) # update LeMonADE executables
            EXECDIR="${OPTARG}";;
        \?) # sonstiges
            help $NONZERO_EXITCODE
    esac
done

## ARGUMENTS
# none

## SCRIPT
# check that correct arguments have been passed to the method
check
# parse csv file, get simulation id and directory
# submit simulations
# either loop through locally and submit serial
# or logon to computing cluster and exectue
