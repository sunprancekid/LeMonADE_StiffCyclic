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
# boolean that determines if a checkfile has been specified
declare -i BOOL_CHECKFILE=0
# boolean that determines test status of script execution
declare -i BOOL_TEST=0
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
	echo -e " -t           | submit one job as test in order to make sure everthing works properly."
# 	echo -e " -u           | upload local default directory to linux cluster."
# 	echo -e " -d           | sync local default directory with linux cluster."
	echo -e "\n"
	echo -e " ## SCRIPT PARAMETERS ##"
	echo -e " -j << ARG >> | MANDATORY: specify job title."
	echo -e " -p << ARG >> | specify path to simulation job directory hirearchy (default is 01_raw_data/\${JOB})."
	echo -e " -f << ARG >> | file which contains simulation parameters (default is 01_raw_data/\${JOB}/\${JOB}.csv)"
	echo -e " -l << ARG >> | specify remote linux cluster (default is ${LINUXSERV})."
	echo -e " -e << ARG >> | specify path to LeMonADE executables (default is ${EXECDIR})."
	echo -e " -c << ARG >> | check for file; if this file exists in the simulation directory, do not execute."
	echo -e "\n"
	## TODO Inform user about assumptions
	# first column if file contains simulation id
	# second column in file contains simulation directory
	# the simulation directory and job directorey (JOBDIR) should already exist
	# in the directory, the executable contains the instructions for executing, and is named after the simid

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
        JOBDIR="01_raw_data/${JOB}/"
    else
        # the path has been specified, assign whatever was passed through the option
        JOBDIR="${OPTPATH}"
    fi
    # check if the path exists. if it does not, make it
    if [ ! -d $JOBDIR ]; then
        echo "./submit_batch_job.sh::ERROR:: job directory ${JOBDIR} does not exist."
        help $NONZERO_EXITCODE
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

    # specify checkfile
    if [[ $BOOL_CHECKFILE -eq 1 ]]; then
        CHECKFILE="${OPTCHECKFILE}"
    fi

    # report instruvtions to user, if verbose execution
    if [[ $VERBOSE -eq 1 ]]; then
        ## TODO add verbose execution
    fi
}

## OPTIONS
while getopts "hvstj:p:f:l:e:c:" option; do
	case $option in
		h) # print script parameters to CLT
			help 0 ;;
        v) # exectue script verbosely
            declare -i BOOL_VERB=1 ;;
		s) # submit simulations to SLURM
			declare -i BOOL_SUBMIT_SLURM=1 ;;\
        t) # test submit status
            declare -i BOOL_TEST=1;;
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
        e) # update LeMonADE executables directory
            EXECDIR="${OPTARG}";;
        c) # specify checkfile
            declare -i BOOL_CHECKFILE=1
            OPTCHECKFILE="${OPTARG}";;
        \?) # sonstiges
            help $NONZERO_EXITCODE
    esac
done

## ARGUMENTS
# none

## SCRIPT
# check that correct arguments have been passed to the method
# initalize variables for script execution
check

# determine how to submit simulation jobs
if [[ $BOOL_SUBMIT_SLURM -eq 1 ]]; then
    # submit simulations to linux cluster
    echo "./submit_batch_job.sh::TODO:: implement method for submiting simulations on remote linux cluster."
    exit 0
else
    # execute jobs locally
    # parse csv file, get simulation id and directory
    declare -i N_LINES=$(wc -l < PARMFILE)
    for i in $(seq 2 $N_LINES); do
        # get the simulation id (file column in file)
        SIMID=$(head -n ${i} ${SIMPARAM_FILE} | tail -n 1 | cut -d , -f 1)
        # get the simulation directory (second column in file)
        SIMDIR=$(head -n ${i} ${SIMPARAM_FILE} | tail -n 1 | cut -d , -f 2)
        # if a checkfile has been specified
        if [[ $BOOL_CHECKFILE -eq 1 ]]; then
            # if the check file exists in the simulation directory
            if [ -f $CHECKFILE ]; then
                # skip execution of this set of parameters
                # continue to the next iteration of the loop
                continue
            fi
        fi

        # move to the simulation directory
		CURRDIR=$(echo $PWD)
        cd "${JOBDIR}${SIMDIR}"
        if [[ $BOOL_VERB -eq 1 ]]; then
            echo "Moving to ${JOBDIR}, executing ./${SIMID} .."
        fi
        # execute the simulation executable
        if [ ! -f ${SIMID}.sh ]; then
            echo "./submit_batch_job.sh::ERROR:: simulation executable ${SIMID}.sh does not exist."
        else
            ./${SIMID}.sh
        fi
        # return to the main directory
        cd $CURRDIR
    done
fi

# submit simulations
# either loop through locally and submit serial
# or logon to computing cluster and exectue
