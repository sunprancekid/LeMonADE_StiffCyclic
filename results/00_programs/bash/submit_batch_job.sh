#!/bin/bash
set -e

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
## PARAMETERS -- SLURM
# string represnting the maximum job time length
MAX_SLURM_TIME="05:00"

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
    if [[ $BOOL_FILE -eq 0 ]]; then
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
#     if [[ $VERBOSE -eq 1 ]]; then
#         ## TODO add verbose execution
#     fi
}

# method for generating a standardized slurm script
gen_slurm_script () {

    ## PARAMETERS
    # name of file contains slurm submission instructions
    local FILENAME="${SIMID}.slurm.sub"
    # directory that the file is stored in
    local FILEPATH="${JOBDIR}${SIMDIR}"

    ## ARGUMENTS
    # none

    ## SCRIPT
    # none

    echo "#!/bin/bash" > $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "#SBATCH -J ${SIMID}.%j.slurm" >> $FILEPATH$FILENAME
    echo "#SBATCH --nodes=1" >> $FILEPATH$FILENAME  # number of nodes
    echo "#SBATCH --ntasks=1" >> $FILEPATH$FILENAME   # number of processor cores (i.e. tasks)
#     echo "#SBATCH --time=${MAX_SLURM_TIME}" >> $FILEPATH$FILENAME # max simulation length
    echo "#SBATCH --error=${SIMID}.%j.err" >> $FILEPATH$FILENAME
    echo "#SBATCH --output=${SIMID}.%j.out" >> $FILEPATH$FILENAME
    # echo "#SBATCH --mail-type=FAIL" >> $FILEPATH$FILENAME
    # echo "#SBATCH --mail-user=dorsey@ipfdd.de" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### PARAMETERS ###" >> $FILEPATH$FILENAME
    echo "echo \"Current directory: \${SLURM_SUBMIT_DIR}\"" >> $FILEPATH$FILENAME
    echo "EXECDIR=/beetmp/dorsey/tmp/ # location of slurm execution" >> $FILEPATH$FILENAME
    echo "cd \${EXECDIR}" >> $FILEPATH$FILENAME
    echo "mkdir -p ${SIMID} " >> $FILEPATH$FILENAME
    echo "cd ${SIMID}" >> $FILEPATH$FILENAME
    echo "echo \"Running ${SIMID} on host \$(hostname) in \$(pwd)\"" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### load modules, executables" >> $FILEPATH$FILENAME
    echo "# module load slurm/15.08.8" >> $FILEPATH$FILENAME
    echo "# module load gcc/6.1.0" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### copy everything from the submit directory to the execute directory ###" >> $FILEPATH$FILENAME
    echo "cp $(pwd)/${JOBDIR}${SIMDIR}${SIMID}.sh ." >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### run job ####" >> $FILEPATH$FILENAME
    echo "echo \"Job start time is \$(date)." >> $FILEPATH$FILENAME
    echo "srun ./${SIMID}.sh -p /beetmp/dorsey/LeMonADE_StiffCyclic/build/bin/ > /dev/null 2>&1" >> $FILEPATH$FILENAME
    echo "echo \"Job end time is \$(date)." >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    ## TODO :: add instructions for copying results back
    echo "### copy back results, delete everthing ###" >> $FILEPATH$FILENAME
    echo "cp * $(pwd)/${JOBDIR}${SIMDIR}" >> $FILEPATH$FILENAME
    echo "cd ../" >> $FILEPATH$FILENAME
    echo "rm -rf ${SIMID}" >> $FILEPATH$FILENAME
}

## OPTIONS
while getopts "hvstj:p:f:l:e:c:" option; do
	case $option in
		h) # print script parameters to CLT
			help 0 ;;
        v) # exectue script verbosely
            declare -i BOOL_VERB=1 ;;
		s) # submit simulations to SLURM
			declare -i BOOL_SUBMIT_SLURM=1 ;;
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
    # login to linux server, update library and compile programs
    ssh $LINUXSERV "cd /beetmp/dorsey/LeMonADE_StiffCyclic/; git pull; ./results/00_programs/bash/build_LeMonADE.sh -l /beetmp/dorsey/LeMonADE/;"
    # make the job directory, which contains that submission scripts
    ssh $LINUXSERV "cd /beetmp/dorsey/; mkdir -p sub/${JOB};"
    # loop through parameter file, write slurm script to directory hirearchy
    declare -i N_LINES=$(wc -l < $PARMFILE)
    for i in $(seq 2 $N_LINES); do
        # get the simulation id (file column in file)
        SIMID=$(head -n ${i} ${PARMFILE} | tail -n 1 | cut -d , -f 1)
        # get the simulation directory (second column in file)
        SIMDIR=$(head -n ${i} ${PARMFILE} | tail -n 1 | cut -d , -f 2)

        # if a checkfile has been specified
        if [[ $BOOL_CHECKFILE -eq 1 ]]; then
            # if the check file exists in the simulation directory
            if [ -f ${JOBDIR}${SIMDIR}$CHECKFILE ]; then
                # skip execution of this set of parameters
                # continue to the next iteration of the loop
                if [[ $BOOL_VERB -eq 1 ]]; then
                    echo "checkfile ${JOBDIR}${SIMDIR}${CHECKFILE} already exists, moving to next set of simulation parameters."
                fi
                continue
            fi
        fi

        if [[ $BOOL_VERBOSE -eq 1 ]]; then
            echo "Submitting ${SIMID} in to ${LINUXSERV}."
        fi
        # write slurm script to local directory
        gen_slurm_script
        # copy slurm script to submit subdirectory on host cluster, submit
        ssh $LINUXSERV "cd /beetmp/dorsey/sub/${JOB}; cp $(pwd)/${JOBDIR}${SIMDIR}${SIMID}.slurm.sub .; sbatch ${SIMID}.slurm.sub;"
        if [[ $TEST_BOOL -eq 1 ]]; then
            # if testing the script, exit the program
            exit
        fi
    done
else
    # execute jobs locally
    # parse csv file, get simulation id and directory
    declare -i N_LINES=$(wc -l < $PARMFILE)
    for i in $(seq 2 $N_LINES); do
        # get the simulation id (file column in file)
        SIMID=$(head -n ${i} ${PARMFILE} | tail -n 1 | cut -d , -f 1)
        # get the simulation directory (second column in file)
        SIMDIR=$(head -n ${i} ${PARMFILE} | tail -n 1 | cut -d , -f 2)
        if [[ $BOOL_VERB -eq 1 ]]; then
            echo "Running ${SIMID} in ${JOBDIR}${SIMDIR}"
        fi

        # if a checkfile has been specified
        if [[ $BOOL_CHECKFILE -eq 1 ]]; then
            # if the check file exists in the simulation directory
            if [ -f ${JOBDIR}${SIMDIR}$CHECKFILE ]; then
                # skip execution of this set of parameters
                # continue to the next iteration of the loop
                if [[ $BOOL_VERB -eq 1 ]]; then
                    echo "checkfile ${JOBDIR}${SIMDIR}${CHECKFILE} already exists, moving to next set of simulation parameters."
                fi
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
