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
# boolean that determintes if the job is being submitted on a chtc cluster
declare -i BOOL_SUBMIT_CHTC=0
# boolean that determines if a checkfile has been specified
declare -i BOOL_CHECKFILE=0
# boolean that determines test status of script execution
declare -i BOOL_TEST=0
# boolean that determines if simulation state should be overwritten
declare -i BOOL_OVERWRITE=0
# default linux server used for remote jobs
LINUXSERV="gandalf"
# default path to parameterized linux executables for job
EXECDIR="00_programs/build/bin/"
# default singularity image used with chtc jobs
DEFAULT_SINGULARITY="l.sif"
# maximum MCS steps for equilibriation period on CHTC systems
declare -i MAX_MCS_EQUIL=1000000000 # 1 billion
# maximum MCS steps for running loop on CHTC systems
declare -i MAX_MCS_RUN=500000000 # 500 million
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
    echo -e " -o           | restart simulations from beginning; otherwise actual state of simulations are saved, if they exists."
	echo -e " -s           | submit job to SLURM via remote linux cluster."
    echo -e " -z           | submit job on CHTC cluster (most be logged in)."
	echo -e " -t           | submit one job as test in order to make sure everthing works properly."
    # echo -e " -u           | upload local default directory to linux cluster."
    # echo -e " -d           | sync local default directory with linux cluster."
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
    echo "echo \"Job start time is \$(date).\"" >> $FILEPATH$FILENAME
    echo "srun ./${SIMID}.sh -g -p /beetmp/dorsey/LeMonADE_StiffCyclic/build/bin/ > /dev/null 2>&1" >> $FILEPATH$FILENAME
    echo "srun ./${SIMID}.sh -e -p /beetmp/dorsey/LeMonADE_StiffCyclic/build/bin/ > /dev/null 2>&1" >> $FILEPATH$FILENAME
    echo "srun ./${SIMID}.sh -r -p /beetmp/dorsey/LeMonADE_StiffCyclic/build/bin/ > /dev/null 2>&1" >> $FILEPATH$FILENAME
    echo "echo \"Job end time is \$(date).\"" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    # add instructions for copying results back
    echo "### copy back results, delete everthing ###" >> $FILEPATH$FILENAME
    echo "cp * $(pwd)/${JOBDIR}${SIMDIR}" >> $FILEPATH$FILENAME
    echo "cd ../" >> $FILEPATH$FILENAME
    echo "rm -rf ${SIMID}" >> $FILEPATH$FILENAME
    echo "cp ../sub/${JOB}/${SIMID}.*.err $(pwd)/${JOBDIR}${SIMDIR} " >> $FILEPATH$FILENAME
    echo "cp ../sub/${JOB}/${SIMID}.*.out $(pwd)/${JOBDIR}${SIMDIR} " >> $FILEPATH$FILENAME
}

# method for generating standardized chtc scripts for bfm jobs
gen_chtc_scripts () {

    ## PARAMETERS
    # path from current working directory to simulation directory
    D="${JOBDIR}${SIMDIR}"
    # list of sub directories to generate inside the main directory
    SUBDIR=( "node" "node/init" "out" "sub" "sub/exec" "anal" )
    # name of the executable file
    EXEC_NAME="${SIMID}.sh"

    ## ARGUMENTS
    # none

    ## SCRIPT

    # establish the subdag path
    SUBDAG="${SIMID}.spl"
    SUBDAG_PATH="${D}${SUBDAG}"
    if [[ -f "$SUBDAG_PATH" ]]; then 
        # if the subdag already exists, remove it
        rm "$SUBDAG_PATH"
    fi

    # generate subdirectories
    for sd in "${SUBDIR[@]}"; do 
        if [[ -d ${D}${sd} ]]; then
            # if the directory already exists, remove it
            rm -r ${D}${sd}
        fi
        mkdir -p "${D}${sd}/"
    done

    # establish the initialization node
    gen_chtc_init

    # establish the equilibriation node
    gen_chtc_equil

    # establish the running node
    # gen_chtc_run

}

# method for generating initialization node
gen_chtc_init () {

    ## PARAMETERS
    # name of file containing submission instructions
    local SUB_NAME="sub/init.sub"
    # path to file containing submission instructions
    local SUB_PATH="${D}${SUB_NAME}"
    # number of retries attempted by init node
    declare -i NUM_RETRY=5
    # name of the executable file
    # EXEC_NAME="sub/exec/${SIMDIR}.sh"

    ## PARAMETERS - FILES
    # configuration file
    local SIM_CONFIG="config_gen.dat"

    ## PARAMETERS - SUBMISSION INTRUCTIONS
    # memory to request
    local REQUEST_MEMORY="500MB"
    # disk space to request
    local REQUEST_DISK="1GB"
    # directory that output files are remapped to
    local REMAP="node/init/"
    # list of files with remapping instructions
    local RMP_SIM_CONFIG="${SIM_CONFIG}=${REMAP}${SIM_CONFIG}"
    # list of files that should be transfered to the execute node
    local TRANSFER_INPUT_FILES="${EXEC_NAME}"
    # list of files that should be transfered from the execute node
    local TRANSFER_OUTPUT_FILES="${SIM_CONFIG}"
    # list of remap instructions for each output file
    local TRANSFER_OUTPUT_REMAPS="${RMP_SIM_CONFIG}"

    ## ARGUMENTS
    # none

    # SCRIPT
    # write submission script
    echo "executable = ${EXEC_NAME}" > $SUB_PATH
    echo "arguments = -g -p /LeMonADE_StiffCyclic/build/bin/" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "+SingularityImage = \"/home/mad/tmp/${DEFAULT_SINGULARITY}\"" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "should_transfer_files = YES" >> $SUB_PATH
    echo "transfer_input_files = ${TRANSFER_INPUT_FILES}" >> $SUB_PATH
    echo "transfer_output_files = ${TRANSFER_OUTPUT_FILES}" >> $SUB_PATH
    # echo "transfer_output_remaps = \"${TRANSFER_OUTPUT_REMAPS}\"" >> $SUB_PATH
    echo "when_to_transfer_output = ON_SUCCESS" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "log = out/init.log" >> $SUB_PATH
    echo "error = out/init.err" >> $SUB_PATH
    echo "output = out/init.out" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "request_cpus = 1" >> $SUB_PATH
    echo "request_disk = ${REQUEST_DISK}" >> $SUB_PATH
    echo "request_memory = ${REQUEST_MEMORY}" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "on_exit_hold = (ExitCode != 0)" >> $SUB_PATH
    echo "requirements = (HAS_GCC == true) && (Mips > 30000)" >> $SUB_PATH
    # echo "requirements = HasSingularity" >> $SUB_PATH
    echo "+ProjectName=\"NCSU_Hall\"" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "queue" >> $SUB_PATH

    # if overwrite is off and the generation config file already exists, can skip adding this node to the subdag
    if [[ $BOOL_OVERWRITE -eq 0 ]]; then
        if [[ -f ${D}${SIM_CONFIG} ]]; then
            # skip to next node without adding to subdag
            return
        fi
    fi

    # add node to subdag
    echo "JOB ${SIMID}_init ${SUB_NAME}" >> $SUBDAG_PATH
    echo "RETRY ${SIMID}_init ${NUM_RETRY}" >> $SUBDAG_PATH
}

# method for generation equilibriation node
gen_chtc_equil () {

    ## PARAMETERS
    # name of file containing submission instructions
    local SUB_NAME="sub/equil.sub"
    # path to file containing submission instructions
    local SUB_PATH="${D}${SUB_NAME}"
    # number of retries attempted by init node
    declare -i NUM_RETRY=5
    # name of the executable file
    # EXEC_NAME="sub/exec/${SIMDIR}.sh"

    ## PARAMETERS - FILES
    # configuration files
    local SIM_CONFIG_EQUIL="config_equil.dat" # equilibrium configuration
    local SIM_CONFIG_GEN="config_gen.dat" # initial configuration, from previous node
    # contains end-to-end instructions
    local SIM_E2E="RE2E.dat"

    ## PARAMETERS - SUBMISSION INTRUCTIONS
    # memory to request
    local REQUEST_MEMORY="500MB"
    # disk space to request
    local REQUEST_DISK="1GB"
    # directory that output files are remapped to
    local REMAP="node/equil/"
    # list of files with remapping instructions
    local RMP_SIM_CONFIG_EQUIL="${SIM_CONFIG_EQUIL}=${REMAP}${SIM_CONFIG_EQUIL}"
    local RMP_SIM_RE2E="${SIM_RE2E}=${REMAP}${SIM_RE2E}"
    # list of files that should be transfered to the execute node
    local TRANSFER_INPUT_FILES="${EXEC_NAME}, ${SIM_CONFIG_GEN}"
    # list of files that should be transfered from the execute node
    local TRANSFER_OUTPUT_FILES="${SIM_CONFIG_EQUIL}, ${SIM_RE2E}"
    # list of remap instructions for each output file
    local TRANSFER_OUTPUT_REMAPS="${RMP_SIM_CONFIG_EQUIL}"

    ## ARGUMENTS
    # none

    ## SCRIPT
    # write submission script
    echo "executable = ${EXEC_NAME}" > $SUB_PATH
    echo "arguments = -e -p /LeMonADE_StiffCyclic/build/bin/ -n ${MAX_MCS_EQUIL}" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "+SingularityImage = \"/home/mad/tmp/${DEFAULT_SINGULARITY}\"" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "should_transfer_files = YES" >> $SUB_PATH
    echo "transfer_input_files = ${TRANSFER_INPUT_FILES}" >> $SUB_PATH
    echo "transfer_output_files = ${TRANSFER_OUTPUT_FILES}" >> $SUB_PATH
    # echo "transfer_output_remaps = \"${TRANSFER_OUTPUT_REMAPS}\"" >> $SUB_PATH
    echo "when_to_transfer_output = ON_SUCCESS" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "log = out/equil.log" >> $SUB_PATH
    echo "error = out/equil.err" >> $SUB_PATH
    echo "output = out/equil.out" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "request_cpus = 1" >> $SUB_PATH
    echo "request_disk = ${REQUEST_DISK}" >> $SUB_PATH
    echo "request_memory = ${REQUEST_MEMORY}" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "on_exit_hold = (ExitCode != 0)" >> $SUB_PATH
    echo "requirements = (HAS_GCC == true) && (Mips > 30000)" >> $SUB_PATH
    # echo "requirements = HasSingularity" >> $SUB_PATH
    echo "+ProjectName=\"NCSU_Hall\"" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "queue" >> $SUB_PATH

    # if overwrite is off and the generation config file already exists
    # can skip adding this node to the subdag
    if [[ $BOOL_OVERWRITE -eq 0 ]]; then
        if [[ -f ${D}${SIM_CONFIG_EQUIL} ]]; then
            # skip to next node without adding to subdag
            return
        fi
    fi

    # add node to subdag
    echo "JOB ${SIMID}_equil ${SUB_NAME}" >> $SUBDAG_PATH
    echo "RETRY ${SIMID}_equil ${NUM_RETRY}" >> $SUBDAG_PATH

    # if overwrite is off and the previous config file exists
    if [[ $BOOL_OVERWRITE -eq 0 && -f ${D}${SIM_CONFIG_GEN} ]]; then
        # init node is not running, do not need to establish parent chile relationship
        return
    else
        # establish parent child relationship
        echo "PARENT ${SIMID}_init CHILD ${SIMID}_equil" >> ${SUBDAG_PATH}
    fi
}

# node for running bfm model, generating statistics
gen_chtc_run () {

    ## PARAMETERS
    # name of file containing submission instructions
    local SUB_NAME="sub/run.sub"
    # path to file containing submission instructions
    local SUB_PATH="${D}${SUB_NAME}"
    # number of retries attempted by init node
    declare -i NUM_RETRY=500
    # name of the executable file
    # EXEC_NAME="sub/exec/${SIMDIR}.sh"


    ## PARAMETERS - FILES
    # configuration files
    local SIM_CONFIG_EQUIL="config_equil.dat" # equilibrium configuration
    # contains end-to-end instructions
    local SIM_RE2E="RE2E.dat"

    ## PARAMETERS - SUBMISSION INTRUCTIONS
    # memory to request
    local REQUEST_MEMORY="500MB"
    # disk space to request
    local REQUEST_DISK="1GB"
    # directory that output files are remapped to
    local REMAP="node/run/"
    # list of files with remapping instructions
    local RMP_SIM_CONFIG_EQUIL="${SIM_CONFIG_EQUIL}=${REMAP}${SIM_CONFIG_EQUIL}"
    local RMP_SIM_RE2E="${SIM_RE2E}=${REMAP}${SIM_RE2E}"
    # list of files that should be transfered to the execute node
    local TRANSFER_INPUT_FILES="${SIM_CONFIG_EQUIL}, ${SIM_RE2E}"
    # list of files that should be transfered from the execute node
    local TRANSFER_OUTPUT_FILES="${SIM_CONFIG_EQUIL}, ${SIM_RE2E}"
    # list of remap instructions for each output file
    local TRANSFER_OUTPUT_REMAPS="${RMP_SIM_CONFIG_EQUIL}"


    ## ARGUMENTS
    # none


    ## SCRIPT
    # write submission script
    echo "executable = ${EXEC_NAME}" > $SUB_PATH
    echo "arguments = -r -p /LeMonADE_StiffCyclic/build/bin/" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "+SingularityImage = \"/home/mad/tmp/${DEFAULT_SINGULARITY}\"" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "should_transfer_files = YES" >> $SUB_PATH
    echo "transfer_input_files = ${TRANSFER_INPUT_FILES}" >> $SUB_PATH
    echo "transfer_output_files = ${TRANSFER_OUTPUT_FILES}" >> $SUB_PATH
    # echo "transfer_output_remaps = \"${TRANSFER_OUTPUT_REMAPS}\"" >> $SUB_PATH
    echo "when_to_transfer_output = ON_SUCCESS" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "log = out/run.log" >> $SUB_PATH
    echo "error = out/run.err" >> $SUB_PATH
    echo "output = out/run.out" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "request_cpus = 1" >> $SUB_PATH
    echo "request_disk = ${REQUEST_DISK}" >> $SUB_PATH
    echo "request_memory = ${REQUEST_MEMORY}" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "on_exit_hold = (ExitCode != 0)" >> $SUB_PATH
    echo "requirements = (HAS_GCC == true) && (Mips > 30000)" >> $SUB_PATH
    # echo "requirements = HasSingularity" >> $SUB_PATH
    echo "+ProjectName=\"NCSU_Hall\"" >> $SUB_PATH
    echo "" >> $SUB_PATH
    echo "queue" >> $SUB_PATH
}

## OPTIONS
while getopts "hzvstoj:p:f:l:e:c:" option; do
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
        o) # overwrite simulation state
            BOOL_OVERWRITE=1;;
        c) # specify checkfile
            declare -i BOOL_CHECKFILE=1
            OPTCHECKFILE="${OPTARG}";;
        z) # submit simulaitons on chtc
            declare -i BOOL_SUBMIT_CHTC=1 ;;
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
elif [[ $BOOL_SUBMIT_CHTC -eq 1 ]]; then 
    # else if submitting jobs on CHTC
    # generate dagman file
    DAG="${JOB}.dag"
    echo "" > $DAG # overwrite / initialize dag
    # loop through simulation parameters, add splice subdags into main dagman header
    declare -i N_LINES=$(wc -l < $PARMFILE)
    for i in $(seq 2 $N_LINES); do
        # get the simulation id (file column in file)
        SIMID=$(head -n ${i} ${PARMFILE} | tail -n 1 | cut -d , -f 1)
        # get the simulation directory (second column in file)
        SIMDIR=$(head -n ${i} ${PARMFILE} | tail -n 1 | cut -d , -f 2)
        # inform user
        if [[ $BOOL_VERB -eq 1 ]]; then 
            echo "Generating CHTC files for: ${SIMID} .."
        fi
        # write splce to dag
        echo -e "SPLICE ${SIMID} ${SIMID}.spl DIR ${JOBDIR}${SIMDIR}" >> $DAG
        gen_chtc_scripts
        # if test bool, just run one simulation
        if [[ $BOOL_TEST -eq 1 ]]; then
            exit 0
        fi
    done
    # TODO add job load maximum ..
else
    ## TODO :: compile executables locally, copy to local directory and delete after job finishes
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
