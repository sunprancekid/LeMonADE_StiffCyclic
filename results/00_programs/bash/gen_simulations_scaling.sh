#!/bin/bash
set -e

## Matthew Dorsey
## 05.06.2024
## Programs for generating simulation parameters, directory hirarchy, and execution instructions


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
# array containing N to test
PARM_N=( 32 48 64 78 96 128 256 )
# array containing the different executables to test
PARM_EX=( "linearChainReal" "linearChainIdeal" )
# number of times to replicate each unique set of simualation conditions
declare -i PARM_R=3
# default directory for upload / download, generating parameters
MAINDIR="01_raw_data"
# default job name
JOB="scaling"
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

# generate simulation parameters, write to file
gen_simparm() {

	## PARAMETERS
	# path to simulation directories
	PATH_SIMPARM="${MAINDIR}/${JOB}/"
	# file to write simulation parameters to
	FILE_SIMPARM="${PATH_SIMPARM}${JOB}.csv"
	# header used for the simulation parameter file
	HEADER_SIMPARM="id,mod,N,R,path"

	## ARGUMENTS
	# none

	## SCRIPT
	# if the path does not exist, make the directories
	if [ ! -d $PATH_SIMPARM ]; then
		mkdir $PATH_SIMPARM
	fi

	# write head to file, write parameters to file
	echo $HEADER_SIMPARM > $FILE_SIMPARM
	for m in "${PARM_EX[@]}"; do
		for n in "${PARM_N[@]}"; do
			for r in $(seq 1 $PARM_R); do
				# generate the simulation directory
				SIMID="${m}_N${n}R${r}"
				SIMDIR="${m}/N${n}/R${r}/"
				mkdir -p ${PATH_SIMPARM}${SIMDIR}
				# write files directory, generate simulation executables
				cp $EXECDIR${m} ${PATH_SIMPARM}${SIMDIR}
				echo "./$m ${n} ${N_MCS} ${N_runs} ${t_equilibrium} > ${SIMID}.txt" > ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
				chmod u+x ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
				# add parameters to parm file
				echo "${SIMID},${m},${n},${r},${SIMDIR}" >> $FILE_SIMPARM
			done
		done
	done
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
## TODO add overwrite feature
# upload directory to linux cluster

# download directory from linux cluster

# generate parameters save to file
if [ $BOOL_GEN -eq 1 ]; then
	gen_simparm
fi

# parse simulation parameters from file, perform protocol as instructed
SIMPARAM_FILE="${MAINDIR}/${JOB}/${JOB}.csv"
if [ ! -f $SIMPARAM_FILE ]; then
	# if the file does not exist, inform user and abort
	echo "Must generate simulation parameters for ${JOB} in ${MAINDIR} before submitting simulations."
	exit $NONZERO_EXITCODE
fi
declare -i N_LINES=$(wc -l < $SIMPARAM_FILE)
for i in $(seq 2 $N_LINES); do
	# get the simulation directory
	SIMDIR=$(head -n ${i} ${SIMPARAM_FILE} | tail -n 1 | cut -d , -f 5)
	# get the simulation id
	SIMID=$(head -n ${i} ${SIMPARAM_FILE} | tail -n 1 | cut -d , -f 1)
	if [ ! -f ${MAINDIR}/${JOB}/${SIMDIR}/RE2E.dat ]; then
		# move to the simulation directory
		CURRDIR=$(echo $PWD)
		cd "${MAINDIR}/${JOB}/${SIMDIR}"
		echo $PWD
		# exectue the simulation
		if [ $BOOL_SUB -eq 1 ]; then
			./${SIMID}.sh
		fi
		# return to the main directory
		cd $CURRDIR
	fi
done
