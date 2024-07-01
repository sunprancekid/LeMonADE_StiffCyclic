#!/bin/bash
set -e

## Matthew Dorsey
## 2024.06.11
##  script for generating constant force extension simulations for real and ideal polymer chains



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
# number of unique forces to test
declare -i N_FORCE_VAL=50
# maxmium force to test
MAX_FORCE_VAL="100"
# minimum force to test
MIN_FORCE_VAL=".0001"
# array containing N to test
PARM_N=( 100 )
# array containing whether to test rings structures or not
PARM_RING=( 0 1 ) # TODO add bonded rings
# array containing which potential to test
PARM_CSA=( "FALSE" )
# array containing bending potential strings to test
PARM_BEND=( 0 7 13 )
# number of times to replicate each unique set of simualation conditions
declare -i PARM_R=1
# default directory for upload / download, generating parameters
MAINDIR="01_raw_data"
# default job name
JOB="forceExtension"
# linux server
LINUXSERV="gandalf"
# path to LeMonADE executables
EXECDIR="00_programs/build/bin/"
## PARAMETERS -- SIMULATION
# number of MCSs between each property calculation
declare -i N_MCS=1000000000
# number of time properties are calculation
declare -i save_interval=1000000
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

# log10 function, echos log10 of first argument passed to method
log10 (){

    ## PARAMETERS
    # none

    ## ARGUMENTS
    # number to perform log10 operation on
    local NUM_LOG=$1

    ## SCRIPT
    # perform log10 operation on number
    echo "l(${NUM_LOG})/l(10)" | bc -l
}

# pow10 function, echos 10 to the power of the first argument passed to method
pow10 () {

    ## PARAMETERS
    # none

    ## ARGUMENTS
    # number to perform pow10 operation on
    local NUM_POW=$1

    ## SCRIPT
    # perform pow10 operation
    VAL=$( echo "l(10)" | bc -l )
    VAL=$( echo "(${VAL})*(${NUM_POW})" | bc -l )
    echo "e(${VAL})" | bc -l
}

# used to generate parameters along a logscale
logscale () {

     ## PARAMETERS
     # minimum number on a log10 scale
     MIN_VAL_LOG10=$( log10 ${MIN_FORCE_VAL} )
     # maximum number on a log10 scale
     MAX_VAL_LOG10=$( log10 ${MAX_FORCE_VAL} )

     ## ARGUMENTS
     # integer, ranging from 1 to N_FORCE_VAL
     declare -i NUM=$1

     ## SCRIPT
     # generate the parameter along scale
     scale=$( echo "($NUM / ( ${N_FORCE_VAL} - 1 ))" | bc -l )
     scale=$( echo "(${scale} * (${MAX_VAL_LOG10} - ${MIN_VAL_LOG10}) + ${MIN_VAL_LOG10})" | bc -l )
     scale=$( pow10 "$scale" )
     echo $(printf "%8.5f\n" "${scale}")
}

# generate slurm submission script for simulation
gen_slurm () {

    ## PARAMETERS
    # name of file contains slurm submission instructions
    local FILENAME="${SIMID}.slurm.sub"
    # directory that the file is stored in
    local FILEPATH="${PATH_SIMPARM}${SIMDIR}"

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
    echo "#SBATCH --output=OneLinearChainIdeal_N100_PerX512_PerYZ128_ConstantForce_f0_0.01.%j.out" >> $FILEPATH$FILENAME
    # echo "#SBATCH --mail-type=FAIL" >> $FILEPATH$FILENAME
    # echo "#SBATCH --mail-user=dorsey@ipfdd.de" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### PARAMETERS ###" >> $FILEPATH$FILENAME
    echo "SUBDIR=\$(pwd) # location of slurm submission" >> $FILEPATH$FILENAME
    echo "EXECDIR=/beetmp/dorsey/tmp/ # location of slurm execution" >> $FILEPATH$FILENAME
    echo "hostname" >> $FILEPATH$FILENAME
    echo "cd \${EXECDIR}" >> $FILEPATH$FILENAME
    echo "mkdir ${SIMID} " >> $FILEPATH$FILENAME
    echo "cd ${SIMID}" >> $FILEPATH$FILENAME
    echo "echo Running on host \$(hostname)" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### load modules" >> $FILEPATH$FILENAME
    echo "# module load slurm/15.08.8" >> $FILEPATH$FILENAME
    echo "# module load gcc/6.1.0" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### copy everything from the submit directory to the execute directory ###" >> $FILEPATH$FILENAME
    echo "cp \${SUBDIR}/* ." >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    echo "### run job ####" >> $FILEPATH$FILENAME
    echo "srun ./${SIMID}.sh > ${SIMID}.wrap.out 2>&1" >> $FILEPATH$FILENAME
    echo "" >> $FILEPATH$FILENAME
    ## TODO :: add instructions for copying results back
    echo "### copy back results, delete everthing ###" >> $FILEPATH$FILENAME
    echo "cp * \${SUBDIR}/" >> $FILEPATH$FILENAME
    echo "cd ../" >> $FILEPATH$FILENAME
    echo "rm -rf ${SIMID}" >> $FILEPATH$FILENAME
}

# generate simulation parameters
gen_simparm() {

	## PARAMETERS
	# path to simulation directories
	PATH_SIMPARM="${MAINDIR}/${JOB}/"
	# file to write simulation parameters to
	FILE_SIMPARM="${PATH_SIMPARM}${JOB}.csv"
	# header used for the simulation parameter file
	HEADER_SIMPARM="id,path,pot,R,N,K,F"

	## ARGUMENTS
	# none

	## SCRIPT
	# if the path does not exist, make the directories
	if [ ! -d $PATH_SIMPARM ]; then
		mkdir $PATH_SIMPARM
	fi

	# write head to file, write parameters to file
	echo $HEADER_SIMPARM > $FILE_SIMPARM
	for c in "${PARM_CSA[@]}"; do

        if [ "${c}" == "TRUE" ]; then
            C="CSA"
            C_FLAG=""
        else
            C="CA"
            C_FLAG="-c "
        fi

        for r in "${PARM_RING[@]}"; do
            if [ $r -eq 0 ]; then
                # chain
                R="CHAIN"
            elif [ $r -eq 1 ]; then
                # single ring
                R="RING"
            elif [ $r -eq 2 ]; then
                R="RINGx2"
            fi
            for n in "${PARM_N[@]}"; do

                if [ $r != 0 ]; then
                    n=$( echo "${n} * 2" | bc -l) # double number of monomers for ring so same length as chain
                fi

                for k in "${PARM_BEND[@]}"; do
                    for f in $(seq 0 $(($N_FORCE_VAL-1))); do
                        # generate the simulation directory
                        SIMID="${R}_${C}_N${n}K${k}F${f}"
                        SIMDIR="${R}/${C}/N${n}/K${k}/F${f}/"
                        mkdir -p ${PATH_SIMPARM}${SIMDIR}
                        # use the force integer to calculate the real force value passed to the simulation executable
                        FORCE_VAL=$( logscale $f )
                        # generate flags for simulation executables
                        GENFLAGS="-n ${n} -f ${FORCE_VAL} -b 512"
                        SIMFLAGS="-e ${t_equilibrium} -n ${N_MCS} -s ${save_interval}"
                        if [ $k != 0 ]; then
                            GENFLAGS="${GENFLAGS} -k ${k}"
                            if [ "${C}" == "CA" ]; then
                                GENFLAGS="${GENFLAGS} -c"
                            fi
                        fi
                        if [ $r -eq 1 ]; then
                            GENFLAGS="${GENFLAGS} -r"
                        elif [ $r -eq 2 ]; then
                            GENFLAGS="${GENFLAGS} -r -m 1"
                        fi
                        # write files directory, generate simulation executables
                        echo "#!/bin/bash" > ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo "set -e" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo "PATH=\"./\"" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo "while getopts \"p:\" option; do" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo -e "\tcase \$option in " >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo -e "\t\tp)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo -e "\t\t\tPATH=\${OPTARG}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo -e "\tesac" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo "done" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo "\${PATH}generatePolymerBFM ${GENFLAGS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        echo "\${PATH}simulatePolymerBFM ${SIMFLAGS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        chmod u+x ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                        # add parameters to parm file
                        echo "${SIMID},${SIMDIR},${C},${r},${n},${k},${FORCE_VAL}" >> $FILE_SIMPARM
                    done
                done
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
	SIMDIR=$(head -n ${i} ${SIMPARAM_FILE} | tail -n 1 | cut -d , -f 2)
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
