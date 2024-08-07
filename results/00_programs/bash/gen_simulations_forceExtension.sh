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
# boolean determining if the script should generate force values along a log scale
declare -i BOOL_LOGSCALE=0
## PARAMETERS -- JOB
# number of unique forces to test
declare -i N_FORCE_VAL=40
# maxmium force to test
MAX_FORCE_VAL="2"
# minimum force to test
MIN_FORCE_VAL=".0001"
# array containing N to test
PARM_N=( 100 )
# array containing whether to test rings structures or not
PARM_RING=( 0 1 2 ) # 2
# array containing which potential to test
PARM_CSA=( "TRUE" )
# array containing bending potential strings to test
PARM_BEND=( 0 1 5 10 30 )
# array containing different force vectors to test
PARM_FORCEVEC=( "100" ) # "100"
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
declare -i t_equilibrium=500000000


## FUNCTIONS
# display script options
help () {
	echo -e "\nScript for generating polymer simulations with BFM model on linux clusters.\nUSAGE: ./gen_scaling simulations.sh << FLAGS >>\n"
	echo -e "\n"
	echo -e " ## SCRIPT PROTOCOL ##"
	echo -e " -g           | generate simulation parameters / directories ."
	echo -e " -l           | generate simulation parameters according to logscale (no negatives) ."
	echo -e " -n           | number of simulation parameters to test (default is ${N_FORCE_VAL})."
	echo -e " -m           | maximum force to test (default is ${MAX_FORCE_VAL})."
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

# used to generate parameters along a log scale
linscale() {
    ## PARAMETERS
    # minimum number on a linear scale
    MIN_VAL_LIN="-${MAX_FORCE_VAL}"
    # maximum value along a linear scale
    MAX_VAL_LIN=${MAX_FORCE_VAL}

    ## ARGUMENTS
    # integer ranging from 1 to N_FORCE_VAL
    declare -i NUM=$1

    ## SCRIPT
    # generate the parameter along the linear scale
    scale=$( echo "($NUM / ( ${N_FORCE_VAL} - 1 ))" | bc -l )
    scale=$( echo "(${scale} * (${MAX_VAL_LIN} - ${MIN_VAL_LIN}) + ${MIN_VAL_LIN})" | bc -l )
    echo $(printf "%8.5f\n" "${scale}")
}

# generate simulation parameters
gen_simparm() {

	## PARAMETERS
	# path to simulation directories
	PATH_SIMPARM="${MAINDIR}/${JOB}/"
	# file to write simulation parameters to
	FILE_SIMPARM="${PATH_SIMPARM}${JOB}.csv"
	# header used for the simulation parameter file
	HEADER_SIMPARM="id,path,pot,R,N,K,F,FV"

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
                        for fv in "${PARM_FORCEVEC[@]}"; do

                            # generate the simulation directory
                            SIMID="${R}_${C}_N${n}K${k}FV${fv}F${f}"
                            SIMDIR="${R}/${C}/N${n}/K${k}/FV${fv}/F${f}/"
                            mkdir -p ${PATH_SIMPARM}${SIMDIR}

                            # use the force integer to calculate the real force value passed to the simulation executable
                            if [[ $LOGSCALE -eq 1 ]]; then
                                FORCE_VAL=$( logscale $f )
                            else
                                FORCE_VAL=$( linscale $f )
                            fi

                            # generate flags for simulation executables
                            GENFLAGS="-n ${n} -f ${FORCE_VAL} -v ${fv} -b 512"
                            SIMFLAGS="-e ${t_equilibrium} -n ${N_MCS} -s ${save_interval} -a -d"
                            if [ $k != 0 ]; then
                                GENFLAGS="${GENFLAGS} -k ${k}"
                                if [ "${C}" == "CA" ]; then
                                    GENFLAGS="${GENFLAGS} -c"
                                fi
                            fi

                            if [ $r -eq 1 ]; then
                                GENFLAGS="${GENFLAGS} -r"
                            elif [ $r -eq 2 ]; then
                                GENFLAGS="${GENFLAGS} -r -m 2"
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
                            echo "${SIMID},${SIMDIR},${C},${r},${n},${k},${FORCE_VAL},${fv}" >> $FILE_SIMPARM
                        done
                    done
                done
            done
        done
	done
}

## OPTIONS
while getopts "hgln:m:" option; do
	case $option in
		h) # print script parameters to CLT
			help
			exit 0 ;;
        g) # generate simulation parameters
            declare -i BOOL_GEN=1 ;;
		l) # generate parameters along logscale
			declare -i BOOL_LOGSCALE=1 ;;
        n) # number of force extension values to test
            declare -i N_FORCE_VAL=${OPTARG};;
        m) # maxmimum force value to test
            MAX_FORCE_VAL="${OPTARG}";;
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
