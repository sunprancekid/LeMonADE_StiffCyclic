#!/bin/bash
set -e

## Matthew Dorsey
## 2024.06.11
## script for generating constant force extension simulations for real and ideal polymer chains



## PARAMETERS
## PARAMETERS -- CONSTANTS
PI_CON=$(echo "scale=10; 4*a(1)" | bc -l)
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
# boolean for ideal simulations
declare -i BOOL_IDEAL=0
## PARAMETERS -- JOB
# number of unique forces to test
declare -i N_FORCE_VAL=200
# maxmium force to test
MAX_FORCE_VAL="10"
# minimum force to test
MIN_FORCE_VAL=".001"
# array containing N to test
PARM_N=( 100 )
# array containing whether to test rings structures or not
PARM_RING=( 0 1 2 )
# array containing which potential to test
PARM_CSA=( "TRUE" )
# array containing persistence lengths to test
PARM_LP=( 0 2 5 9 )
# array containing bending potential strings to test
# PARM_BEND=( 0 1 5 10 30 )
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
declare -i N_MCS=10000000000 # 10 billion
# number of time properties are calculation
declare -i save_interval=1000000
# number of MCSs before equilibrium properties are calculated
declare -i t_equilibrium=1000000000 # 1 billion


## FUNCTIONS
# display script options
help () {
	echo -e "\nScript for generating polymer simulations with BFM model on linux clusters.\nUSAGE: ./gen_scaling simulations.sh << FLAGS >>"
	echo -e "\n"
	echo -e " ## SCRIPT PROTOCOL ##"
	echo -e " -g           | generate simulation parameters / directories ."
	echo -e " -l           | generate simulation parameters according to logscale (no negatives)."
    echo -e " -i           | simulate ideal chain (deafult is real)."
	echo -e " -j << ARG >> | rename job."
	echo -e " -n << ARG >> | number of simulation parameters to test (default is ${N_FORCE_VAL})."
	echo -e " -m << ARG >> | maximum force to test (default is ${MAX_FORCE_VAL})."
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
     scale=$( echo "(($NUM - 1) / ( ${N_FORCE_VAL} - 1 ))" | bc -l )
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
    scale=$( echo "(($NUM - 1 ) / ( ${N_FORCE_VAL} - 1 ))" | bc -l )
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
	HEADER_SIMPARM="id,path,pot,mcs_run,mcs_equil,mcs_si,R,N,K,F,FV"

	## ARGUMENTS
	# none

	## SCRIPT
	# if the path does not exist, make the directories
	if [ ! -d $PATH_SIMPARM ]; then
		mkdir $PATH_SIMPARM
	fi

    # determine flag used for simulating either a real or an ideal chain
    EXCLUDED_VOLUME_FLAG="Real"
    if [[ $BOOL_IDEAL -eq 1 ]]; then
        EXCLUDED_VOLUME_FLAG="Ideal"
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
                for l in "${PARM_LP[@]}"; do
                    for f in $(seq 0 $(($N_FORCE_VAL))); do
                        for fv in "${PARM_FORCEVEC[@]}"; do

                            # generate the simulation directory
                            SIMID="${R}_${C}_N${n}LP${l}FV${fv}F${f}${log_tag}"
                            SIMDIR="${R}/${C}/N${n}/LP${l}/FV${fv}/F${f}/"
                            mkdir -p ${PATH_SIMPARM}${SIMDIR}

                            # for the first incrementation of the force val, get the equilibrium end to end distance (i.e. no force)
                            FORCE_FLAG=" "
                            FORCE_VAL="0"
                            if [[ $f -gt 0 ]]; then
                                # use the force integer to calculate the real force value passed to the simulation executable
                                if [[ $BOOL_LOGSCALE -eq 1 ]]; then
                                    FORCE_VAL=$( logscale $f )
                                else
                                    FORCE_VAL=$( linscale $f )
                                fi
                                FORCE_FLAG=" -f ${FORCE_VAL} "
                            fi

                            # generate flags for simulation executables
                            GENFLAGS="-n ${n}${FORCE_FLAG}-v ${fv} -b 512"
                            SIMFLAGS="-e ${t_equilibrium} -n ${N_MCS} -s ${save_interval} -a"
                            BENDFLAG=""
                            if [ $l != 0 ]; then
                                # calculate the bending parameter constant from the assigned persistence length
                                # according to the type of potential assigned to the simulation
                                if [ "${C}" == "CA" ]; then
                                    BENDFLAG=" -c"
                                    k="${l}"
                                else
                                    k=$( echo "(${l} ^ 2) / ${PI_CON}" | bc -l )
                                fi
                                BENDFLAG="${BENDFLAG} -k ${k}"
                            else
                                k="0."
                            fi

                            # topology flag
                            TOPFLAG=""
                            if [ $r -eq 1 ]; then
                                TOPFLAG="-r"
                            elif [ $r -eq 2 ]; then
                                TOPFLAG="-r -m 2"
                            fi

                            ## TODO :: add loop int

                            # write files directory, generate simulation executables
                            echo "#!/bin/bash" > ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "set -e" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "PATH=\"./\"" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "CHECKFILE_EQUIL=\"config_equil.dat\"" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "CHECKFILE_GEN=\"config_gen.dat\"" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i GEN_BOOL=0" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i EQUIL_BOOL=0" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i RUN_BOOL=0" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i MAX_MCS_BOOL=0" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i SAVE_MCS=${save_interval}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i EQUIL_MCS=${t_equilibrium}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i RUN_MCS=${N_MCS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "declare -i MAX_MCS=100000000000" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo "while getopts \"gerp:n:\" option; do" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tcase \$option in " >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tp)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tPATH=\${OPTARG};;" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tg)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i GEN_BOOL=1;;" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\te)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i EQUIL_BOOL=1;;" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tr)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i RUN_BOOL=1;;" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tn)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i MAX_MCS_BOOL=1" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i MAX_MCS=\${OPTARG};;" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tesac" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "done" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "if [ \$MAX_MCS_BOOL -eq 1 ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tif [ \$EQUIL_BOOL -eq 1 ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tif [ \$EQUIL_MCS -gt \$MAX_MCS ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i EQUIL_MCS=\$MAX_MCS" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tif [ \$RUN_BOOL -eq 1 ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tif [ \$RUN_MCS -gt \$MAX_MCS ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\tdeclare -i RUN_MCS=\$MAX_MCS" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "fi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "if [ \$GEN_BOOL -eq 1 ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tif [ ! -f \$CHECKFILE_GEN ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t# if the inital state does not exist, create it" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\${PATH}generatePolymerBFM ${GENFLAGS} -o \${CHECKFILE_GEN} ${BENDFLAG} ${TOPFLAG}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "fi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "if [ \$EQUIL_BOOL -eq 1 ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tif [ ! -f \${CHECKFILE_EQUIL} ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t# if the equilibrium state does not exist, create it" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\${PATH}simulate${EXCLUDED_VOLUME_FLAG}PolymerBFM -e \${EQUIL_MCS} -n \${EQUIL_MCS} -s \${EQUIL_MCS} -a -m -f \${CHECKFILE_GEN} -o \${CHECKFILE_EQUIL}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "fi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "if [ \$RUN_BOOL -eq 1 ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tif [ ! -f \${CHECKFILE_EQUIL} ]; then" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t# if the equilibrium file does not exist, there is an error" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\techo \"\${CHECKFILE_EQUIL} does not exist on run call.\" 1>&2" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\texit 1" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\telse" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t# iterate through run calls, accumulate properties" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\t\t\${PATH}simulate${EXCLUDED_VOLUME_FLAG}PolymerBFM -n \${RUN_MCS} -s \${SAVE_MCS} -a -m -f \${CHECKFILE_EQUIL} -o \${CHECKFILE_EQUIL}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "\tfi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            echo -e "fi" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            # echo "\${PATH}generatePolymerBFM ${GENFLAGS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            # echo "\${PATH}simulateRealPolymerBFM ${SIMFLAGS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            chmod u+x ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
                            # add parameters to parm file
                            echo "${SIMID},${SIMDIR},${C},${N_MCS},${t_equilibrium},${save_interval},${r},${n},${k},${FORCE_VAL},${fv}" >> $FILE_SIMPARM
                        done
                    done
                done
            done
        done
	done
}

## OPTIONS
while getopts "hglij:n:m:" option; do
	case $option in
		h) # print script parameters to CLT
			help
			exit 0 ;;
        j) # change job name
            JOB="${OPTARG}";;
        g) # generate simulation parameters
            declare -i BOOL_GEN=1 ;;
		l) # generate parameters along logscale
			declare -i BOOL_LOGSCALE=1 ;;
        i) # simulate ideal chain
            declare -i BOOL_IDEAL=1;;
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

    # if the subdirectory does not exist, create it
    if [[ ! -d "01_raw_data" ]]; then
        mkdir "01_raw_data"
    fi

	# differentiate jobs with lin and log tags
	if [[ $BOOL_LOGSCALE -eq 1 ]]; then
        log_tag="log"
    else
        log_tag="lin"
	fi
	# tag job
	JOB="${JOB}_${log_tag}"
	gen_simparm
fi
