#!/bin/bash
set -e

## Matthew Dorsey
## 05.06.2024
## Programs for generating simulation parameters, directory hirarchy, and execution instructions


## PARAMETERS
## PARAMETERS -- CONSTANTS
PI_CON=$(echo "scale=10; 4*a(1)" | bc -l)
## PARAMETERS -- BOOLEAN
# boolean determining if the script should be executed verbosely
declare -i BOOL_VERB=0
# boolean determining if the script should generate the simulation parameters
declare -i BOOL_GEN=0
# boolean for ideal simulations
declare -i BOOL_IDEAL=0
## PARAMETERS -- JOB
# array containing N to test
PARM_N=( 200 )
# array containing whether to test rings structures or not
PARM_RING=( "FALSE" )
# array containing which potential to test
PARM_CSA=( "FALSE" "TRUE" )
# array containing persistance lengths to test for either potential
# assigned persistance length is used to select bending potential strength
# according to scaling equation unique to each potential'
PARM_LP=( "1" "1.2" "1.4" "1.6" "1.8" "2.0" "2.2" "2.4" "2.6" "2.8" "3" "3.2" "3.4" "3.6" "3.8" "4.0" "4.2" "4.4" "4.6" "4.8" "5" "7" "9" "11" "13" "15" "20" )
# default directory for upload / download, generating parameters
MAINDIR="01_raw_data"
# default job name
JOB="bendingPARM"
# path to LeMonADE executables
EXECDIR="00_programs/build/bin/"
## PARAMETERS -- SIMULATION
# number of MCSs between each property calculation
declare -i N_MCS=2000000000
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
	echo -e " -h           | display script options, exit 0."
	echo -e " -v           | execute script verbosely."
	echo -e " -g           | generate simulation parameters / directories ."
	echo -e " -i           | simulate ideal polymers (otherwise excluded volume interactions are included)."
	echo -e "\n"
	echo -e " ## SCRIPT PARAMETERS ##"
	echo -e " -j << ARG >> | specify job title (default is ${JOB})."
	echo -e " -p << ARG >> | specify the local directory (default is ${MAINDIR})."
	echo -e "\n"
}

# generate simulation parameters, write to file
gen_simparm() {

	## PARAMETERS
	# path to simulation directories
	PATH_SIMPARM="${MAINDIR}/${JOB}_${IDEAL_TAG}/"
	# file to write simulation parameters to
	FILE_SIMPARM="${PATH_SIMPARM}${JOB}_${IDEAL_TAG}.csv"
	# header used for the simulation parameter file
	HEADER_SIMPARM="id,path,pot,N,R,k"

	## ARGUMENTS
	# none

	## SCRIPT
	# if the path does not exist, make the directories
	if [ ! -d $PATH_SIMPARM ]; then
		mkdir $PATH_SIMPARM
	fi

	# write head to file, write parameters to file
	echo $HEADER_SIMPARM > $FILE_SIMPARM
	for r in ${PARM_RING[@]}; do
		for n in "${PARM_N[@]}"; do
			for c in ${PARM_CSA[@]}; do

				if [ "${c}" == "TRUE" ]; then
					C="CSA"
					C_FLAG=""
				else
					C="CA"
					C_FLAG="-c "
				fi

				for l in "${PARM_LP[@]}"; do
					# generate the simulation directory
					if [ "${r}" == "TRUE" ]; then
						SIMID="${C}_RING_N${n}LP${l}"
						SIMDIR="${C}/RING/N${n}/LP${l}/"
					else
						SIMID="${C}_CHAIN_N${n}LP${l}"
						SIMDIR="${C}/CHAIN/N${n}/LP${l}/"
					fi
					if [ "${c}" == "TRUE" ]; then
						C="CSA"
						C_FLAG=""
						K_STRING=$( echo "(${l} ^ 2) / ${PI_CON}" | bc -l )
					else
						C="CA"
						C_FLAG="-c "
						K_STRING="${l}"
					fi
					mkdir -p ${PATH_SIMPARM}${SIMDIR}
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
					if [ "${r}" == "TRUE" ]; then
						echo "\${PATH}generatePolymerBFM ${C_FLAG}-r -n ${n} -k ${K_STRING}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
					else
						echo "\${PATH}generatePolymerBFM ${C_FLAG}-n ${n} -k ${K_STRING}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
					fi
					echo "\${PATH}simulate${IDEAL_TAG}PolymerBFM -e ${t_equilibrium} -n ${N_MCS} -s ${save_interval} -q -d -c -a -g -m" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
					chmod u+x ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
					# add parameters to parm file
					echo "${SIMID},${SIMDIR},${C},${n},${r},${K_STRING}" >> $FILE_SIMPARM
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
		j) # update job name
			JOB="${OPTARG}" ;;
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
	# check for ideal polymer simulation
	if [[ $BOOL_IDEAL -eq 1 ]]; then
		IDEAL_TAG="ideal"
	else
		IDEAL_TAG="real"
	fi
	gen_simparm
fi
