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
# overwrite / delete all existing simulations
declare -i BOOL_OVER=0
## PARAMETERS -- JOB
# array containing N to test
PARM_N=( 100 )
# array containing whether to test rings structures or not
# PARM_RING=( 0 1 2 )
PARM_RING=( 0 )
# array containing which potential to test
PARM_CSA=( "TRUE" )
# array containing persistance lengths to test for either potential
# assigned persistance length is used to select bending potential strength
# according to scaling equation unique to each potential'
PARM_LP=( 0 2 5 9 )
# number of points sampled by hysteresis calculator
# NOTE :: this is a fixed value in the hysteresis analyzer
# TODO :: integrate number of sampling points in the script and the number of the executable
declare -i N_SAMPLE=500
# number of times that the hysteresis loops through one period once equilibrated
declare -i N_PERIOD_EQUILIBRIUM=45
# number of times that the hysteresis loops through one period to reach equilibrium
declare -i N_PERIOD_EQUILIBRIATE=5
# number of different period points to test
declare -i N_PERIOD_VAL=40
# maximum period value to test
declare -i MAX_PERIOD_VAL=100000000
# minimum period value to test
declare -i MIN_PERIOD_VAL=10000
# force base value
VAL_FO="0.0"
# force amplitude
VAL_FA="0.5"
# force vector applied by oscillatory force
FORCEVEC="100"
## PARAMETERS -- SIMULATION
# default directory for upload / download, generating parameters
MAINDIR="01_raw_data"
# default job name
JOB="hysteresis"
# path to LeMonADE executables
EXECDIR="00_programs/build/bin/"


## FUNCTIONS
# display script options
help () {

	echo -e "\nScript for generating polymer simulations with BFM model on linux clusters.\nUSAGE: ./gen_scaling simulations.sh << FLAGS >>\n"
	echo -e "\n"
	echo -e " ## SCRIPT PROTOCOL ##"
	echo -e " -h           | display script options, exit 0."
	echo -e " -v           | execute script verbosely."
	echo -e " -g           | generate simulation parameters / directories ."
	echo -e " -o           | overwrite existing simulations, restart parameter file."
	echo -e "\n"
	echo -e " ## SCRIPT PARAMETERS ##"
	echo -e " -j << ARG >> | specify job title (default is ${JOB})."
	echo -e " -p << ARG >> | specify the local directory (default is ${MAINDIR})."
	## TODO :: add option for ring type
	## TODO :: add option for bending potential
	## TODO :: add option for base force
	## TODO :: add option for force amplitude
	## TODO :: add option for persistence length
	## TODO :: add option for number of sampling points
	## TODO :: add option for min period
	## TODO :: add option for max period
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
     MIN_VAL_LOG10=$( log10 ${MIN_PERIOD_VAL} )
     # maximum number on a log10 scale
     MAX_VAL_LOG10=$( log10 ${MAX_PERIOD_VAL} )

     ## ARGUMENTS
     # integer, ranging from 1 to N_FORCE_VAL
     declare -i NUM=$1

     ## SCRIPT
     # generate the parameter along scale
     scale=$( echo "($NUM / ( ${N_PERIOD_VAL} - 1 ))" | bc -l )
     scale=$( echo "(${scale} * (${MAX_VAL_LOG10} - ${MIN_VAL_LOG10}) + ${MIN_VAL_LOG10})" | bc -l )
     scale=$( pow10 "$scale" )
     echo $(printf "%8.0f\n" "${scale}")
}

# method for rounding any whole number
# to a number divisable by $N_SAMPLE
round () {

	## PARAMETERS
	# none

	## ARGUMENTS
	# first argument: integer to round
	declare -i rnd=$1
	# second argument: whole number that the rounded number should be divisiable by
	declare -i wn=$2


	## SCRIPT
	# use mod to determine the remainder
	declare -i rnd_dn=$(bc -l <<< "oldscale=scale; scale=0; ${rnd}%${wn}; scale=oldscale") # remained to remove if rounding down
	declare -i rnd_up=$( echo "(${wn} - ${rnd_dn})" | bc -l) # remainder to add if rounding up
	# determing wether to round up or down
	if (( rnd_dn > rnd_up )); then
		# if rnd_dn is greater than rnd_up
		# rnd is closer to the round up value so round up
		declare -i wn=$( echo "(${rnd} + ${rnd_up})" | bc -l )
	else
		# otherwise, round down
		declare -i wn=$( echo "(${rnd} - ${rnd_dn})" | bc -l )
	fi
	# return to user
	echo $wn

}

# generate simulation parameters, write to file
gen_simparm() {

	## PARAMETERS
	# path to simulation directories
	PATH_SIMPARM="${MAINDIR}/${JOB}/"
	# file to write simulation parameters to
	FILE_SIMPARM="${PATH_SIMPARM}${JOB}.csv"
	# header used for the simulation parameter file
	HEADER_SIMPARM="id,path,pot,N,R,k,fo,fA,T,fv"

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

		if [ $r -eq 0 ]; then
			# chain
			R="CHAIN"
			R_FLAG=""
		elif [ $r -eq 1 ]; then
			# single ring
			R="RING"
			R_FLAG=" -r"
		elif [ $r -eq 2 ]; then
			R="RINGx2"
			R_FLAG=" -r -m 2"
		fi

		for n in "${PARM_N[@]}"; do

# 			if [ $r != 0 ]; then
# 				n=$( echo "${n} * 2" | bc -l) # double number of monomers for ring so same length as chain
# 			fi

			for c in ${PARM_CSA[@]}; do

				if [ "${c}" == "TRUE" ]; then
					C="CSA"
					C_FLAG=""
				else
					C="CA"
					C_FLAG=" -c"
				fi

				for l in "${PARM_LP[@]}"; do

					# use the persistence length to determine the bending constant according to the potential
					if [ "${c}" == "TRUE" ]; then
						K_STRING=$( echo $( printf "%0.4f\n" $( echo "(${l} ^ 2) / ${PI_CON}" | bc -l )))
					else
						K_STRING="${l}"
					fi
					K_FLAG=" -k ${K_STRING}"

                    for p in $(seq 0 $(($N_PERIOD_VAL-1))); do

						# loop through points, generate on a log scale
						PERIOD_VAL=$( logscale $p )
						# round the period value so that it integrates with the sampling frequency
						PERIOD_VAL=$( round $PERIOD_VAL $N_SAMPLE )

						# generate the simulation directory
						SIMDIR="${R}/${C}/N${n}/LP${l}/FO${VAL_FO}/FA${VAL_FA}/T${PERIOD_VAL}/FV${FORCEVEC}/"
						SIMID="${R}${C}_N${n}LP${l}_FV${FORCEVEC}FO${VAL_FO}FA${VAL_FA}T${PERIOD_VAL}"

						# adjust the length of the simulation according the period
						declare -i TOTAL=$( echo "(${PERIOD_VAL} * (${N_PERIOD_EQUILIBRIATE} + ${N_PERIOD_EQUILIBRIUM}))" | bc -l )
						declare -i EQUILIBRIATE=$( echo "(${PERIOD_VAL} * ${N_PERIOD_EQUILIBRIATE})" | bc -l )

						# generate flags for simulation executables
						GENFLAGS="-n ${n} -v ${FORCEVEC} -f ${VAL_FO} -a ${VAL_FA} -p ${PERIOD_VAL} -b 512 ${R_FLAG}${K_FLAG}${C_FLAG}"
						SIMFLAGS="-e ${EQUILIBRIATE} -n ${TOTAL} -b"
						# TODO :: add sample frequency to execute flag

						# generate directory, write files with parameters
						mkdir -p ${PATH_SIMPARM}${SIMDIR}
						# write files directory, generate simulation executables
						echo "#!/bin/bash" > ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo "set -e" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo "PATH=\"./\"" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo "while getopts \"p:\" option; do" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo -e "\tcase \$option in" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo -e "\t\tp)" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo -e "\t\t\tPATH=\${OPTARG}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo -e "\tesac" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo "done" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo "\${PATH}generatePolymerBFM ${GENFLAGS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						echo "\${PATH}simulateRealPolymerBFM ${SIMFLAGS}" >> ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						chmod u+x ${PATH_SIMPARM}${SIMDIR}${SIMID}.sh
						# add parameters to parm file
						# "id,path,pot,N,R,k,fo,fA,T,fv"
						echo "${SIMID},${SIMDIR},${C},${n},${r},${K_STRING},${VAL_FO},${VAL_FA},${PERIOD_VAL},${FORCEVEC}" >> $FILE_SIMPARM

                    done
				done
			done
		done
	done
}


## OPTIONS
while getopts "hvgosudj:l:p:" option; do
	case $option in
		h) # print script parameters to CLT
			help
			exit 0 ;;
        v) # exectue script verbosely
            declare -i BOOL_VERB=1 ;;
        g) # generate simulation parameters
            declare -i BOOL_GEN=1 ;;
		o) # restart job accumulation
			declare -i BOOL_OVER=1 ;;
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
## TODO :: implement overwrite feature (this rewrite the commands for any directoraies that already exist)
## TODO :: add restart feature (scorched earth, delete everything and restart)
# generate parameters save to file
if [ $BOOL_GEN -eq 1 ]; then
	gen_simparm
fi
