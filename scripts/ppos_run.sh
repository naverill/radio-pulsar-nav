#!/bin/bash
# alias tempo2 "/usr/local/tempo2/bin/tempo2"
# setenv TEMPO2 /usr/local/tempo2

INPUT_SIM_DIR=~/workspace/rpnav/rpnav/simulate/inputs
OUTPUT_SIM_DIR=~/workspace/rpnav/rpnav/simulate/outputs
SIM_SCRIPT=navigate_7d
SIM_DIR=${OUTPUT_SIM_DIR}/${SIM_SCRIPT}_sim
SIM_RESULTS_DIR=${SIM_DIR}/output/real_0
runscript=${SIM_DIR}/scripts/runScripts_master

export PSR1="J0835-4510"
export PSR2="J1939+2134"

cd $OUTPUT_SIM_DIR

for ((i = 44; i < 100; i++));
do
    rm -rf $SIM_RESULTS_DIR

    # Rerun simulation  
    /home/naverill/external/ptasimulate/ptaSimulate ${INPUT_SIM_DIR}/${SIM_SCRIPT}.input
    

    # Insert aliasing tempo2 at the start of the script
    echo "alias tempo2 /usr/local/tempo2/bin/tempo2
setenv TEMPO2 /usr/local/tempo2
" | \
        cat - ${runscript} >> temp \
        && mv temp ${runscript}

    tcsh ${runscript}

    # Transform simulated par files into TDB
    /usr/local/tempo2/bin/tempo2 -gr transform ${SIM_RESULTS_DIR}/${PSR1}.par ${SIM_RESULTS_DIR}/${PSR1}.tdb.par tdb
    /usr/local/tempo2/bin/tempo2 -gr transform ${SIM_RESULTS_DIR}/${PSR2}.par ${SIM_RESULTS_DIR}/${PSR2}.tdb.par tdb


    OUT_RESULTS_DIR=$(printf "${OUTPUT_SIM_DIR}/${SIM_SCRIPT}_sim/output/S%02d" $i)
    mkdir -p ${OUT_RESULTS_DIR}
    simFile=$(printf "results_S%02d.csv" $i)
    echo "Iteration,Longitude(deg),Latitude(deg),X(m),Y(m),Z(m),Error,Step Size" > ${outResDir}/${simFile}
    for ((j = 0; j < 100; j++));
    do
        resFile=$(printf "results_S%02dI%02d.csv"  $i $j)

        /usr/local/tempo2/bin/tempo2 -gr pulsar_positioning \
            -f ${SIM_RESULTS_DIR}/${PSR1}.tdb.par ${SIM_RESULTS_DIR}/${PSR1}.tim \
            -f ${SIM_RESULTS_DIR}/${PSR2}.tdb.par ${SIM_RESULTS_DIR}/${PSR2}.tim \
            -a grde -e chi -n -r -s ${OUT_RESULTS_DIR}/${resFile}
            
        tail -n 1 ${OUT_RESULTS_DIR}/${resFile} >> ${OUT_RESULTS_DIR}/${simFile}
    done
done