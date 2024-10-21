#!/bin/bash

#************ CONFIGURATION PARAMETERS *************
# Location of the directory containing the .input file
INPUT_SIM_DIR=~/workspace/rpnav/rpnav/simulate/inputs

# Location of the directory to save the results
OUTPUT_SIM_DIR=~/workspace/rpnav/rpnav/simulate/outputs

# Name of simulation script
# Tells the script to look for base simulation configuration of the form ${SIM_SCRIPT}_sim.input
SIM_SCRIPT=parkes

export PSR1="J0835-4510"
export PSR2="J1939+2134"

# Set total observation time in days
obsTime=(0.168 0.33 0.5 1 2 3 5 8 13 21 25 28 31)
#***************************************************

cd $OUTPUT_SIM_DIR

# Collect positioning results for set of integration times
for t in $obsTime
do
    simName=${SIM_SCRIPT}_${t}d
    simInput=${INPUT_SIM_DIR}/${simName}.input
    simDir=${OUTPUT_SIM_DIR}/${simName}_sim
    simResDir=${simDir}/output/real_0
    runscript=${simDir}/scripts/runScripts_master

    # Copy existing sim config
    cp ${INPUT_SIM_DIR}/${SIM_SCRIPT}.input ${simInput}

    # update sim name    
    sed -i "0,/name:.*/s/name:.*/name: "${simName}"_sim/" ${simInput}

    # Get start time from file and update end time
    start=$(sed -n 's/^start://p' ${simInput})
    start=$(bc <<< $start)
    end=$(bc <<< "${start} + ${t}")
    sed -i "s/finish:.*/finish: "${end}"/" ${simInput}

    # Get name of observatory
    # Assumes input file takes the format
    # <obsRun>
    # name: PARKES
    # tel: ...
    observer=$(awk '// {if (lastLine == "<obsRun>"){print $NF}lastLine = $0}' ${simInput})

    for ((i = 0; i < 50; i++));
    do
        rm -rf $simResDir

        # Rerun simulation  
        /home/naverill/external/ptasimulate/ptaSimulate ${simInput}

        # Insert aliasing tempo2 at the start of the script
        echo "alias tempo2 /usr/local/tempo2/bin/tempo2
    setenv TEMPO2 /usr/local/tempo2
    " | \
            cat - ${runscript} >> temp \
            && mv temp ${runscript}

        tcsh ${runscript}

        # Transform simulated par files into TDB
        /usr/local/tempo2/bin/tempo2 -gr transform ${simResDir}/${PSR1}.par ${simResDir}/${PSR1}.tdb.par tdb
        /usr/local/tempo2/bin/tempo2 -gr transform ${simResDir}/${PSR2}.par ${simResDir}/${PSR2}.tdb.par tdb

        # Create output directory
        outResDir=$(printf "${simDir}/output/S%02d" $i)
        mkdir -p ${outResDir}
        simFile=$(printf "results_S%02d.csv" $i)

        # Create results file
        echo "Iteration,Longitude(deg),Latitude(deg),X(m),Y(m),Z(m),Error,Step Size" > ${outResDir}/${simFile}
        for ((j = 0; j < 100; j++));
        do
            echo ""
            echo $simName
            echo ""
            resFile=$(printf "results_S%02dI%02d.csv"  $i $j)

            # Run gradient descent algorithm
            /usr/local/tempo2/bin/tempo2 -gr pulsar_positioning \
                -f ${simResDir}/${PSR1}.tdb.par ${simResDir}/${PSR1}.tim \
                -f ${simResDir}/${PSR2}.tdb.par ${simResDir}/${PSR2}.tim \
                -a grde -e chi -n -r -s ${outResDir}/${resFile} -o $observer
                
            # Add final position and error to results file
            tail -n 1 ${outResDir}/${resFile} >> ${outResDir}/${simFile}
        done
    done
done
