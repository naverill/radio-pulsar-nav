#!/bin/bash

#************ CONFIGURATION PARAMETERS *************
# Location of the directory containing the .input file
INPUT_SIM_DIR=~/workspace/rpnav/rpnav/simulate/inputs

# Location of the directory to save the results
OUTPUT_SIM_DIR=~/workspace/rpnav/rpnav/simulate/outputs

# Name of simulation script
# Tells the script to look for base simulation configuration of the form ${SIM_SCRIPT}_sim.input
SIM_SCRIPT=woodchester_strong

# Set total observation time in days
# obsTime=(0.168 0.33 0.5 1 2 3 5 8 13 21 25 28 31)
# obsTime=(0.5 1 2 3 5 8 13 21 25 28 31)
obsTime=(21 25 28 31)

nSims=50
nIter=50

psrList=("J0835-4510" "J1017-7156" "J1024-0719" "J1600-3053" "J1732-5049" "J1909-3744" "J2129-5721" "J2241-5236" "J0613-0200" "J0711-6830" "J1022+1001" "J1045-4509" "J1125-6014" "J1446-4701" "J1545-4550" "J1603-7202" "J1643-1224" "J1713+0747" "J1730-2304" "J1744-1134" "J1824-2452A" "J1832-0836" "J1857+0943" "J1939+2134" "J2124-3358" "J2145-0750") 
psrNum=${#psrList[@]}
err="chi"
#***************************************************

cd $OUTPUT_SIM_DIR

# Collect positioning results for set of integration times
for t in "${obsTime[@]}"
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
    

    for ((i =0; i < $nSims; i++));
    do
        rm -rf $simResDir

        psr1N=$(($RANDOM % psrNum))
        psr2N=$(($RANDOM % psrNum))

        # Update pulsar names in sim script
        PSR1="${psrList[$psr1N]}"
        PSR2="${psrList[$psr2N]}"
        echo $PSR1
        echo $PSR2
        sed -i "0,/psr: name=.*/s/psr: name=.*/psr: name="${PSR1}"/" ${simInput}
        sed -i "/psr: name=.*/ {n; :a; /psr: name=.*/! {N; ba;}; s/psr: name=.*/psr: name="${PSR2}"/; :b; n; $! bb}" ${simInput}
        
        sed -i -E "0,/observe: psr=J[0-9]*[+-][0-9]*[A]?/s/observe: psr=J[0-9]*[+-][0-9]*[A]?/observe: psr="${PSR1}"/" ${simInput}
        sed -i -E "/observe: psr=J[0-9]*[+-][0-9]*[A]?/ {n; :a; /observe: psr=J[0-9]*[+-][0-9]*[A]?/! {N; ba;}; s/observe: psr=J[0-9]*[+-][0-9]*[A]?/observe: psr="${PSR2}"/; :b; n; $! bb}" ${simInput}

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
        cp ${simInput} ${outResDir}
        simFile=$(printf "results_S%02d.csv" $i)

        # Create results file
        echo "Iteration,Longitude(deg),Latitude(deg),X(m),Y(m),Z(m),Error,Step Size" > ${outResDir}/${simFile}
        for ((j = 0; j < $nIter; j++));
        do
            resFile=$(printf "results_S%02dI%02d.csv"  $i $j)
            echo ""
            echo $simName
            echo $resFile
            echo $PSR1
            echo $PSR2
            echo ""

            # Run gradient descent algorithm
            /usr/local/tempo2/bin/tempo2 -gr pulsar_positioning \
                -f ${simResDir}/${PSR1}.tdb.par ${simResDir}/${PSR1}.tim \
                -f ${simResDir}/${PSR2}.tdb.par ${simResDir}/${PSR2}.tim \
                -a grde -e ${err} -n -r -s ${outResDir}/${resFile} -o $observer
                
            # Add final position and error to results file
            tail -n 1 ${outResDir}/${resFile} >> ${outResDir}/${simFile}
        done
    done
done
