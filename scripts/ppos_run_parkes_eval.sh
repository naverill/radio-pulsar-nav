#!/bin/bash

#************ CONFIGURATION PARAMETERS *************
# Location of the directory containing the .input file
INPUT_DIR=~/workspace/rpnav/rpnav/simulate/inputs/ppta_data/averaged

# Location of the directory to save the results
OUTPUT_DIR=~/workspace/rpnav/rpnav/simulate/outputs
RES_DIR="${OUTPUT_DIR}/parkes_real"

# Set total observation time in days

# psrList=("J1017-7156" "J1024-0719" "J1600-3053" "J1732-5049" "J1909-3744" "J2129-5721" "J2241-5236" "J0613-0200" "J0711-6830" "J1022+1001" "J1045-4509" "J1125-6014" "J1446-4701" "J1545-4550" "J1603-7202" "J1643-1224" "J1713+0747" "J1730-2304" "J1744-1134" "J1824-2452A" "J1832-0836" "J1857+0943" "J1939+2134" "J2124-3358" "J2145-0750") 
# psrList=("J1024-0719" "J1600-3053" "J1732-5049" "J1909-3744" "J2129-5721" "J2241-5236" "J0613-0200" "J0711-6830" "J1022+1001" "J1045-4509" "J1125-6014" "J1446-4701" "J1545-4550" "J1603-7202" "J1643-1224" "J1713+0747" "J1730-2304" "J1744-1134" "J1824-2452A" "J1832-0836" "J1857+0943" "J1939+2134" "J2124-3358" "J2145-0750") 
psrList=("J0613-0200" "J0711-6830" "J1017-7156" "J1022+1001"  "J1024-0719") 
obsTime=(24)
# obsTime=(744 672 600 504 312 192 124 72)

psrNum=${#psrList[@]}
nIter=50
nSim=10
#***************************************************


cd $OUTPUT_DIR

truncDir=${INPUT_DIR}/trunc
mkdir $truncDir
# Collect positioning results for set of integration times
len=$((psrNum-1))

for t in "${obsTime[@]}"
do
    for i in $(seq  0 $len); 
    do
        PSR1="${psrList[i]}"

        for j in $(seq  0 $len); 
        do
            PSR2="${psrList[j]}"
            # Create output directory
            simFile=$(printf "results_${PSR1}_${PSR2}_${t}h.csv")
            outResDir=$(printf "${RES_DIR}/${PSR1}_${PSR2}_${t}h/S${%02d}" $i)
            mkdir -r ${outResDir}

            # Avoid repeating same pulsar or pair
            if [ $i -ge $j ]; then
                continue
            fi

            PSR1File="${INPUT_DIR}/${PSR1}.avg.tim"
            PSR2File="${INPUT_DIR}/${PSR2}.avg.tim"
            PSR1TruncFile="${truncDir}/${PSR1}.avg.tim"
            PSR2TruncFile="${truncDir}/${PSR2}.avg.tim"
            PSR1Len=$(awk 'END { print NR }' $PSR1File)
            PSR2Len=$(awk 'END { print NR }' $PSR2File)
                
            for j in $(seq  0 $nSim); 
            do
                PSR1Start=$(($RANDOM % $PSR1Len))
                PSR1End=$(($PSR1Start + $t))
                while [ $PSR1Start -le 2 ] || [ $PSR1End -ge $PSR1Len ]; do
                    echo "INVALID PSR1"
                    PSR1Start=$(($RANDOM % $PSR1Len))
                    PSR1End=$(($PSR1Start + $t))
                done

                PSR2Start=$(($RANDOM % $PSR2Len))
                PSR2End=$(($PSR2Start + $t))

                while [ $PSR2Start -le 2 ] || [ $PSR2End -ge $PSR2Len ]; do
                    echo "INVALID PSR2"
                    PSR2Start=$(($RANDOM % $PSR2Len))
                    PSR2End=$(($PSR2Start + $t))
                done


                psr1Fend=$((PSR1End + 1))
                psr2Fend=$((PSR2End + 1))
                head -n 2 ${PSR1File} > $PSR1TruncFile
                sed -n "${PSR1Start},${PSR1End}p;${psr1Fend}q" $PSR1File >> $PSR1TruncFile

                head -n 2 ${PSR2File} > $PSR2TruncFile
                sed -n "${PSR2Start},${PSR1End}p;${psr2Fend}q" $PSR2File >> $PSR2TruncFile

                cp $PSR1TruncFile ${outResDir}
                cp $PSR2TruncFile ${outResDir}

                # Create results file
                echo "Iteration,Longitude(deg),Latitude(deg),X(m),Y(m),Z(m),Error,Step Size" > ${outResDir}/${simFile}
                for ((t = 0; t < $nIter; t++));
                do
                    resFile=$(printf "results_${PSR1}_${PSR2}_${t}h_I%02d.csv" $t)
                    echo ""
                    echo $resFile
                    echo ""

                    # Run gradient descent algorithm
                    /usr/local/tempo2/bin/tempo2 -gr pulsar_positioning \
                        -o "PARKES" \
                        -f ${INPUT_DIR}/${PSR1}.par ${PSR1TruncFile} \
                        -f ${INPUT_DIR}/${PSR2}.par ${PSR2TruncFile} \
                        -a grde -e chi -n -r -s ${outResDir}/${resFile} 
                        
                    # Add final position and error to results file
                    tail -n 1 ${outResDir}/${resFile} >> ${outResDir}/${simFile}
                done
            done
        done
    done
done