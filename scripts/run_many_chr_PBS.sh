#!/bin/bash

#setting variables

while getopts s: option
do
        case "${option}"
        in
                s) SPECIES=${OPTARG};;
        esac
done



CHR_ARRAY=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 scaffolds)

for CHR in ${CHR_ARRAY[*]}
do
	echo "/g/data1/xe2/projects/euc_hybrid_adapt/workspace/scripts/Raijin_angsd_step2_arg.sh -s $SPECIES -c $CHR" | qsub -N $SPECIES.$CHR.saf -P gh6 -l jobfs=1GB,ncpus=1,mem=2G,walltime=2:00:00
done

