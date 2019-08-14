#!/bin/bash

#Run pairwise analysis on all pairs

indir=.
suffix=.saf.idx
#ls -r $indir/*$suffix>prefix.list

while read line1
do
while read line2
do
	if [ "$line1" = "$line2" ]; then
	   break
	fi
	# prefix1=`basename -s $suffix $line1` # works on EucBox
	# prefix2=`basename -s $suffix $line2`
	prefix1=`basename $line1 $suffix` #syntax for Raijin
	prefix2=`basename $line2 $suffix`
	echo "Raijin_2dsfs_angsd_arg.sh -1 $prefix1 -2 $prefix2 -c Chr01" | qsub -N $prefix1.$prefix2.2dsfs
	# bash Raijin_2dsfs_angsd_test_arguments.sh -1 $prefix1 -2 $prefix2 -c Chr01
done<prefix.list
done<prefix.list

