#!/bin/bash


labels=results/Final_results/Microhaplot/labels.txt
SAM_files=results/map/SAM/

#extract all the folder names to an array
sampleFiles=`ls -1 $SAM_files| awk '{ print $1 }'`

echo $sampleFiles

for sampleName in $sampleFiles
do
	Sample=$(basename -- "$sampleName")
	extension="${Sample##.*.*}"
	Sample="${Sample%.*.*}"

	echo -ne "$sampleName\t$Sample\t"NA"\n" >> $labels
done

