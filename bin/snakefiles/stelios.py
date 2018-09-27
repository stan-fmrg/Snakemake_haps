#!/bin/bash

SAM_files=/results/sam/


#extract all the folder names to an array
sampleFiless=`ls -1 $SAM_files| awk '{ print $1 }'`

echo $sampleFiles

