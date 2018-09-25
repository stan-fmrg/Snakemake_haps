#!/bin/bash
#multiQC.sh


input=$(realpath $1)
output=$(realpath $2)


multiqc -f -m fastqc $input -o $output