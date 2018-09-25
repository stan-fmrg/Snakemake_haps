#!/bin/bash
FASTA=$1
samtools faidx $1
awk '{ print $1 "\t0\t" $2}' $1.fai
