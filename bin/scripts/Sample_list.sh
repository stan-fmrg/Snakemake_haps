#!/bin/bash
#Sample_list.sh




#Pass the relative path of the raw data from build_hap_config.py script to Sample_list.sh
echo "Relative path of raw Data: "$1

#The raw Data Folder as variable 
rawDataFolder=$1

#To extract all folder names with samples into an array and use them as input files for the mapping and variant calling
sampleFolders=`ls -1 $rawDataFolder | awk '{ print $1 }'`



###################################

#These values can be changed depending on the sequencing technology and whether different adapters must be used.

adapters_paired=TruSeq3-PE-2.fa #Available adapters can be found in data/adapters
adapters_single=TruSeq3-SE.fa #Available adapters can be found in data/adapters

PHRED=phred33 #PHRED can be either phred64 or phred33, however, phred64 is for older technologies and more recent Illumina sequencers use phred33

####################################

workingDir=`pwd`
echo "current working dir is $workingDir"


#make a .csv file for the paired-end samples
cd data/config/
pe_samples=pe_samples.csv
#check if this file exists:
if [ -f $pe_samples ];
then
	#if it does, delete it
   rm $pe_samples
fi
#then create new samples file
touch $pe_samples
cd $workingDir



#make a .csv file for the single-end samples
cd data/config/
se_samples=se_samples.csv
#check if this file exists:
if [ -f $se_samples ];
then
	#if it does, delete it
   rm $se_samples
fi
#then create new samples file
touch $se_samples
cd $workingDir


#Then iterate over these
#iterate over these

for sampleName in $sampleFolders
do
	echo -e "\n++++++++processing sample folder $sampleName"

	#go into the sample folder so we can list files there
	cd $rawDataFolder/$sampleName
	echo "in folder `pwd`"
	if [ $(find . \( -name '*_R2_*' -or -name '*_2_*' -or -name '*_R_*' -or -name '*.R.*' -or -name '*.2.*' -or -name '*.R2.*' \)) ]
		then
		SampleR1=$(find \( -name '*_R1_*' -or -name '*_1_*' -or -name '*_F_*' -or -name '*.F.*' -or -name '*.1.*' -or -name '*.R1.*' \))
		SampleR2=$(find \( -name '*_R2_*' -or -name '*_2_*' -or -name '*_R_*' -or -name '*.R.*' -or -name '*.2.*' -or -name '*.R2.*' \))
		Sample=$(basename -- "$SampleR1")
		extension="${Sample##.*.*}"
		Sample="${Sample%.*.*}"
		echo -ne "$Sample,$sampleName/$SampleR1,$sampleName/$SampleR2,$adapters_paired,$PHRED\n" >> $workingDir/data/config/$pe_samples
		cd $workingDir
	else
		SampleR1=$(find \( -name '*_R1_*' -or -name '*_1_*' -or -name '*_F_*' -or -name '*.F.*' -or -name '*.1.*' -or -name '*.R1.*' \))
		Sample=$(basename -- "$SampleR1")
		extension="${Sample##.*.*}"
		Sample="${Sample%.*.*}"
		echo -ne "$Sample,$sampleName/$SampleR1,$adapters_single,$PHRED\n" >> $workingDir/data/config/$se_samples
		cd $workingDir
	fi
done

sed -i '1i sample,forward,reverse,adapter,phred' data/config/$pe_samples
sed -i '1i sample,forward,adapter,phred' data/config/$se_samples
