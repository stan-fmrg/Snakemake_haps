# A snakemake workflow for high-throughput amplicon sequencing


This is a snakemake workflow for generating a list of variants (Single Nucleotide Polymorphisms- SNPs) from amplicon sequencing data and the production of haplotypes using these variants with rad_haplotyper, a program designed to produce SNP haplotypes and is included in this workflow.  
The workflow is designed to handle both single-end and paired-end sequencing data, and it gives the user the flexibility to filter the identified SNPs under stringend conditions and use these for haplotype calling, or use the overlapping SNPs with an existing list of high-quality variants and use only these for haplotype calling.
Additionally, the workflow has the option for the pipeline to produce only the variants without having the haplotypes called, to reduce the analysis time when only the variants are needed.


## Overview

The standard snakemake_haps workflow essentially performs the following steps: 

##### 1) Quality control of the raw reads.
##### 2) Quality trimming and adapter removal of raw reads.
##### 3) Mapping of the processed reads against the reference sequence.
##### 4) Duplicate reads based on the resulting alignments are marked. 
##### 5) The resulting alignments are sorted and indexed.
##### 6) Variants are called using FreeBayes
##### 7) Variants are filtered with very stringent parameters to identify high quality SNPs
##### 8) Haplotypes are called using rad_haplotyper


## Getting started

If you are new to [Snakemake], you might want to start by working the available [tutorial for beginners] with instructions for how to install Snakemake, as well as Miniconda Python 3 and setting an environment with all required software.

## 1. Installation

"Snakemake_haps" can be downloaded and installed from a Github repository ("link"), however, since "snakemake_haps" has many dependencies, it is suggested to use one of the [Anaconda] or [Miniconda] repositories to avoid compiling dependecies one at a time which is time consuming.

### 1.1 Downloading and Installing (Ana|mini)conda:
- [Anaconda] - Free
- [Miniconda] - Free

#### For installing Anaconda:
##### - Steps for downloading and installing Anaconda:

```bash
# Curl is needed so if not installed you can install it by typing:
sudo apt install curl	

# Installing Anaconda
curl -o /tmp/Anaconda.sh https://repo.continuum.io/archive/Anaconda3-5.1.0-Linux-x86_64.sh && bash /tmp/Anaconda.sh
```

##### or

##### Download the available installer version of Anaconda from [Anaconda]

```bash
#Change directory to where the installer is
cd ~/Path/to/Anaconda.sh
#Run the installer by typing 
bash Anaconda.sh
```

#### For installing Miniconda Pyhon 3 distribution:
##### - Steps for downloading and installing Miniconda3:

```bash
# If Curl is not installed you can install it by typing:
$sudo apt install curl	

# Installing Miniconda3t
$curl -o /tmp/miniconda.sh https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash /tmp/miniconda.sh
```

##### or

##### Download the available installer version of Miniconda3 from [Miniconda]

```bash
#Change directory to where the installer is
$cd ~/Path/to/Miniconda3.sh
#Run the installer by typing 
$bash Miniconda3-latest-Linux-x86_64.sh
```

###### - Answer yes to the question whether conda shall be put into your PATH.

### 1.2 Installing Snakemake workflow manager:

##### Via Bioconda channel:

```bash
#Installing Snakemake from the Bioconda channel:
$conda install -c bioconda snakemake
```
##### With a working Python >=3.5 setup, installation of Snakemake can be performed via PyPi:

```bash
#Installing Snakemake via PyPi:
$pip3 install snakemake
```

##### or

```bash
easy_install3 snakemake
```

### 1.3 Downloading and installing "Snakemake_haps" files:

#### Description

This repository contains scripts that are part of a complete workflow for DNA-seq analysis and mapping reads to a reference genome for the identification of SNPs and designing of haplotypes. 

#### Clone the repository using the Command Line:

```bash
#Download files
git clone https://github.com/chollenbeck/snakemake_haps.git

#Change directory to snakemake_haps
cd snakemake_haps
```

##### or

#### Visit ("link") and select "Clone or download" option and select "Download Zip"
#### unzip the zip. file to the preffered location and then change directory:

```bash
#Change directory to snakemake_haps
cd snakemake_haps
```

### 1.4 Create an environment

#### For using "Snakemake_haps" it is suggested to create an environment dedicated to "Snakemake_haps". This environemnt is a directory that contains all dependencies for "Snakemake_haps" and it can be activated before running the workflow and deactivatedd after it finishes. The environment configuration for "Snakemake_haps" is located in /bin/install/environment.yml and can be created as follows:

```bash
#Installing available environment
conda env create --name *environment_name* --file /bin/install/environment.yml

#Activate environment
source activate *environment_name* 
```

#### and deactivate it by using 
```bash
#Deactivating environment
source deactivate
```

### 1.5 Quick run of "Snakemake_haps"

For running "Snakemake_haps with default parameters, the following set of commands must be executed.

1) Change home directory where "Snakemake_haps" is downloaded:

```bash
#Change directory
cd /home/Snakemake_haps_directory/
```

2) Run build_hap_config.py script to create the config.yaml file with all analysis parameters and paths of the raw fastq files. 

```
#Run the build_hap_config.py script
python bin/scripts/build_hap_config.py -g data/genome/"name_of_ref_seq.fasta" -d /home/Path/of/Fastq/files
```

3) Run "Snakemake_haps"

```
#Start workflow by running snakemake
snakemake
``` 

_A more detailed explanation of the above steps are given in following "Tutorial"._

### 2 Tutorial

This next part, will use the available sample data files in [data/fastq_raw] to illustrate how to run the workflow using both paired-end and single-end reads. The single-end data comes from a sequencing (using MiSeq technology) of two King scallop samples (*Pecten maximus*), while the paired-end reads are coming from Atlantic salmon (*Salmo salar*) sequencing. For testing purposes, the sample files can be found in the "Snakemake_haps" data/fastq_raw directory/Single_end_reads and data/fastq_raw directory/Paired_end_reads directory. Their corresponding reference sequences can be found 

#### 2.1 Steps for running "Snakemake_haps" with default parameters

- After activating the newly created environment, the working directory must be set where the "Snakemake_haps" workflow is .

```bash
#Change directory
cd /home/Snakemake_haps_directory/
```

##### 2.1.1 To use the default workflow parameters:

- Raw data must be FASTQ sequencing files and stored in the data/fastq_raw/ directory in the "Snakemake_haps" folder.
- A reference sequence must also be provided and stored in the data/genome/ directory in the "Snakemake_haps" folder.

##### 2.1.2 Create the config.yaml file with default parameters, which consists the list of the samples and their paths in /data/fastq_raw/ and the reference sequence path in data/genome/, as well as mapping and variant filtering parameters.

```
#Run the build_hap_config.py script
python bin/scripts/build_hap_config.py -g /data/genome/"name_of_ref_seq.fasta"
```

This python script will create a .yaml file consisting all the available parameters needed for the workflow to succesfully run, as well as all the paths of the FASTQ sample files with their designation as paired-end or single-end reads.

**- The config.yaml file is then stored in the home directory of "Snakemake_haps" and it can be opened by using a typical text editor. It is advisable for the user, before running the workflow, to open the config file and have a visual inspection of the
sequencing read path, as well as the default parameters of mapping, variant calling and variant filtering.** 
**- Other options that can be controled and changed from the config.yaml file are the option of including the QC step in the workflow and the removing the haplotype calling step.**

##### 2.1.3 After creating config.yaml file with all the necessary parameters and the correct paths of the raw sequencing data, the workflow can be executed. 

```
#Start workflow by running snakemake
snakemake
``` 

#### 2.2 Steps for running "Snakemake_haps" only for variant calling. 

The "Snakemake_haps" program is coming with the option to be used only for variant calling when the designing of haplotypes are not needed. To configure "Snakemake_haps" to proceed with the analysis without the [rad_haplotyper] step, [build_hap_config.py] script (described in 2.1.2) must be used with a changed option ##-q## as False (Default value = True). 

```
#Run the build_hap_config.py script
python bin/scripts/build_hap_config.py -g /data/genome/"name_of_ref_seq.fasta -q False"
```
That creates a [config.yaml] file with the [rad_haplotyper] de-activated with all other parameters remaining the same.

Then "Snakemake_haps" workflow can be activated as desctibed in 2.1.3

```
#Start workflow by running snakemake
snakemake
``` 

#### 2.3 Steps for running "Snakemake_haps" without QC step **(Not suggested)**. 

In the same way as before, "Snakemake_haps" can be used for variant analysis and haplotype calling, however it is possitble to skip the QC step of the sequencing reads, which can take a huge amount of time. To do that, [build_hap_config.py] script has to be run by using the option **-p** as False (Default value = False).

```
#Run the build_hap_config.py script
python bin/scripts/build_hap_config.py -g /data/genome/"name_of_ref_seq.fasta -p False"
```
That creates a [config.yaml] file with the [rad_haplotyper] de-activated with all other parameters remaining the same.

Then "Snakemake_haps" workflow can be activated as desctibed in 2.1.3

```
#Start workflow by running snakemake
snakemake
``` 

#### 2.4 Steps for running "Snakemake_haps" without both QC and haplotype calling:

For using "Snakemake_haps" only for variant calling with no QC and haplotype calling, the two previiously explained options  **-p** and **-q** can be combined when [build_hap_config.py] is run.

```
#Run the build_hap_config.py script
python bin/scripts/build_hap_config.py -g /data/genome/"name_of_ref_seq.fasta -p False -q False"
```

Then "Snakemake_haps" workflow can be activated as desctibed in 2.1.3

```
#Start workflow by running snakemake
snakemake
``` 

#### 2.4 Other parameters that can be changed in [build_hap_config.py]

As mentioned before, [build_hap_config.py] is the script that creates the [config.yaml] that contains all the necessary parameters for "Snakemake_haps" to run, and there is a number of other options that can be configured. some of these are   

- **"MANDATORY"** **"-g", "--genome"**: This option must be included by the user when [build_hap_config.py] is run, since the complete relative path and name of the reference sequence is needed for the workflow. If the sequence is not in data/genome directory, the correct relative path can also be added by using the option **-g** followed by the name of the .fasta file (Default value: data/genome/).
- **"-n", "--run_name"**: When there is the need to use a different name for each run (Default value: haps_run).
- **"-a", "--adapter_dir"**: If adapters for adapter removal are not located in the directory data/adapters, the new relative path for the adapters can be added by using the option **-a** (Default value: data/adapters).
- **"-s", "--split_genome"**: This option splits the reference genome into n pieces for parallel variant calling (Default value: 5)
- **"-v", "--diff_vcf"**: This option allows the "Snakemake_haps" workflow use a VCF file with variants that the user has previsouly identified as of high-quality and only the overlapping SNPs between this VCF and the identified SNPs in the VCF file of the "Snakemake_haps" workflow will be used for the haplotype calling. The relative path of the VCF file can be added to the workflow with option **-v** (Default value: "Null").

#### 2.5 Parameters that can be changed in [config.yaml]

As mentioned before, [config.yaml] is the created file from running [build_hap_config.py] script and contains all of the parameters for each step of the workflow, including mapping, variant calling, and variant filtering. These parameters, even though they work very good with sequencing data coming from Miseq runs, they can be non-ideal for other sequencing data or for the needs of each project, and be edited by opening [config.yaml] by using a text editor.  

- **"Trimmomatic read trimming"**: Trimmomatic performs a variety of useful trimming steps for illumina paired-end and single ended reads. The trimming tasks that it performs include the removal of adapters, leading and trailing low quality bases, the scanning of the read with a 4-base wide sliding window and cutting when the average quality per base dropw below a threshold, and the removal of low quality reads. The parameters for these tasks, that can be also found in the [config.yaml], are the following:   
	- **LEADING:[20]**  Cut bases off the start of a read, if below a threshold quality
	- **TRAILING:[20]** Cut bases off the end of a read, if below a threshold quality
	- **AVGQUAL:[30]** Drop the read if the average quality is below the specified level
	- **MINLEN:[100]** Drop the read if it is below a specified length
	- **SLIDINGWINDOW:[4:15]**  Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below a threshold.

The values of these parameters can be easily modified by changing the values inside te brackets and then save the changes. For more information about Trimmomatic and the available options please go to [Trimmomatic].

- **"BWA-MEM"**: BWA-MEM algorithm performs the local alignment of the reads to a reference sequence. The parameters taken for the mapping of the reads in this workflow are the default, with the only the number of threads used being changed from default for better performance in multithreading mode. However, there are more options available in [BWA-MEM] that can be used and added in pconfig.yaml[] file.
	- **"-t" INT**: Number of threads [8]
	- **"-B" INT**: Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4] (Default)
	- **"-O" INT**: Gap open penalty. [6] (Default)

The values of these parameters can be easily modified by changing the values inside te brackets and then save the changes. For more information about BWA-MEM algorithm and the available options please go to [BWA-MEM].

- **"Samtools Fixmate"**: This step fill in mate coordinates, ISIZE and mate related flags from a name-sorted alignment. This step is essential for the duplicate read marking and removal with samtools markdup.
	- **"-r"**: Remove secondary and unmapped reads.
	- **"-m"**: Add ms (mate score) tags. These are used by markdup to select the best reads to keep.
	
 - **"Samtools MarkDup"**: Mark duplicate alignments from a coordinate sorted file that has been run through fixmate with the -m option. This program relies on the MC and ms tags that fixmate provides.
	- **"-s"**: Print some basic stats.
	Extra options:
	- **"-r"**: Remove duplicate reads. (Default: Null)
	- **"-S"**: Mark supplementary reads of duplicates as duplicates. (Default: Null)
	
To activate these options, just add **-r** and/or **-S** next to **-s** in [config.yaml] file and then save it. 

- **"Samtools View"**: Samtools view command is a very versatile command from the Samtools package where its main function is to convert the binary form of the alignments into a human readable format. However, Samtools view has other functions, such as count the total lines of a BAM file and exctract reads that fulfill specific characteristics, like paired-reads mapped properly, unmapped reads and mates that are unmaped. For more information about the SAM flags that correspond to ecah property, visit [SAMformat] webpage.
For this workflow, the parameters taken for filtering the reads that were aligned are:	
	- **"-f" _INT_**: Keep only reads that have the _INT_ present in the FLAG field [3] (Default value: Keep paired reads  and mapped in proper way (for paired-end reads only))
	- **"-F" _INT_** : Discard reads that have the _INT_ present in the FLAG field [780] (Default value: Discard reads that are unmapped, the mate is unmmaped, the read has split hits, and read fails quality checks.
	
The values of these parameters can be easily modified by removing them or adding new ones and then save the changes. For more information about Samtools fixmate, markdup, and view and their available options please go to [Samtools].

- **"FreeBayes Variant Calling"**: FreeBayes is a genetic variant detector, specifically designed to identify small polymorphisms, IndDels, MNPs and compelx events (composite insertion and substitution events). It uses short-read alignments (BAM files), which in this workflow all BAM files for a number of samples are merged into one BAM file, and a reference genome (in FASTA format) are used to determine the most-likely combination of genotypes for the population at each position in the reference:
 	- **--haplotype-length**: [0] Simply annotate observation counts of SNPs and indels:
	- **--no-complex**: Do not generate information about complex events
	- **--use-duplicate-reads**: Duplicate reads are included in the analysis
	
The values of these parameters can be easily modified by removing them or adding new ones and then save the changes. For more information about FreeBayes please go to [FreeBayes].

- **"VCFtools Variant Filtering"**: VCFtools is a program package and a toolset that is designed to work wth VCF file and perform a number of different operations on VCF files. On this wokrflow, VCFtools is used to filter out specific variants based on their available tags. These filters are:
	- **--minDP > 10x**: Includes only genotypes greater than or equal to the "--minDP" value.
	- **--minGQ: > 20**: Exclude all genotypes with a quality below the threshold specified. 
	- **--Max-missing: 0.25**: Exclude sites on the basis of the proportion of missing data (defined to be between 0 and 1, where 0 allows sites that are completely missing and 1 indicates no missing data allowed).
	- **--MAF > 0.05**: Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value

The values of these parameters can be easily modified by removing them or adding new ones and then save the changes. For more information about VCFtools please go to [VCFtools].

- **"Rad_haplotyper"**: The filtered set of SNPs from previous step are used for haplotype calling using rad_haplotyper tool (Willis et al., 2017), which is a program designed to produce SNP haplotypes. Haplotyping in this program is based on the idea that single- and paired-end reads capture the phase of linked SNPs. Apart from the filtered VCF with high quality SNPs, a BAM file with the aligned reads is also required. 
This program works first by iterating through each locus attempting to construct haplotypes from the raw reads and then by comparing the haplotypes to the bases appeared in the VCF file, haplotypes that are unique only on the raw reads are removed while haplotypes that can be matched with the bases called in the VCF are kept. The filter options that are taken for this step are:
	- **--miss_cutoff: 0.75**": Proportion of missing data cutoff for removing loci from the fnal output.
	- **--hap_rescue: 3**": Specify a rescue parameter that controls the behavior of the cript when dealing with loci that have more observed haplotypes han are possible given the genotypes.
	
	The values of these parameters can be easily modified by removing them or adding new ones and then save the changes. For more information about rad_haplotyper please go to [rad_haplotyper].


### 3 Results & produced files

The results of the analysis are stored in multiple locations, based on which step of the workflow they were produced (QC, mapping, variant calling, haplotype calling). However, the main results of the workflow are stored in the results/Final_results/ directory.
The files that are produced during the workflow contain figures and statistis about the reads and their alignemnt against the reference sequence, as well as two VCF files with the raw and filtered variants from the variant analysis and the haplotypes that are designed using [rad_haplotyper]. A detailed description of these files is going to be given over the next part. 

#### 3.1 - Quality_control: The quality control directory contains the statistics about the quality of the reads before and after the quality trimming and adapter removal step of the workflow by using FastQC. The FASTQ files received from the sequencing can have various artefacts, which prior to variant calling must be removed. These artefacts can either be low quality bases at the end of the reads and contamination from Illumina adapters.  

##### The main functions of FastQC are:
###### -Import of data from FastQ files
###### -Providing a quick overview to tell you in which areas there may be problems
###### -Summary graphs and tables to quickly assess your data
###### -Export of results to an HTML based permanent report
###### -Rather than looking at quality scores for each individual read, FastQC looks at quality collectively across all reads within a sample. 

##### All figures and statistics are then stored in results/final_results/Quality_control/QC_reads and results/final_results/Quality_control/Raw_reads directories, depending on whether the FastQC was performed before or after the quality trimming step.

**As mentioned earlier, the FastQC step can be avoided (even though it is not suggested), thus, the final FastQC files will be missing at the end of the workflow.** 

#### 3.2 BAM_files: This directory contains the merged sequence alignment data from the mapping of all samples against the refernece sequence. The BAM format is the compressed binary version of the SAM file which is produced from the alignent of each sample against the reference sequece, and the conversion from SAM to BAM can save valuable space when storing a large number of samples.
- This merged BAM file is the main input file for the variant analysis with FreeBayes.

#### 3.3 Variants: After completing the alignment step where all the processed reads are mapped against the reference sequence, Freebayes proceeds with the identification of variants. The identified calls are then stored in a VCF format file where it can be further filtered to remove low quality variants and false positive SNPs. The results from this step are stored in two different directories in results/Final_results/Variants/Raw_Variants or results/Final_results/Variants/Filtered_Variants, based on whether they are the raw variants or after filtering.
##### - Inside the Filtered_Variants directory, a .txt format file with p-values for each site from a Hardy-Weinberg Equilibrium test are reported, as well as the Observed numbers of Homozygotes and Heterozygotes and the corresponding Expected numbers under HWE.   
##### - In case a set of pre-defined high-quality SNPs have been used for the selection of overlapping SNPs with the newly identified SNPs, a new VCF file will be stored in the results/Final_results/Variants/Filtered_Variants directory, where only the common SNPs between the two datasets will be present. This file can then be used for the haplotype calling.

#### 3.4 Haplotypes: The results from using the previously indentified SNPs to create haplotypes are stored in results/Final_results/Haplotypes. Various files with information about haplotypes can be found in this directory such as:

###### - codes.haps.gen: This file contains a complete list of haplotypes that passed the quality control of [rad_haplotyper] for each locus.
###### - haps.gen: This file has the list of the loci that have haplotypes that passed the quality control of [rad_haplotyper].
###### - stats.out: This file lists the number of sites at the locus (SNPs or indels), the number of haplotype alleles observed in the dataset, the number of individuals haplotyped, the total number of individuals in the data set, and the proportion of individuals successfully haplotyped per locus. It also indicates the status of each locus as either passed/failed, and if it failed, the possible reason for failure. Possibilities include: possible paralogs, possible low-coverage/genotyping errors, previously missing genotypes, or a complex polymorphism that is difficult to haplotype. See the section on haplotype calling for more details about how this works.
**This (tidy) file can be loaded into R and can be used to further filter loci from the data set based on user-defined cut-off values.**
###### - ind_stats.out: This file lists number of loci that are possible paralogs, have low coverage/errors, missing genotypes, number of failed loci, the total number of loci, and the proportion of successful loci on a per individual basis. This file can be used to identify and remove problematic individuals from the final data set.

#### 3.5 Documents: This directory contains all the files that contain statistics about the reads and the QC step, as well as the mapping of the reads. Some of the information that can be found in these files is:

###### - The total number of reads per sample,
###### - The number of mapped and unmapped reads per sample,
###### - The number of duplicate reads,
###### - The percentage of properly paired if the reads are paired-end,
###### - The minimum and maximum length of the reads,
###### - Their average quality,
and
###### - The total number of reads aligned to every locus

 
Snakefiles **(These files should not be changed)**

In snakemake, snakefiles are the files that contain the rules that denote how to create output files from input files. For this workflow, snakefiles are categorised based on their function and on which step of the workflow they are activated.

-[raw.py] creates links with the reference genome and the fastq files for downstream analysis by the workflow.

-[folders.py] creates folders for all output files produced by each step of the workflow.

-[qc.py] This snakefile includes the adapter removal by Trimmomatic which removes Illumina adapters that can cause a problem downstream analysis and to the identification of high quality SNPs. Additionally, Trimmomatic trims low quality regions of reads. Apart from Trimmomatic, qc.py contains a quality control step using fastQC tool to confirm the high quality of trimmed reads producing figures and tables with statistics which are explained in [fastQC] help page.

-[map.py] This script contains the full workflow for mapping the reads to a reference genome by using [BWA-MEM] mapper (Li, 2009) and [FreeBayes] (Garrison et al., 2012) variant caller. 
- Steps
	1) Indexing of the reference genome.
	2) Mapping the trimmed reads against the reference genome and sorting the output [BAM] file.
	3) Apply [SAMtools] -fixmate and -markdup utilities to sorted BAM files for the removal of potential PCR duplicates. The output file of this step must be indexed.
	4) For variant calling, the output BAM files from mapping each sample to the reference genome must be merged and indexed.

-[variant_calling.py] contains the workflow for the identification of variants from the produced BAM files. Initially this script splits the reference sequence into  multiple components for faster variant analysis. FreeBayes is the variant caller which requires as input files the indexed merged BAM file and produes a VCF file with raw variants. A filtering procedure is included in this script to remove low quality variants from the final variant set.

-[haplotype.py] builds SNP haplotypes using the filtered variants from previous steps, providing a test for paralogy.

You can read more about the method in the available manual of [rad_haplotyper] and the following publication:

Willis, S. C., Hollenbeck, C. M., Puritz, J. B., Gold, J. R. and Portnoy, D. S. (2017), Haplotyping RAD loci: an efficient method to filter paralogs and account for physical linkage. Mol Ecol Resour, 17: 955–965. doi:10.1111/1755-0998.12647 [link]

-[clean.py] removes all uneccessary files and folders created thoughout the analysis


## Data

This repository includes 2 directories with fastq files of 3 samples each, one of which is single-end fastq files and the other is paired-end fastq files.These files can be used as exampels with the tutorial of this manual. Additionally, two reference sequences for the two types of raw reads are available, one for the single-end and one for the paired-end sample tutorial.

**Single-end reads**
[data/fastq_raw/Single_end_sample_fastq/]
- Sample_1
	- `S1_001` 
		- `S1_S1_L001_R1_001.fastq.gz`
- Sample_2
	- `S2_001` 
		- `S2_S2_L001_R1_001.fastq.gz`
- Sample_3
	- `S3_001` 
		- `S3_S3_L001_R1_001.fastq.gz`
			
		
**Paired-end reads**
[data/fastq_raw/Single_end_sample_fastq/]

- Sample_1
	- `S1_001` 
		- `FR14992874_S1_L001_R1_001.fastq.gz`: Containing the forward reads of Sample 1
		- `FR14992874_S1_L001_R2_001.fastq.gz`: Containing the reverse reads of Sample 1
		
- Sample_2
	- `S2_001` 
		- `FR14992882_S2_L001_R1_001.fastq.gz`: Containing the forward reads of Sample 2
		- `FR14992882_S2_L001_R2_001.fastq.gz`: Containing the reverse reads of Sample 2
		
- Sample_3
	- `S3_001` 
		- `FR14992890_S3_L001_R1_001.fastq.gz`: Containing the forward reads of Sample 3
		- `FR14992890_S3_L001_R2_001.fastq.gz`: Containing the reverse reads of Sample 3
		
		
## Scripts

-[config.yaml] contains links with configuration files and parameters that are used in the workflow.

-[Snakefile] contains the rules that are used in the workflow, written as python script and stored in home directory. This file creates all the links with scripts hat contain the rules of the workflow. Executing [Snakefile] by using the command "snakemake", starts running the workflow specified.

-[build_hap_config.py] creates the [config.yaml] with the default parameters for mapping and filtering, as well as the links to the raw data of the samples.






## References

Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012

Li H. and Durbin R. (2009). Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25(14), 1754-1760.




[Snakfile]: https://github.com/chollenbeck/snakemake_haps/bin/snakefiles/
[config.yaml]: https://github.com/chollenbeck/snakemake_haps/blob/master/config.yaml
[data/fastq_raw]: https://github.com/chollenbeck/snakemake_haps/tree/master/data/fastq_raw
[data/genome]: https://github.com/chollenbeck/snakemake_haps/tree/master/data/genome
[raw.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/raw.py
[folders.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/folders.py
[map.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/map.py
[qc.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/qc.py
[variant_calling.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/variant_calling.py
[haplotype.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/haplotype.py
[clean.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/snakefiles/clean.py
[bin/snakefiles]: https://github.com/chollenbeck/snakemake_haps/tree/master/bin/snakefiles 
[Snakemake]: https://snakemake.readthedocs.io/en/stable/
[tutorial for beginners]: http://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html
[fastqc]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
[BAM]: http://samtools.github.io/hts-specs/SAMv1.pdf
[SAMtools]: http://samtools.sourceforge.net/
[BWA-MEM]: http://bio-bwa.sourceforge.net/bwa.shtml
[FreeBayes]: https://github.com/ekg/freebayes
[rad_haplotyper]: https://github.com/chollenbeck/rad_haplotyper
[build_hap_config.py]: https://github.com/chollenbeck/snakemake_haps/blob/master/bin/scripts/build_hap_config.py
[Trimmomatic]: http://www.usadellab.org/cms/?page=trimmomatic
[SAMformat]: http://www.samformat.info/sam-format-flag





#Snakemake_haps
