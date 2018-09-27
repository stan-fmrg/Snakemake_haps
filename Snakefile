
#### DO NOT CHANGE THIS FILE ####


shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

# Imports
import os
import argparse
import sys

os.chmod('bin/scripts/fa_to_bed.sh', 0o755)

os.chmod('bin/scripts/Sample_list.sh', 0o755)

os.chmod('bin/scripts/flash_reads.sh', 0o755)

os.chmod('bin/scripts/Make_labels.sh', 0o755)
###########################################################################################################################################################

var = []


if config["illumina_pe"] is None:
	SAMPLES_SE_RAW = config["illumina_se"] if config["illumina_se"] is not None else []
elif config["illumina_se"] is None:
	SAMPLES_PE_RAW = config["illumina_pe"] if config["illumina_pe"] is not None else []



if config["illumina_se"] is None:
	SAMPLES_PE = [x for x in SAMPLES_PE_RAW]
elif config["illumina_pe"] is None:
	SAMPLES_SE = [x for x in SAMPLES_SE_RAW]



if config["illumina_se"] is None:
	if SAMPLES_PE == []:
		sys.exit("Paired-end reads cannot be found")
	else:
		pass
elif config["illumina_pe"] is None:
    if SAMPLES_SE == []:
        sys.exit("Single-end reads cannot be found")
    else:
        pass


BLOCK_THREADS = 99999
ALL_THREADS = 24

snakefiles = "bin/snakefiles/"

include: snakefiles + "folders.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "map.py"
include: snakefiles + "qc.py"
include: snakefiles + "variant_calling.py"
include: snakefiles + "haplotype.py"
include: snakefiles + "BAM_to_SAM.py"


#########################################################################################################

if config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is False and config["Validation"] is None and config["Flash"] is False:
	rule all:
		input:
		  	expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}"],
		    	sample = SAMPLES_PE, pair = "1 2".split())
		  

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is False and config["Validation"] is None and config["Flash"] is False:
	rule all:
		input:
			expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat"],
				sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_pe"] is None and config["QC_Tests"] is True and config["Haplotyper"] is False and config["Validation"] is None:
	rule all:
		input:
			expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.fastqc", fastQC + "{sample}.se"],
		  		sample = SAMPLES_SE)

elif config["illumina_pe"] is None and config["QC_Tests"] is False and config["Haplotyper"] is False and config["Validation"] is None:
	rule all:
		input:
			expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_DOCS_GENERALSTATS + "{sample}.generalStats"],
				sample = SAMPLES_SE)


#########################################################################################################


elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is True and config["Validation"] is None and config["Flash"] is False:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}"],
		    sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is True and config["Validation"] is None and config["Flash"] is False:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat"],
			sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_pe"] is None and config["QC_Tests"] is True and config["Haplotyper"] is True and config["Validation"] is None:
	rule all:
		input:
			HAPLOTYPE_DIR + "haplotypes.done",
			expand([MICROHAPLOT + "filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.fastqc",fastQC + "{sample}.se"],
		  		sample = SAMPLES_SE)

elif config["illumina_pe"] is None and config["QC_Tests"] is False and config["Haplotyper"] is True and config["Validation"] is None:
	rule all:
		input:
			HAPLOTYPE_DIR + "haplotypes.done",
			expand([MICROHAPLOT + "filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat"],
				sample = SAMPLES_SE)


#########################################################################################################

elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is True and config["Validation"] is not None and config["Flash"] is False:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
		    sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is True and config["Validation"] is not None and config["Flash"] is False:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
		  	sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_pe"] is None and config["QC_Tests"] is True and config["Haplotyper"] is True and config["Validation"] is not None:
	rule all:
		input:
			HAPLOTYPE_DIR + "haplotypes.done",
			expand([MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.fastqc", fastQC + "{sample}.se", FINAL_VARIANTS + "Overlapping_SNPs_se", FINAL_VARIANTS + "Overlapping_SNPs_se.vcf"],
				sample = SAMPLES_SE)

elif config["illumina_pe"] is None and config["QC_Tests"] is False and config["Haplotyper"] is True and config["Validation"] is not None:
	rule all:
		input:
			HAPLOTYPE_DIR + "haplotypes.done",
			expand([MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat",  FINAL_VARIANTS + "Overlapping_SNPs_se", FINAL_VARIANTS + "Overlapping_SNPs_se.vcf"],
				sample = SAMPLES_SE)

#########################################################################################################


elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is False and config["Validation"] is not None and config["Flash"] is False:
	rule all:
		input:
		  	expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
		    	sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is False and config["Validation"] is not None and config["Flash"] is False:
	rule all:
		input:
			expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
				sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_pe"] is None and config["QC_Tests"] is True and config["Haplotyper"] is False and config["Validation"] is not None:
	rule all:
		input:
			expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.fastqc", fastQC + "{sample}.se", FINAL_VARIANTS + "Overlapping_SNPs_se", FINAL_VARIANTS + "Overlapping_SNPs_se.vcf"],
		  		sample = SAMPLES_SE)

elif config["illumina_pe"] is None and config["QC_Tests"] is False and config["Haplotyper"] is False and config["Validation"] is not None:
	rule all:
		input:
			expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_VARIANTS + "Overlapping_SNPs_se", FINAL_VARIANTS + "Overlapping_SNPs_se.vcf"],
				sample = SAMPLES_SE)

#########################################################################################################


elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is False and config["Validation"] is None and config["Flash"] is True:
	rule all:
		input:
		  	expand([FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}"],
		    	sample = SAMPLES_PE, pair = "1 2".split())
		  

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is False and config["Validation"] is None and config["Flash"] is True:
	rule all:
		input:
			expand([FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat"],
				sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is True and config["Validation"] is None and config["Flash"] is True:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}"],
		    sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is True and config["Validation"] is None and config["Flash"] is True:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat"],
			sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is True and config["Validation"] is not None and config["Flash"] is True:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
		    sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is True and config["Validation"] is not None and config["Flash"] is True:
	rule all:
		input:
		  HAPLOTYPE_DIR + "haplotypes.done",
		  expand([MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf", MICROHAPLOT + "{sample}.sorted.sam", MICROHAPLOT + "labels.txt", MAP_SAM + "{sample}.sorted.sam", FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
		  	sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is True and config["Haplotyper"] is False and config["Validation"] is not None and config["Flash"] is True:
	rule all:
		input:
		  	expand([ FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_MULTIQC_QUALIMAP, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", fastQC_raw + "{sample}.raw.{pair}.fastqc", fastQC + "{sample}.{pair}", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
		    	sample = SAMPLES_PE, pair = "1 2".split())

elif config["illumina_se"] is None and config["QC_Tests"] is False and config["Haplotyper"] is False and config["Validation"] is not None and config["Flash"] is True:
	rule all:
		input:
			expand([ FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed", FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_VARIANTS + "Hardy_Weinberg_values", FINAL_DOCS_MULTIQC_IDXSTATS, FINAL_DOCS_MULTIQC_FLAGSTAT, FINAL_DOCS_GENERALSTATS + "{sample}.generalStats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_VARIANTS + "Overlapping_SNPs_pe", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf", FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"],
				sample = SAMPLES_PE, pair = "1 2".split())
