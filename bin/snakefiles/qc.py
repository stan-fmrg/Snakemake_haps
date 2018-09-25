############################################################################################################

if config["illumina_se"] is None and config["QC_Tests"] is True:
    rule fastqc_raw_pe:
        input:
            RAW_DIR + "{sample}_{pair}.fq.gz",
        output:
            fastQC_raw + "{sample}.raw.{pair}.fastqc"
        shell:
            """
            fastqc --outdir {fastQC_raw} {input} > {output}
            """

elif config["illumina_pe"] is None and config["QC_Tests"] is True:
    rule fastqc_raw_se:
        input:
            RAW_DIR + "{sample}.se.fq.gz"
        output:
            fastQC_raw + "{sample}.raw.fastqc"
        shell:
            """
            fastqc --outdir {fastQC_raw} {input} > {output}
            """
else:
    pass


############################################################################################################
#"""Run trimmomatic in paired end mode to eliminate Illumina adaptors and remove
 #low quality regions and reads."""
############################################################################################################

if config["illumina_se"] is None:
    rule qc_trimmomatic_pe:
        """
        Run trimmomatic in paired end mode to eliminate Illumina adaptors and
        remove low quality regions and reads.
        Inputs _1 and _2 are piped through gzip/pigz.
        Outputs _1 and _2 are piped to gzip/pigz (level 9).
        Outputs _3 and _4 are compressed with the builtin compressor from
        Trimmomatic. Further on they are catted and compressed with gzip/pigz
        (level 9).
        Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
        header. It is done posterior to the trimming since the output comes
        slower than the input is read.
        Number of threads used:
            4 for trimmomatic
            2 for gzip inputs
            2 for gzip outputs
            Total: 8
        """
        input:
            forward = RAW_DIR + "{sample}_1.fq.gz",
            reverse = RAW_DIR + "{sample}_2.fq.gz"
        output:
            forward     = QC_DIR + "{sample}.final.1.fq.gz",
            reverse     = QC_DIR + "{sample}.final.2.fq.gz",
            unpaired_pe    = protected(QC_DIR + "{sample}.final.unpaired.fq.gz")
        params:
            unpaired_1  = QC_DIR + "{sample}_3.fq.gz",
            unpaired_2  = QC_DIR + "{sample}_4.fq.gz",
            adapter     = lambda wildcards: config["illumina_pe"][wildcards.sample]["adapter"],
            phred       = lambda wildcards: config["illumina_pe"][wildcards.sample]["phred"],
            trimmomatic_params = config["Read Trimming"]["trimmomatic_params"]
        log:
            QC_DOC + "trimmomatic_pe_{sample}.log"
        benchmark:
            QC_DOC + "trimmomatic_pe_{sample}.json"
        threads:
            24 # I've been able to work with pigz and 24 trimmomatic threads.
        shell:
            """
            trimmomatic PE \
                -threads {threads} \
                -{params.phred} \
                {input.forward} {input.reverse} \
                {output.forward} {params.unpaired_1} {output.reverse} {params.unpaired_2} \
                ILLUMINACLIP:{params.adapter}:2:30:10 \
                {params.trimmomatic_params} \
                2> {log}

                zcat {params.unpaired_1} {params.unpaired_2} |
                cut -f 1 -d " " |
                pigz -9 > {output.unpaired_pe}

            rm {params.unpaired_1} {params.unpaired_2}
            """

elif config["illumina_pe"] is None:
    rule qc_trimmomatic_se:
        """
        Run trimmomatic in paired end mode to eliminate Illumina adaptors and
        remove low quality regions and reads.
        Inputs _1 and _2 are piped through gzip/pigz.
        Outputs _1 and _2 are piped to gzip/pigz (level 9).
        Outputs _3 and _4 are compressed with the builtin compressor from
        Trimmomatic. Further on they are catted and compressed with gzip/pigz
        (level 9).
        Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
        header. It is done posterior to the trimming since the output comes
        slower than the input is read.
        Number of threads used:
            4 for trimmomatic
            2 for gzip inputs
            2 for gzip outputs
            Total: 8
        """
        input:
            single = RAW_DIR + "{sample}.se.fq.gz"
        output:
            single	= QC_DIR + "{sample}.se.final.fq.gz"
        params:
            adapter	= lambda wildcards: config["illumina_se"][wildcards.sample]["adapter"],
            phred	= lambda wildcards: config["illumina_se"][wildcards.sample]["phred"],
            trimmomatic_params = config["Read Trimming"]["trimmomatic_params"]
        log:
            QC_DOC + "trimmomatic_se_{sample}.log"
        benchmark:
            QC_DOC + "trimmomatic_se_{sample}.json"
        threads:
            24 # I've been able to work with pigz and 24 trimmomatic threads.
        shell:
            """
            trimmomatic SE \
            	-threads {threads} \
            	-{params.phred} \
            	{input.single} \
            	{output.single} \
            	ILLUMINACLIP:{params.adapter}:2:30:10 \
            	{params.trimmomatic_params} \
            	2> {log}
            """



############################################################################################################

if config["illumina_se"] is None and config["QC_Tests"] is True:
    rule fastqc_pe:
        input:
            QC_DIR + "{sample}.final.{pair}.fq.gz"
        output:
            fastQC + "{sample}.{pair}"
        shell:
            """
            fastqc --outdir {fastQC} {input} > {output}
            """

elif config["illumina_pe"] is None and config["QC_Tests"] is True:
    rule fastqc_se:
        input:
            QC_DIR + "{sample}.se.final.fq.gz"
        output:
            fastQC + "{sample}.se"
        shell:
            """
            fastqc --outdir {fastQC} {input} > {output}
            """
else:
	pass

############################################################################################################

if config["illumina_se"] is None and config["Flash"] is True:
    rule flash_pe_reads:
        input:
            forward = QC_DIR + "{sample}.final.1.fq.gz",
            reverse = QC_DIR + "{sample}.final.2.fq.gz"
        output:
            name = FLASHED_READS + "{sample}.flashed",
            file_name = FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz"
        params:
            scripts_dir = config["scripts_dir"],
            flash_params = config["FLASH"]["FLASH_params"]
        threads:
            4
        shell:
            "flash {params.flash_params} {input.forward} {input.reverse} -o {output.name} -z"
            #"{params.scripts_dir}/flash_reads.sh {input.forward} {input.reverse} {output.name}"
            "> {output.name} "   


elif config["illumina_pe"] is None:
    pass


############################################################################################################

if config["illumina_se"] is None:
    rule qc_results_pe:
        '''Generate only the resulting files, not the reports'''
        input:
            pe_files = expand(
                [QC_DIR + "{sample}.final.{pair}.fq.gz", fastQC_raw + "{sample}.raw.{pair}.fastqc", FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", FLASHED_READS + "{sample}.flashed"],
                sample = SAMPLES_PE,
                pair = "1 2".split()
            )

elif config["illumina_pe"] is None:
    rule qc_results_se:
        '''Generate only the resulting files, not the reports'''
        input:
            se_files = expand(
				[QC_DIR + "{sample}.se.final.fq.gz", fastQC_raw + "{sample}.raw_fastqc.zip"],
                sample = SAMPLES_SE
            )

############################################################################################################



if config["illumina_se"] is None and config["QC_Tests"] is True:
    rule raw_results_fastqc_pe:
        """Checkpoint to generate all the links for raw data"""
        input:
            expand(
                fastQC_raw + "{sample}.raw.{pair}.fastqc.{extension}", 
                sample = SAMPLES_PE,
                pair = "1 2".split(),
                extension = "zip html".split()
                )
elif config["illumina_pe"] is None and config["QC_Tests"] is True:
    rule raw_results_fastqc_se:
        """Checkpoint to generate all the links for raw data"""
        input:
            expand(
                fastQC_raw + "{sample}.raw.fastqc.{extension}",
                sample = SAMPLES_SE,
                extension = "zip html".split()
           )