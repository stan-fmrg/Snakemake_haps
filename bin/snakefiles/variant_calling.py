rule variant_call_split_genome:
    """
    Split the reference genome into separate components for parallel variant calling
    """
    input:
        genome = RAW_DIR + "genome.fasta"
    output:
        expand(
            RAW_DIR + "genome.{index}.fasta",
            index = map(lambda n: "%01d" % n, range(0, int(config["variant_call_params"]["split_genome"])))
            )
    threads:
        1
    params:
        split = config["variant_call_params"]["split_genome"]
    log:
        VARIANT_CALL_DOC + "variant_call_split_fasta.log"
    benchmark:
        VARIANT_CALL_DOC + "variant_call_split_fasta.json"
    shell:
        "pyfasta split -n {params.split} {input} 2> {log}"


rule variant_call_make_bed:
    """
    Generate a BED file for each of the split FASTA files
    """
    input:
        RAW_DIR + "genome.{index}.fasta"
    output:
        VARIANT_CALL_DIR + "genome.{index}.bed"
    params:
        scripts_dir = config["scripts_dir"]
    threads:
        2
    log:
        VARIANT_CALL_DOC + "variant_call_make_bed.{index}.log"
    benchmark:
        VARIANT_CALL_DOC + "variant_call_make_bed.{index}.json"
    shell:
        "{params.scripts_dir}/fa_to_bed.sh {input}"
        "> {output} "
        "2> {log}"


#########################################################################################################

if config["illumina_se"] is None and config["Flash"] is False:
    rule variant_call_subset_bam_pe:
        """
        Generate a BAM file for each of the genome partitions
        """
        input:
            bam = MAP_DIR + "merged_properly_paired.bam",
            bai = MAP_DIR + "merged_properly_paired.bam.bai",
            bed = VARIANT_CALL_DIR + "genome.{index}.bed"
        output:
            bam = VARIANT_CALL_DIR + "merged.{index}.bam",
            bai = VARIANT_CALL_DIR + "merged.{index}.bam.bai"
        params:
            scripts_dir = config["scripts_dir"]
        threads:
            2
        log:
            VARIANT_CALL_DOC + "variant_call_make_bam.log"
        benchmark:
            VARIANT_CALL_DOC + "variant_call_make_bam.json"
        shell:
            "samtools view -L {input.bed} -b {input.bam} "
            "> {output.bam} "
            "2> {log} && "
            "samtools index {output.bam}"

elif config["illumina_se"] is None and config["Flash"] is True:
    rule variant_call_flashed_bam_pe:
        """
        Generate a BAM file for each of the genome partitions
        """
        input:
            bam = MAP_DIR + "merged.bam",
            bai = MAP_DIR + "merged.bam.bai",
            bed = VARIANT_CALL_DIR + "genome.{index}.bed"
        output:
            bam = VARIANT_CALL_DIR + "merged.{index}.bam",
            bai = VARIANT_CALL_DIR + "merged.{index}.bam.bai"
        params:
            scripts_dir = config["scripts_dir"]
        threads:
            2
        log:
            VARIANT_CALL_DOC + "variant_call_make_bam.log"
        benchmark:
            VARIANT_CALL_DOC + "variant_call_make_bam.json"
        shell:
            "samtools view -L {input.bed} -b {input.bam} "
            "> {output.bam} "
            "2> {log} && "
            "samtools index {output.bam}"


elif config["illumina_pe"] is None:
    rule variant_call_subset_bam_se:
        """
        Generate a BAM file for each of the genome partitions
        """
        input:
            bam = MAP_DIR + "merged.bam",
            bai = MAP_DIR + "merged.bam.bai",
            bed = VARIANT_CALL_DIR + "genome.{index}.bed"
        output:
            bam = VARIANT_CALL_DIR + "merged.{index}.bam",
            bai = VARIANT_CALL_DIR + "merged.{index}.bam.bai"
        params:
            scripts_dir = config["scripts_dir"]
        threads:
            2
        log:
            VARIANT_CALL_DOC + "variant_call_make_bam.log"
        benchmark:
            VARIANT_CALL_DOC + "variant_call_make_bam.json"
        shell:
            "samtools view -L {input.bed} -b {input.bam} "
            "> {output.bam} "
            "2> {log} && "
            "samtools index {output.bam}"

#########################################################################################################


rule variant_call_run_freebayes:
    """
    Run FreeBayes
    """
    input:
        bam = VARIANT_CALL_DIR + "merged.{index}.bam",
        bai = VARIANT_CALL_DIR + "merged.{index}.bam.bai",
        fasta = RAW_DIR + "genome.{index}.fasta"
    output:
        VARIANT_CALL_DIR + "snps.{index}.vcf"
    params:
        fb_params = config["FreeBayes_variant_calling"]["freebayes_params"]
    threads:
        4
    log:
        VARIANT_CALL_DOC + "freebayes.{index}.log"
    benchmark:
        VARIANT_CALL_DOC + "freebayes.{index}.json"
    shell:
        "freebayes {params.fb_params} -f {input.fasta} {input.bam} > {output}"


#########################################################################################################


rule variant_call_bgzip_vcf:
    """
    Compress the VCF files with bgzip
    """
    input:
        VARIANT_CALL_DIR + "snps.{index}.vcf"
    output:
        VARIANT_CALL_DIR + "snps.{index}.vcf.gz"
    threads:
        1
    log:
        VARIANT_CALL_DOC + "bgzip_{index}.log"
    benchmark:
        VARIANT_CALL_DOC + "bgzip_{index}.json"
    shell:
        "bgzip -c {input} > {output}"


#########################################################################################################


rule variant_call_index_vcf:
    """
    Index the VCF files with tabix
    """
    input:
        VARIANT_CALL_DIR + "snps.{index}.vcf.gz"
    output:
        VARIANT_CALL_DIR + "snps.{index}.vcf.gz.tbi"
    threads:
        2
    log:
        VARIANT_CALL_DOC + "vcf_index_{index}.log"
    benchmark:
        VARIANT_CALL_DOC + "vcf_index_{index}.json"
    shell:
        "tabix {input}"


#########################################################################################################


rule variant_call_merge_vcf:
    """
    Merge all VCF files into a single file with bcftools merge
    """
    input:
        vcf = expand(
            VARIANT_CALL_DIR + "snps.{index}.vcf.gz",
            index = map(lambda n: "%01d" % n, range(0, int(config["variant_call_params"]["split_genome"])))),
	    index = expand(
            VARIANT_CALL_DIR + "snps.{index}.vcf.gz.tbi",
            index = map(lambda n: "%01d" % n, range(0, int(config["variant_call_params"]["split_genome"]))))
    output:
        VARIANTS_1 = VARIANT_CALL_DIR + "raw_snps.vcf.gz",
        VARIANTS_2 = RAW_VARIANTS + "raw_snps.vcf.gz"
    threads:
        10
    log:
        VARIANT_CALL_DOC + "bcftools_merge.log"
    benchmark:
        VARIANT_CALL_DOC + "bcftools_merge.json"
    shell:
        """
        bcftools concat -Oz {input.vcf} > {output.VARIANTS_1} 2> {log}
        cp {output.VARIANTS_1} {output.VARIANTS_2}
        """

#########################################################################################################


rule variant_call_filter_vcf:
    """
    Basic site quality filtering for the VCF file
    """
    input:
        VARIANT_CALL_DIR + "raw_snps.vcf.gz"
    output:
        FILTERED_VARIANTS_1 = VARIANT_CALL_DIR + "filtered_snps.recode.vcf",
        FILTERED_VARIANTS_2 = FINAL_VARIANTS + "filtered_snps.recode.vcf",
        HARDY_WEINGBERG_VARIANTS = FINAL_VARIANTS + "Hardy_Weinberg_values"
    threads:
        2
    params:
        vcftools_params = config["VCFtools"]["variant_filtering"]
    log:
        VARIANT_CALL_DOC + "variant_call_filter_vcf.log"
    benchmark:
        VARIANT_CALL_DOC + "variant_call_filter_vcf.json"
    shell:
        """
        vcftools --gzvcf {input} {params.vcftools_params} -c > {output.FILTERED_VARIANTS_1}
        vcftools --gzvcf {output.FILTERED_VARIANTS_1} --hardy -c > {output.HARDY_WEINGBERG_VARIANTS}
        cp {output.FILTERED_VARIANTS_1} {output.FILTERED_VARIANTS_2}
        """

#########################################################################################################

if config["illumina_se"] is None and config["Validation"] is not None:
    rule variant_validation_vcf:
        """
        Use a previous variant dataset to compare, extract and store common high-value SNPs between the two datasets
        """
        input:
            FILTERED_VARIANTS = FINAL_VARIANTS + "filtered_snps.recode.vcf",
            INPUT_VCF_FILE = config["Validation"]
        output:
            FINAL_VARIANTS + "Overlapping_SNPs"
        log:
            VARIANT_CALL_DOC + "variant_validation.log"
        benchmark:
            VARIANT_CALL_DOC + "variant_validation.json"
        shell:
            """
            vcftools --vcf {input.FILTERED_VARIANTS} --diff {input.INPUT_VCF_FILE} --diff-site -c > {output}
            """

else:
    pass

#########################################################################################################

if config["illumina_se"] is None and config["Validation"] is not None:
    rule intersect_vcf:
        input:
            FILTERED_VARIANTS = FINAL_VARIANTS + "filtered_snps.recode.vcf",
            INPUT_VCF_FILE = config["Validation"]
        output:
            FINAL_VARIANTS + "Overlapping_SNPs.vcf"
        log:
            VARIANT_CALL_DOC + "bedtools_intersect.log"
        benchmark:
            VARIANT_CALL_DOC + "bedtools_intersect.json"
        shell:
            """
            bedtools intersect -a {input.FILTERED_VARIANTS} -b {input.INPUT_VCF_FILE} -header > {output}
            """

else:
    pass



#########################################################################################################


#if config["illumina_se"] is None and config["Validation"] is True:
#	rule variant_validation_vcf_pe:
#		input:
#			FILTERED_VARIANTS = FINAL_VARIANTS + "filtered_snps_pe.recode.vcf",
#			INPUT_VCF_FILE = config["validation_file"]
#		output:
#			FINAL_VARIANTS + "Common_SNPs.txt"
#		log:
#			VARIANT_CALL_DOC + "variant_validation_pe.log"
#		benchmark:
##			VARIANT_CALL_DOC + "variant_validation_pe.json"
#		shell:
#			"SnpSift concordance -v {input.INPUT_VCF_FILE} {input.FILTERED_VARIANTS} > {output}" 
#			
#elif config["illumina_pe"] is None and config["Validation"] is True:
#	rule variant_validation_vcf_se:
#		"""
#		Use a previous variant dataset to compare, extract and store common high-value SNPs between the two datasets
#		"""
#		input:
##			FILTERED_VARIANTS = FINAL_VARIANTS + "filtered_snps_se.recode.vcf",
#			INPUT_VCF_FILE = config["validation_file"]
#		output:
#			FINAL_VARIANTS + "Common_SNPs.txt"
#		log:
#			VARIANT_CALL_DOC + "variant_validation_se.log"
#		benchmark:
#			VARIANT_CALL_DOC + "variant_validation_se.json"
#		shell:
#			"SnpSift concordance -v {input.INPUT_VCF_FILE} {input.FILTERED_VARIANTS} > {output}"	
#
#else:
#	pass

