
#rule flashed_reads:
#    input:
#        expand(FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz", sample = SAMPLES_PE)
#########################################################################################################

rule map_bwa_index_reference:
    input:
        RAW_DIR + "genome.fasta"
    output:
        RAW_DIR + "genome.fasta.bwt"
    log:
        MAP_DOC + "bwa_index_reference.log"
    benchmark:
        MAP_DOC + "bwa_index_reference.json"
    threads:
        4
    shell:
        "bwa index {input} &> {log}"

#########################################################################################################

if config["illumina_se"] is None and config["Flash"] is False:
   rule map_bwa_sample_pe:
        input:
            reference = RAW_DIR + "genome.fasta",
            index = RAW_DIR + "genome.fasta.bwt",
            forward = QC_DIR + "{sample}.final.1.fq.gz",
            reverse = QC_DIR + "{sample}.final.2.fq.gz"
        output:
            temp(MAP_DIR + "{sample}.unsorted.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}\tPL:Illumina",
            bwa_params = config["BWA-Mapping"]["BWA_mem_params"],
            samtools_params = config["Samtools_mapping_params"]["samtools_view_params"]
        log:
            MAP_DOC + "bwa_{sample}_pe.log"
        benchmark:
            MAP_DOC + "bwa_{sample}_pe.json"
        shell:
            "(bwa mem -R '{params.rg}' {params.bwa_params} {input.reference} {input.forward} {input.reverse} | "
            "samtools view {params.samtools_params} - "
            "> {output}) "
            "2> {log}"

elif config["illumina_se"] is None and config["Flash"] is True:
   rule map_bwa_sample_pe_flashed:
        input:
            reference = RAW_DIR + "genome.fasta",
            index = RAW_DIR + "genome.fasta.bwt",
            Flashed_read = FLASHED_READS + "{sample}.flashed.extendedFrags.fastq.gz"
        output:
            temp(MAP_DIR + "{sample}.unsorted.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}\tPL:Illumina",
            bwa_params = config["BWA-Mapping"]["BWA_mem_params"],
            samtools_params = config["Samtools_mapping_params"]["samtools_view_params"]
        log:
            MAP_DOC + "bwa_{sample}.log"
        benchmark:
            MAP_DOC + "bwa_{sample}.json"
        shell:
            "(bwa mem -R '{params.rg}' {params.bwa_params} {input.reference} {input.Flashed_read} | "
            "samtools view {params.samtools_params} - "
            "> {output}) "
            "2> {log}"


elif config["illumina_pe"] is None:
    rule map_bwa_sample_se:
        input:
            reference = RAW_DIR + "genome.fasta",
            index = RAW_DIR + "genome.fasta.bwt",
            single = QC_DIR + "{sample}.se.final.fq.gz"
        output:
            temp(MAP_DIR + "{sample}.unsorted.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}\tPL:Illumina",
            bwa_params = config["BWA-Mapping"]["BWA_mem_params"],
            samtools_params = config["Samtools_mapping_params"]["samtools_view_params"]
        log:
            MAP_DOC + "bwa_{sample}.log"
        benchmark:
            MAP_DOC + "bwa_{sample}.json"
        shell:
            "(bwa mem -R '{params.rg}' {params.bwa_params} {input.reference} {input.single} | "
            "samtools view {params.samtools_params} - "
            "> {output}) "
            "2> {log}"

#########################################################################################################


rule map_sort_sample:
    input:
        MAP_DIR + "{sample}.unsorted.bam"
    output:
        temp(MAP_DIR + "{sample}.Rmdup.bam")
    log:
        MAP_DOC + "sort_{sample}.log"
    benchmark:
        MAP_DOC + "sort_{sample}.json"
    threads:
        4
    shell:
        "samtools sort -n -@ {threads} "
            "-T $(mktemp --dry-run) "
            "-O bam {input} "
        "> {output} "
        "2> {log}"

#########################################################################################################


rule map_fixmate_sample:
    input:
        MAP_DIR + "{sample}.Rmdup.bam"
    output:
        temp(MAP_DIR + "{sample}.fixmate.bam")
    params:
        fixmate_params = config["SamTools_Fixmate"]["fixmate_params"]
    log:
        MAP_DOC + "fixmate_{sample}.log"
    benchmark:
        MAP_DOC + "fixmate_{sample}.json"
    shell:
        "samtools fixmate {params.fixmate_params} {input} {output} "




#########################################################################################################


rule map_sort_2_sample:
    input:
        MAP_DIR + "{sample}.fixmate.bam"
    output:
        temp(MAP_DIR + "{sample}.fixmate_sorted.bam")
    log:
        MAP_DOC + "sort_2_{sample}.log"
    benchmark:
        MAP_DOC + "sort_2_{sample}.json"
    threads:
        4
    shell:
        "samtools sort -@ {threads} "
            "-T $(mktemp --dry-run) "
            "-O bam {input} "
        "> {output} "
        "2> {log}"


#########################################################################################################


rule map_markdup_sample:
    input:
        MAP_DIR + "{sample}.fixmate_sorted.bam"
    output:
        BAM = MAP_DIR + "{sample}.sorted.bam",
        BAI = MAP_DIR + "{sample}.sorted.bam.bai"
    params:
        markdup_params = config["Marking_Duplicates"]["markdup_params"]
    log:
        MAP_DOC + "rmdup_{sample}.log"
    benchmark:
        MAP_DOC + "rmdup_{sample}.json"
    shell:
        """
        samtools markdup {params.markdup_params} {input} {output.BAM}
        samtools index {output.BAM} {output.BAI}
        """

#########################################################################################################


rule map_stats_sample:
    input:
        MAP_DIR + "{sample}.sorted.bam"
    output:
        FINAL_DOCS_GENERALSTATS + "{sample}.generalStats"
    log:
        MAP_DOC + "stats_{sample}.log"
    benchmark:
        MAP_DOC + "stats_{sample}.json"
    shell:
        "samtools stats {input} | "
        "grep ^SN | cut -f 2- > {output}"


#########################################################################################################


rule map_idxstats_sample:
    input:
        MAP_DIR + "{sample}.sorted.bam"
    output:
        idxstats = FINAL_DOCS_IDXSTATS + "{sample}.idxstats",
        flagstat = FINAL_DOCS_FLAGSTAT + "{sample}.flagstat"
    benchmark:
        MAP_DOC + "index_{sample}_idxstats.json"
    shell:
        """
        samtools idxstats {input} > {output.idxstats}
        samtools flagstat {input} > {output.flagstat}
        """


#########################################################################################################



rule map_QualiMap_bam_sample:
    input:
        MAP_DIR + "{sample}.sorted.bam"
    output:
        sample_header = temp(MAP_DIR + "{sample}.header.sam"),
        new_sample_header = temp(MAP_DIR + "{sample}.header.new.sam"),
        new_bam = temp(MAP_DIR + "{sample}.new.bam"),
        new_bam_index = temp(MAP_DIR + "{sample}.new.bam.bai"),
        qualimap_outfile = FINAL_DOCS_QUALIMAP + "{sample}.qualimap"
    benchmark:
        MAP_DOC + "qualimap_{sample}.json"
    shell:
        """
        samtools view -H {input} > {output.sample_header}
        head -n -3 {output.sample_header} > {output.new_sample_header}
        samtools reheader {output.new_sample_header} {input} > {output.new_bam}
        samtools index {output.new_bam} {output.new_bam_index}
        qualimap bamqc -bam {output.new_bam} -outdir {output.qualimap_outfile} -outfile {output.qualimap_outfile} -outformat html
        """

#########################################################################################################

if config["illumina_se"] is None:
    rule map_MultiQC_Samtools_sample_pe:
        input:
            idxstats = expand(FINAL_DOCS_IDXSTATS + "{sample}.idxstats", sample = SAMPLES_PE),
            flagstat = expand(FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", sample = SAMPLES_PE)
        output:
            MultiQC_idxstats = FINAL_DOCS_MULTIQC_IDXSTATS,
            MultiQC_flagstat = FINAL_DOCS_MULTIQC_FLAGSTAT
        shell:
            """
            multiqc -m samtools {input.idxstats} -o {output.MultiQC_idxstats}
            multiqc -m samtools {input.flagstat} -o {output.MultiQC_flagstat}
            """

elif config["illumina_pe"] is None:
    rule map_MultiQC_Samtools_se:
        input:
            idxstats = expand(FINAL_DOCS_IDXSTATS + "{sample}.idxstats", sample = SAMPLES_SE),
            flagstat = expand(FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", sample = SAMPLES_SE)
        output:
            MultiQC_idxstats = FINAL_DOCS_MULTIQC_IDXSTATS,
            MultiQC_flagstat = FINAL_DOCS_MULTIQC_FLAGSTAT
        shell:
            """
            multiqc -m samtools {input.idxstats} -o {output.MultiQC_idxstats}
            multiqc -m samtools {input.flagstat} -o {output.MultiQC_flagstat}
            """

#########################################################################################################

if config["illumina_se"] is None and config["QC_Tests"] is True:
    rule map_MultiQC_Qualimap_sample_pe:
        input:
            QUALIMAP = expand(FINAL_DOCS_QUALIMAP + "{sample}.qualimap", sample = SAMPLES_PE)
        output:
            MULTIQC_QUALIMAP = FINAL_DOCS_MULTIQC_QUALIMAP
        shell:
            """
            multiqc -m qualimap {input.QUALIMAP} -o {output.MULTIQC_QUALIMAP}
            """

elif config["illumina_pe"] is None and config["QC_Tests"] is True:
    rule map_MultiQC_Qualimap_sample_se:
        input:
            QUALIMAP = expand(FINAL_DOCS_QUALIMAP + "{sample}.qualimap", sample = SAMPLES_SE)
        output:
            MULTIQC_QUALIMAP = FINAL_DOCS_MULTIQC_QUALIMAP
        shell:
            """
            multiqc -m qualimap {input.QUALIMAP} -o {output.MULTIQC_QUALIMAP}
            """
else:
	pass

#########################################################################################################

if config["illumina_se"] is None:
    rule map_merge_bam_pe:
        input:
            bam = expand(MAP_DIR + "{sample}.sorted.bam",
                         sample = SAMPLES_PE),
            bai = expand(MAP_DIR + "{sample}.sorted.bam.bai",
                         sample = SAMPLES_PE)
        output:
            BAM_1 = MAP_DIR + "merged.bam",
            BAM_2 = BAM_FILES + "merged.bam",
            BAI_1 = BAM_FILES + "merged.bam.bai",
            BAI_2 = MAP_DIR + "merged.bam.bai",
            MD5_1 = MAP_DIR + "mergedMD5.txt",
            MD5_2 = BAM_FILES + "mergedMD5.txt"
        log:
            MAP_DOC + "merge_bam_pe.log"
        benchmark:
            MAP_DOC + "merge_bam_pe.json"
        shell:
            """
            samtools merge {output.BAM_1} {input.bam} > {log} 2>&1
            cp {output.BAM_1} {output.BAM_2}
            md5sum {output.BAM_1} > {output.MD5_1}
            md5sum {output.BAM_1} > {output.MD5_2}
            samtools index {output.BAM_1} {output.BAI_1}
            cp {output.BAI_1} {output.BAI_2}
            """

elif config["illumina_pe"] is None:
    rule map_merge_bam_se:
        input:
            bam = expand(MAP_DIR + "{sample}.sorted.bam",
                         sample = SAMPLES_SE),
            bai = expand(MAP_DIR + "{sample}.sorted.bam.bai",
                         sample = SAMPLES_SE)
        output:
            BAM_1 = MAP_DIR + "merged.bam",
            BAM_2 = BAM_FILES + "merged.bam",
            BAI_1 = BAM_FILES + "merged.bam.bai",
            BAI_2 = MAP_DIR + "merged.bam.bai",
            MD5_1 = MAP_DIR + "mergedMD5.txt",
            MD5_2 = BAM_FILES + "mergedMD5.txt"
        log:
            MAP_DOC + "merge_bam_se.log"
        benchmark:
            MAP_DOC + "merge_bam_se.json"
        shell:
            """
            samtools merge {output.BAM_1} {input.bam} > {log} 2>&1
            cp {output.BAM_1} {output.BAM_2}
            md5sum {output.BAM_1} > {output.MD5_1}
            md5sum {output.BAM_1} > {output.MD5_2}
            samtools index {output.BAM_1} {output.BAI_1}
            cp {output.BAI_1} {output.BAI_2}
            """

#########################################################################################################


if config["illumina_se"] is None and config["Flash"] is False:
    rule map_samtools_properly_paired_pe:
        input:
            MAP_DIR + "merged.bam"
        output:
            BAM_1 = MAP_DIR + "merged_properly_paired.bam",
            BAM_2 = BAM_FILES + "merged_properly_paired.bam"
        params:
            samtools_params = config["Samtools_mapping_params"]["samtools_view_params"],
            samtools_filtering_params = config["Samtools_read_filtering"]["samtools_view_filtering"]
        log:
            MAP_DOC + "properly_paired_merged.log"
        benchmark:
            MAP_DOC + "properly_paired_merged.json"
        shell:
            """
			samtools view {params.samtools_params} {params.samtools_filtering_params} {input} > {output.BAM_1}
            cp {output.BAM_1} {output.BAM_2}
            """

elif config["illumina_se"] is None and config["Flash"] is True:
    pass


elif config["illumina_pe"] is None:
    pass


#########################################################################################################

if config["illumina_se"] is None:
    rule properly_paired_index_pe:
        input:
            MAP_DIR + "merged_properly_paired.bam"
        output:
            BAI_1 = MAP_DIR + "merged_properly_paired.bam.bai",
            BAI_2 = BAM_FILES + "merged_properly_paired.bam.bai"
        log:
            MAP_DOC + "index_merged_properly_paired.log"
        benchmark:
            MAP_DOC + "index_merged_properly_paired.json"
        shell:
            """
            samtools index {input} {output.BAI_1}
            cp {output.BAI_1} {output.BAI_2}
            """

elif config["illumina_se"] is None and config["Flash"] is True:
    pass

elif config["illumina_pe"] is None:
    pass

#########################################################################################################

if config["illumina_se"] is None:
    rule map:
        input:
            MAP_DIR + "merged.bam.bai",
            expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_DOCS_GENERALSTATS + "{sample}.generalStats"],
                sample = SAMPLES_PE)

elif config["illumina_pe"] is None:
    rule map:
        input:
            MAP_DIR + "merged.bam.bai",
            expand([FINAL_DOCS_IDXSTATS + "{sample}.idxstats", FINAL_DOCS_FLAGSTAT + "{sample}.flagstat", FINAL_DOCS_GENERALSTATS + "{sample}.generalStats"],
                sample = SAMPLES_SE)