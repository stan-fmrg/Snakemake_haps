if config["illumina_se"] is None and config["Haplotyper"] is True:
    rule map_BAM_to_SAM_pe:
        input:
            MAP_DIR + "{sample}.sorted.bam"
        output:
            MAP_SAM + "{sample}.sorted.sam"
        params:
            scripts_dir = config["scripts_dir"]
        shell:
            "samtools view -h -o {output} {input}"


elif config["illumina_pe"] is None and config["Haplotyper"] is True:
    rule map_BAM_to_SAM_se:
        input:
            MAP_DIR + "{sample}.sorted.bam"
        output:
            MAP_SAM + "{sample}.sorted.sam"
        shell:
            "samtools view -h -o {output} {input}"

else:
	pass

##########################################################################################

if config["illumina_se"] is None and config["Haplotyper"] is True:
    rule make_label:
        input:
        	expand(MAP_SAM + "{sample}.sorted.sam", sample = SAMPLES_PE)
        output:
            MICROHAPLOT + "labels.txt"
        params:
            scripts_dir = config["scripts_dir"]
        shell:
            "{params.scripts_dir}/Make_labels.sh {output}"

elif config["illumina_pe"] is None and config["Haplotyper"] is True:
    rule make_label:
        input:
            expand(MAP_SAM + "{sample}.sorted.sam", sample = SAMPLES_SE)
        output:
            MICROHAPLOT + "labels.txt"
        params:
            scripts_dir = config["scripts_dir"]
        shell:
            "{params.scripts_dir}/Make_labels.sh {output}"

else:
	pass

##########################################################################################

if config["Haplotyper"] is True:
    rule make_SAM_link:
        """
        Make symlinks to run "microhaplot"
        """
        input:
            MAP_SAM + "{sample}.sorted.sam"
        output:
            MICROHAPLOT + "{sample}.sorted.sam"
        log:
            HAPLOTYPE_DOC + "make_SAM_link.{sample}.log"
        threads:
            1
        shell:
            "ln -s $(readlink -f {input}) {output} 2> {log} "

else:
    pass 


########################################################################################


if config["Haplotyper"] is True and config["Validation"] is None:
    rule make_VCF_link:
        """
        Make symlinks to run "microhaplot"
        """
        input:
            vcf = FINAL_VARIANTS + "filtered_snps.recode.vcf"
        output:
            vcf = MICROHAPLOT + "filtered_snps.recode.vcf"
        log:
            HAPLOTYPE_DOC + "make_VCF_link.log"
        threads:
            1
        shell:
            "ln -s $(readlink -f {input.vcf}) {output.vcf} 2> {log} "

elif config["Haplotyper"] is True and config["Validation"] is not None:
    rule make_VCF_link:
        """
        Make symlinks to run "microhaplot"
        """
        input:
            vcf = FINAL_VARIANTS + "Overlapping_SNPs.vcf"
        output:
            vcf = MICROHAPLOT + "Overlapping_filtered_snps.recode.vcf"
        log:
            HAPLOTYPE_DOC + "make_VCF_link.log"
        threads:
            1
        shell:
            "ln -s $(readlink -f {input.vcf}) {output.vcf} 2> {log} "

else:
	pass