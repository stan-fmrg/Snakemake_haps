if config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_make_vcf_link_pe:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = FINAL_VARIANTS + "filtered_snps_pe.recode.vcf"
        output:
            vcf = HAPLOTYPE_DIR + "filtered_snps_pe.recode.vcf"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_vcf_link.log"
        benchmark:
            HAPLOTYPE_DOC + "make_vcf_link.json"
        shell:
            "ln -s $(readlink -f {input.vcf}) {output.vcf} 2> {log} "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_make_vcf_link_se:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = FINAL_VARIANTS + "filtered_snps_se.recode.vcf"
        output:
            vcf = HAPLOTYPE_DIR + "filtered_snps_se.recode.vcf"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_vcf_link.log"
        benchmark:
            HAPLOTYPE_DOC + "make_vcf_link.json"
        shell:
            "ln -s $(readlink -f {input.vcf}) {output.vcf} 2> {log} "


elif config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_make_overlapping_vcf_link_pe:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = FINAL_VARIANTS + "Overlapping_SNPs_pe.vcf"
        output:
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_pe.recode.vcf"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_vcf_link.log"
        benchmark:
            HAPLOTYPE_DOC + "make_vcf_link.json"
        shell:
            "ln -s $(readlink -f {input.vcf}) {output.vcf} 2> {log} "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_make_overlapping_vcf_link_se:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = FINAL_VARIANTS + "Overlapping_SNPs_se.vcf"
        output:
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_se.recode.vcf"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_vcf_link.log"
        benchmark:
            HAPLOTYPE_DOC + "make_vcf_link.json"
        shell:
            "ln -s $(readlink -f {input.vcf}) {output.vcf} 2> {log} "

else:
	pass

#########################################################################################################

if config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_make_bam_links_pe:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "filtered_snps_pe.recode.vcf",
            bam = MAP_DIR + "{sample}.sorted.bam",
            bai = MAP_DIR + "{sample}.sorted.bam.bai"
        output:
            bam = HAPLOTYPE_DIR + "{sample}.bam",
            bai = HAPLOTYPE_DIR + "{sample}.bam.bai"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.log"
        benchmark:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.json"
        shell:
            "ln -s $(readlink -f {input.bam}) {output.bam} 2> {log}; "
            "ln -s $(readlink -f {input.bai}) {output.bai} 2> {log}; "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_make_bam_links_se:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "filtered_snps_se.recode.vcf",
            bam = MAP_DIR + "{sample}.sorted.bam",
            bai = MAP_DIR + "{sample}.sorted.bam.bai"
        output:
            bam = HAPLOTYPE_DIR + "{sample}.bam",
            bai = HAPLOTYPE_DIR + "{sample}.bam.bai"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.log"
        benchmark:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.json"
        shell:
            "ln -s $(readlink -f {input.bam}) {output.bam} 2> {log}; "
            "ln -s $(readlink -f {input.bai}) {output.bai} 2> {log}; "

elif config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_make_bam_links_overlapping_SNPs_pe:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_pe.recode.vcf",
            bam = MAP_DIR_PE + "{sample}.sorted.bam",
            bai = MAP_DIR_PE + "{sample}.sorted.bam.bai"
        output:
            bam = HAPLOTYPE_DIR + "{sample}.bam",
            bai = HAPLOTYPE_DIR + "{sample}.bam.bai"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.log"
        benchmark:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.json"
        shell:
            "ln -s $(readlink -f {input.bam}) {output.bam} 2> {log}; "
            "ln -s $(readlink -f {input.bai}) {output.bai} 2> {log}; "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_make_bam_links_overlapping_SNPs_se:
        """
        Make symlinks to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_se.recode.vcf",
            bam = MAP_DIR_SE + "{sample}.sorted.bam",
            bai = MAP_DIR_SE + "{sample}.sorted.bam.bai"
        output:
            bam = HAPLOTYPE_DIR + "{sample}.bam",
            bai = HAPLOTYPE_DIR + "{sample}.bam.bai"
        threads:
            1
        log:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.log"
        benchmark:
            HAPLOTYPE_DOC + "make_bam_links.{sample}.json"
        shell:
            "ln -s $(readlink -f {input.bam}) {output.bam} 2> {log}; "
            "ln -s $(readlink -f {input.bai}) {output.bai} 2> {log}; "


else:
	pass
#########################################################################################################

if config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_make_popmap_pe:
        """
        Make popmap to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "filtered_snps_pe.recode.vcf",
            bam = expand(HAPLOTYPE_DIR + "{sample}.bam",
                         sample = SAMPLES_PE)
        output:
            popmap = HAPLOTYPE_DIR + "popmap"
        threads:
            1
        params:
            hap_dir = HAPLOTYPE_DIR
        log:
            HAPLOTYPE_DOC + "make_popmap.log"
        benchmark:
            HAPLOTYPE_DOC + "make_popmap.json"
        shell:
            "ls {params.hap_dir}/*.bam | sed 's/.bam//' | awk '{{print $1, 1}}' > {output.popmap} "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_make_popmap_se:
        """
        Make popmap to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "filtered_snps_se.recode.vcf",
            bam = expand(HAPLOTYPE_DIR + "{sample}.bam",
                         sample = SAMPLES_SE)
        output:
            popmap = HAPLOTYPE_DIR + "popmap"
        threads:
            1
        params:
            hap_dir = HAPLOTYPE_DIR
        log:
            HAPLOTYPE_DOC + "make_popmap.log"
        benchmark:
            HAPLOTYPE_DOC + "make_popmap.json"
        shell:
            "ls {params.hap_dir}/*.bam | sed 's/.bam//' | awk '{{print $1, 1}}' > {output.popmap} "


elif config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_make_popmap_overlapping_SNPs_pe:
        """
        Make popmap to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_pe.recode.vcf",
            bam = expand(HAPLOTYPE_DIR + "{sample}.bam",
                         sample = SAMPLES_PE)
        output:
            popmap = HAPLOTYPE_DIR + "popmap"
        threads:
            1
        params:
            hap_dir = HAPLOTYPE_DIR
        log:
            HAPLOTYPE_DOC + "make_popmap.log"
        benchmark:
            HAPLOTYPE_DOC + "make_popmap.json"
        shell:
            "ls {params.hap_dir}/*.bam | sed 's/.bam//' | awk '{{print $1, 1}}' > {output.popmap} "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_make_popmap_overlapping_SNPs_se:
        """
        Make popmap to run rad_haplotyper
        """
        input:
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_se.recode.vcf",
            bam = expand(HAPLOTYPE_DIR + "{sample}.bam",
                         sample = SAMPLES_SE)
        output:
            popmap = HAPLOTYPE_DIR + "popmap"
        threads:
            1
        params:
            hap_dir = HAPLOTYPE_DIR
        log:
            HAPLOTYPE_DOC + "make_popmap.log"
        benchmark:
            HAPLOTYPE_DOC + "make_popmap.json"
        shell:
            "ls {params.hap_dir}/*.bam | sed 's/.bam//' | awk '{{print $1, 1}}' > {output.popmap} "


else:
	pass

#########################################################################################################


if config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_run_rad_haplotyper_pe:
        """
        Run rad_haplotyper
        """
        input:
            popmap = HAPLOTYPE_DIR + "popmap",
            vcf = HAPLOTYPE_DIR + "filtered_snps_pe.recode.vcf"
        output:
            output_1 = HAPLOTYPE_DIR + "haplotypes.done",
            output_2 = FINAL_HAPLOTYPES + "haplotypes.done"
        threads:
            20
        params:
            hap_dir = HAPLOTYPE_DIR,
            hap_dir_final_results = FINAL_HAPLOTYPES,
            home_dir = "../..",
            Haplotyping_params = config["Haplotyping_params"]["rad_haplotyper_params"]
        log:
            HAPLOTYPE_DOC + "rad_haplotyper.log"
        benchmark:
            HAPLOTYPE_DOC + "rad_haplotyper.json"
        shell:
            "cd {params.hap_dir}; "
            "rad_haplotyper.pl -v filtered_snps_pe.recode.vcf -x {threads} -g haps.gen -p popmap {params.Haplotyping_params}; "
            "cd {params.home_dir}; "
            "touch {output.output_1}; "
            "cp {output.output_1} {output.output_2}; "
            "cp results/Haplotype/codes.haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/ind_stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/hap_loci.txt results/Final_results/Haplotypes/; "

elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is None:
    rule haplotype_run_rad_haplotyper_se:
        """
        Run rad_haplotyper
        """
        input:
            popmap = HAPLOTYPE_DIR + "popmap",
            vcf = HAPLOTYPE_DIR + "filtered_snps_se.recode.vcf"
        output:
            output_1 = HAPLOTYPE_DIR + "haplotypes.done",
            output_2 = FINAL_HAPLOTYPES + "haplotypes.done"
        threads:
            20
        params:
            hap_dir = HAPLOTYPE_DIR,
            hap_dir_final_results = FINAL_HAPLOTYPES,
            home_dir = "../..",
            Haplotyping_params = config["Haplotyping_params"]["rad_haplotyper_params"]
        log:
            HAPLOTYPE_DOC + "rad_haplotyper.log"
        benchmark:
            HAPLOTYPE_DOC + "rad_haplotyper.json"
        shell:
            "cd {params.hap_dir}; "
            "rad_haplotyper.pl -v filtered_snps_se.recode.vcf -x {threads} -g haps.gen -p popmap {params.Haplotyping_params}; "
            "cd {params.home_dir}; "
            "touch {output.output_1}; "
            "cp {output.output_1} {output.output_2}; "
            "cp results/Haplotype/codes.haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/ind_stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/hap_loci.txt results/Final_results/Haplotypes/; "


elif config["illumina_se"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_run_rad_haplotyper_overlapping_SNPs_pe:
        """
        Run rad_haplotyper
        """
        input:
            popmap = HAPLOTYPE_DIR + "popmap",
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_pe.recode.vcf"
        output:
            output_1 = HAPLOTYPE_DIR + "haplotypes.done",
            output_2 = FINAL_HAPLOTYPES + "haplotypes.done"
        threads:
            20
        params:
            hap_dir = HAPLOTYPE_DIR,
            hap_dir_final_results = FINAL_HAPLOTYPES,
            home_dir = "../..",
            Haplotyping_params = config["Haplotyping_params"]["rad_haplotyper_params"]
        log:
            HAPLOTYPE_DOC + "rad_haplotyper.log"
        benchmark:
            HAPLOTYPE_DOC + "rad_haplotyper.json"
        shell:
            "cd {params.hap_dir}; "
            "rad_haplotyper.pl -v Overlapping_filtered_snps_pe.recode.vcf -x {threads} -g haps.gen -p popmap {params.Haplotyping_params}; "
            "cd {params.home_dir}; "
            "touch {output.output_1}; "
            "cp {output.output_1} {output.output_2}; "
            "cp results/Haplotype/codes.haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/ind_stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/hap_loci.txt results/Final_results/Haplotypes/; "


elif config["illumina_pe"] is None and config["Haplotyper"] is True and config["Validation"] is not None:
    rule haplotype_run_rad_haplotyper_overlapping_SNPs_se:
        """
        Run rad_haplotyper
        """
        input:
            popmap = HAPLOTYPE_DIR + "popmap",
            vcf = HAPLOTYPE_DIR + "Overlapping_filtered_snps_se.recode.vcf"
        output:
            output_1 = HAPLOTYPE_DIR + "haplotypes.done",
            output_2 = FINAL_HAPLOTYPES + "haplotypes.done"
        threads:
            20
        params:
            hap_dir = HAPLOTYPE_DIR,
            hap_dir_final_results = FINAL_HAPLOTYPES,
            home_dir = "../..",
            Haplotyping_params = config["Haplotyping_params"]["rad_haplotyper_params"]
        log:
            HAPLOTYPE_DOC + "rad_haplotyper.log"
        benchmark:
            HAPLOTYPE_DOC + "rad_haplotyper.json"
        shell:
            "cd {params.hap_dir}; "
            "rad_haplotyper.pl -v Overlapping_filtered_snps_se.recode.vcf -x {threads} -g haps.gen -p popmap {params.Haplotyping_params}; "
            "cd {params.home_dir}; "
            "touch {output.output_1}; "
            "cp {output.output_1} {output.output_2}; "
            "cp results/Haplotype/codes.haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/ind_stats.out results/Final_results/Haplotypes/; "
            "cp results/Haplotype/haps.gen results/Final_results/Haplotypes/; "
            "cp results/Haplotype/hap_loci.txt results/Final_results/Haplotypes/; "


else:
	pass



###########################################