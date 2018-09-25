rule raw_make_link_reference:
    """
    Make a link to the reference genome.
    """
    input:
        config["genome"]
    output:
        RAW_DIR + "genome.fasta"
    log:
        RAW_DOC + "make_link_reference.log"
    benchmark:
        RAW_DOC + "make_link_reference.json"
    shell:
        "ln -s $(readlink -f {input}) {output} 2> {log}"



if config["illumina_se"] is None:
    rule raw_make_links_pe_sample:
        """
        Make a link next to the original file, with a prettier name than default.
        """
        input:
            forward = lambda wildcards: config["illumina_pe"][wildcards.sample]["forward"],
            reverse = lambda wildcards: config["illumina_pe"][wildcards.sample]["reverse"]
        output:
            forward = protected(RAW_DIR + "{sample}_1.fq.gz"),
            reverse = protected(RAW_DIR + "{sample}_2.fq.gz")
        log:
            RAW_DOC + "make_links_pe_{sample}.log"
        benchmark:
            RAW_DOC + "make_links_pe_{sample}.json"
        shell:
            "ln -s $(readlink -f {input.forward}) {output.forward} 2> {log};"
            "ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}"

elif config["illumina_pe"] is None:
    rule raw_make_links_se_sample:
        """
        Make a link next to the original file, with a prettier name than default.
        """
        input:
           single = lambda wildcards: config["illumina_se"][wildcards.sample]["single"]
        output:
            single = protected(RAW_DIR + "{sample}.se.fq.gz")
        log:
            RAW_DOC + "make_links_se_{sample}.log"
        benchmark:
            RAW_DOC + "make_links_se_{sample}.json"
        shell:
            "ln -s $(readlink -f {input.single}) {output.single} 2> {log}"


if config["illumina_se"] is None:
    rule raw_results_pe:
        """Checkpoint to generate all the links for raw data"""
        input:
            expand(
                RAW_DIR + "{sample}_{pair}.fq.gz",
                sample = SAMPLES_PE,
                pair = "1 2".split()
                )
elif config["illumina_pe"] is None:
    rule raw_results_se:
        """Checkpoint to generate all the links for raw data"""
        input:
            expand(
                RAW_DIR + "{sample}.se.fq.gz",
                sample = SAMPLES_SE
           )
