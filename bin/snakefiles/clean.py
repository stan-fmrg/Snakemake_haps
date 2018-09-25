rule clean:
    shell:
        "rm -rf "
            "{RAW_DIR} {RAW_DOC} "
            "{MAP_DIR} {MAP_DOC} "
            "{QC_DIR} {QC_DOC} "
            "{ASSEMBLY_QC_DIR} {ASSEMBLY_QC_DOC} "

rule clean_raw:
    shell:
        "rm -rf {RAW_DIR} {RAW_DOC} "

rule clean_map:
    shell:
        "rm -rf {MAP_DIR} {MAP_DOC} "

rule clean_qc:
    shell:
        "rm -rf {QC_DIR} {QC_DOC} "

rule clean_assembly_qc:
    shell:
        "rm -rf {ASSEMBLY_QC_DIR} {ASSEMBLY_QC_DOC} "