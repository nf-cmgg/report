process EXOMECNV_TO_CSV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python:3.14.2--0562a0df3245213a"

    input:
    tuple val(meta), path(samplesheet)
    path(input_files)

    output:
    tuple val(meta), path("*.tsv"), emit: summaries
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    exomecnv_to_csv.py \\
        --samplesheet ${samplesheet} \\
        --output-dir .
    """

    stub:
    """
    touch cnvs_exomedepth_designs_summary.tsv
    touch cnvs_exomedepth_panels_summary.tsv
    """
}