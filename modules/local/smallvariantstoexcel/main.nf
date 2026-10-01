process SMALLVARIANTS_TO_EXCEL {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker.io/library/python:3.14-bookworm"

    input:
    tuple val(meta), path(samplesheet)
    path(api_data)
    path(input_files)
    val threshold_coverage
    val build
    val variant_caller
    val runtype

    output:
    tuple val(meta), path("results_${variant_caller}/**.xlsx"), emit: excels
    tuple val(meta), path("results_${variant_caller}/**.txt"), emit: reports
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    smallvariants_to_excel.py \\
        --api-data ${api_data} \\
        --samplesheet ${samplesheet} \\
        --output-dir "results_${variant_caller}" \\
        --run-name ${meta.id} \\
        --threshold-coverage ${threshold_coverage} \\
        --build "${build}" \\
        --variant-caller ${variant_caller} \\
        --runtype ${runtype}
    """

    stub:
    """
    mkdir -p results_${variant_caller}/_RAWdata
    touch results_${variant_caller}/_RAWdata_${meta.id}.xlsx
    touch results_${variant_caller}/_RAWdata/${meta.id}_variants.txt
    """
}
