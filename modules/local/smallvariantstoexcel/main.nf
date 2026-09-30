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
    tuple val(meta), path("results/**.xlsx"), emit: excels
    tuple val(meta), path("results/**.txt"), emit: reports
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python

    script:
    """
    smallvariants_to_excel.py \\
        --api-data ${api_data} \\
        --samplesheet ${samplesheet} \\
        --output-dir results \\
        --run-name ${meta.id} \\
        --threshold-coverage ${threshold_coverage} \\
        --build "${build}" \\
        --variant-caller ${variant_caller} \\
        --runtype ${runtype}
    """

    stub:
    """
    mkdir -p results/_RAWdata
    touch results/_RAWdata_${meta.id}.xlsx
    touch results/_RAWdata/${meta.id}_variants.txt
    """
}
