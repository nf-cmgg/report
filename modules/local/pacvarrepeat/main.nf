process PACVAR_REPEAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92f382f0f3c2af60bdaa142a038896267c45fa6c34ee6a9b04a3eae5b4bb1639/data'
        : 'community.wave.seqera.io/library/pip_openpyxl_pandas:9248bb616642565e' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.xlsx"), emit: excels

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    pacvar_repeat_xlsx_report.py \\
        --sample_name ${prefix} \\
        --vcf_file ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.pacvar_repeat_report.xlsx
    """
}
