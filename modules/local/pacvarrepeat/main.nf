process PACVAR_REPEAT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pip_openpyxl_pandas:d73ced3a164f683e'
        : 'community.wave.seqera.io/library/pip_openpyxl_pandas:9248bb616642565e' }"

    publishDir "${dir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(vcf), val(dir)

    output:
    tuple val(meta), path("*.xlsx"), emit: excels

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 ${projectDir}/bin/pacvar_repeat_xlsx_report.py \\
        --sample_name ${prefix} \\
        --vcf_file ${vcf}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """ 
    touch ${prefix}.pacvar_repeat_report.xlsx
    """
}
