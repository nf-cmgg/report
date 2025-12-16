process VARCOV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/59cf2a3b7b08e2c3b188333bef269e49a44da3bc46d2bbe46407fc5d59e275dd/data'
        : 'community.wave.seqera.io/library/samtools_pip_cyvcf2_openpyxl_pandas:5f705f003fd04932'}"

    input:
    tuple val(meta), path(vcfs, stageAs: 'vcfs/*'), path(stringtie, stageAs: 'stringtie/*'), path(fusionreport, stageAs: 'fusionreport/*'), path(ctat, stageAs: 'ctat/*'), path(multiqc, stageAs: 'multiqc/*'), path(bams, stageAs: 'bams/*'), path(bais, stageAs: 'bams/*'), val(run_nr)
    path genes
    path fusions
    path mane
    val pipeline_version

    output:
    tuple val(meta), path("*.xlsx"), emit: output
    tuple val("${task.process}"), val('python'), eval("python --version 2>&1 | sed 's/^Python //'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('pandas'), eval("pip freeze | grep pandas | sed 's/pandas==//'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('openpyxl'), eval("pip freeze | grep openpyxl | sed 's/openpyxl==//'"), topic: versions, emit: versions_openpyxl
    tuple val("${task.process}"), val('cyvcf2'), eval("pip freeze | grep cyvcf2 | sed 's/cyvcf2==//'"), topic: versions, emit: versions_cyvcf2
    tuple val("${task.process}"), val('samtools'), eval("samtools --version 2>&1 | head -n 1 | sed 's/^.* //'"), topic: versions, emit: versions_samtools

    script:
    """
    rnafusion_varcov.py \\
        --input vcfs \\
        --stringtie stringtie \\
        --fusionreport fusionreport \\
        --ctat ctat \\
        --multiqc multiqc \\
        --output . \\
        --bams bams \\
        --genes ${genes} \\
        --fusion_whitelist ${fusions} \\
        --mane ${mane} \\
        --run ${run_nr} \\
        --pipeline_version ${pipeline_version}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.xlsx
    """
}
