process HOTCOUNT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(query_file), path(assembled_fastq)

    output:
    tuple val(meta), path("*counts.txt"), emit: counts
    tuple val("${task.process}"), val('hotcount'), val("0.0"), topic: versions, emit: versions_hotcount

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    do_it_gz.sh ${query_file} ${assembled_fastq} > ${prefix}.counts.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.txt
    """
}
