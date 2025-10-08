process MERGE_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/52/52ccce28d2ab928ab862e25aae26314d69c8e38bd41ca9431c67ef05221348aa/data'
        : 'community.wave.seqera.io/library/coreutils_grep_gzip_lbzip2_pruned:838ba80435a629f8'}"

    input:
    tuple val(meta), path(pear_assembled), path(samtools_singleton)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: merged
    tuple val("${task.process}"), val('cat'), eval("cat --version 2>&1 | head -n 1 | sed 's/^.* //'"), topic: versions, emit: versions_cat
    tuple val("${task.process}"), val('gunzip'), eval("gunzip --version 2>&1 | head -n 1 | sed 's/^.* //'"), topic: versions, emit: versions_gunzip
    tuple val("${task.process}"), val('gzip'), eval("gzip --version 2>&1 | head -n 1 | sed 's/^.* //'"), topic: versions, emit: versions_gzip

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip -c ${pear_assembled} > pear_unzipped.fastq
    gunzip -c ${samtools_singleton} > singleton_unzipped.fastq
    cat pear_unzipped.fastq singleton_unzipped.fastq | gzip > ${prefix}.merged.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -c < /dev/null > ${prefix}.merged.fastq.gz
    """
}
