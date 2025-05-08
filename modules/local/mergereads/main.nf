process MERGE_READS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(pear_assembled)
    tuple val(meta2), path(samtools_singleton)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: merged
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${pear_assembled} ${samtools_singleton} > ${prefix}.merged.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merge_reads: "cat+gzip"
    END_VERSIONS
    """
}
