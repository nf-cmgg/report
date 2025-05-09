process HOTCOUNT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(query_file), path(assembled_fastq)

    output:
    tuple val(meta), path("*counts.txt"), emit: counts
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    do_it_gz.sh ${query_file} ${assembled_fastq} > ${prefix}.counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hotcount: \$(echo "1.0.0") # Replace with actual version command if available
    END_VERSIONS
    """
}
