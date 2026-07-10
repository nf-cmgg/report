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
    design="${query_file}"
    assembled_file="${assembled_fastq}"
    sample_name=\$(basename "\${assembled_file}")
    run_date=\$(date)

    {
        echo "## script do_it.sh running in version 0.0"
        echo "## date is \${run_date}"
        echo -n "Sample"
        while IFS= read -r line || [ -n "\${line}" ]; do
            line="\${line%\$'\r'}"
            [ -z "\${line}" ] && continue
            mut="\${line%%=*}"
            echo -n " \${mut}"
        done < "\${design}"
        echo ""

        echo -n "\${sample_name}"
        while IFS= read -r line || [ -n "\${line}" ]; do
            line="\${line%\$'\r'}"
            [ -z "\${line}" ] && continue
            forward="\${line#*=}"
            [ -z "\${forward}" ] && continue
            reversed=\$(printf '%s' "\${forward}" | rev | tr 'ACGT' 'TGCA')

            if [[ "\${assembled_file}" == *.gz ]]; then
                nb=\$(zgrep -c -E "\${forward}|\${reversed}" "\${assembled_file}" || true)
            else
                nb=\$(grep -c -E "\${forward}|\${reversed}" "\${assembled_file}" || true)
            fi
            nb="\${nb:-0}"
            echo -n " \${nb}"
        done < "\${design}"
        echo ""
    } > "${prefix}.counts.txt"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.txt
    """
}
