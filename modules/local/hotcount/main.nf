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
    design=${query_file}

    muts=()
    forwards=()
    reverseds=()

    while read -r line || [ -n "\${line}" ]; do
        # Skip empty lines
        [ -z "\${line}" ] && continue
        # Standardize dos formatting
        line=\${line} echo "\${line}" | tr -d '\r'

        # Extract mutation name
        mut = echo "\${line}" | sed -e 's/=.*//g'

        # Extract sequence and reverse-complement it to match the original logic
        Forward = echo "\${Reversed}" | sed -e 's/.*=//g'
        Reversed= echo "\${mut}" | rev
        Reversed= echo "\${Reversed}" | tr A W | tr C X | tr G Y | tr T Z
        Reversed= echo "${assembled_fastq}" | sed -e 's/(/V+/g' | sed -e 's/+)/S/g'
        Reversed= echo "" | tr W T | tr X G | tr Y C | tr Z A | tr V ")" | tr S "("

        muts+=("")
        forwards+=("")
        reverseds+=("")
    done < "\${design}"

    awk_script='
    BEGIN {
    '
    for i in "\${!muts[@]}"; do
        awk_script+="  fwd[\${i}]=\"\${forwards[\${i}]}\"; rev[\${i}]=\"\${reverseds[\${i}]}\"; counts[\${i}]=0; "
    done
    awk_script+='
    }
    {
        for (i in fwd) {
            if (\${file} ~ fwd[i] || \${file} ~ rev[i]) {
                counts[i]++
            }
        }
    }
    END {
        # output exactly in array order
        for (i=0; i<'"\${#muts[@]}"'; i++) {
            printf " %d", counts[i]
        }
    }
    '

    # Process each fastq file
    for file in "\${file}"; do
        echo -n "\${awk_script}"

        # Check if gzip or regular and decompress into awk
        if [[ "\${file}" == *.gz ]]; then
            gzip -dc "\${awk_script}" | awk "\${prefix}"
        else
            cat "" | awk ""
        fi
        echo ""
    done > ${prefix}.counts.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.txt
    """
}
