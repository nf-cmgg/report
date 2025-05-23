include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { PEAR } from '../../../modules/nf-core/pear/main.nf'
include { MERGE_READS } from '../../../modules/local/mergereads/main.nf'
include { HOTCOUNT } from '../../../modules/local/hotcount/main.nf'
include { MULTIQC } from '../../../modules/nf-core/multiqc/main.nf'

workflow COUNT_READS_AT_TARGET {
    take:
    ch_samplesheet
    queries
    multiqc_config

    main:
    ch_reference = Channel.value( [ [:], file(params.fasta) ] )

    SAMTOOLS_VIEW(
        ch_samplesheet,
        ch_reference,
        [],    // qname
        []   // index_format
    )
    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam,
        ch_reference
    )
    SAMTOOLS_FASTQ(
        SAMTOOLS_SORT.out.bam,
        false
    )
    PEAR(
        SAMTOOLS_FASTQ.out.fastq
    )

    ch_merge_input = PEAR.out.assembled.join(SAMTOOLS_FASTQ.out.singleton, failOnDuplicate:true, failOnMismatch:true)

    MERGE_READS(
        ch_merge_input
    )

    def query_list = file("${queries}/*.txt")

    ch_queries = ch_samplesheet.map { meta, _cram, _crai ->
        def query = query_list.find { file -> file.name.startsWith(meta.design) }
        if(!query) {
            error("Could not find a query file for design ${meta.design} in the query directory (${queries})")
        }
        tuple(meta,query)
    }
    ch_hotcount_input = MERGE_READS.out.merged
        .join(ch_queries, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, fastq, query -> tuple(meta, query, fastq)}

    HOTCOUNT(
        ch_hotcount_input
    )

    ch_multiqc_input = HOTCOUNT.out.counts.collect{_meta, file -> file}

    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    ch_multiqc_custom_config = multiqc_config ?
        Channel.fromPath(multiqc_config, checkIfExists: true) :
        Channel.empty()

    MULTIQC (
        ch_multiqc_input,
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        [],
        []
    )

    emit:
    bam = SAMTOOLS_VIEW.out.bam
    sorted_bam = SAMTOOLS_SORT.out.bam
    multiqc_report = MULTIQC.out.report
}
