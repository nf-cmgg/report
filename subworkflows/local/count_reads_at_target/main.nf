include { SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { PEAR           } from '../../../modules/nf-core/pear/main.nf'
include { MERGE_READS    } from '../../../modules/local/mergereads/main.nf'
include { HOTCOUNT       } from '../../../modules/local/hotcount/main.nf'

workflow COUNT_READS_AT_TARGET {
    take:
    ch_samplesheet
    queries

    main:
    ch_reference = Channel.value([[:], file(params.fasta)])
    def ch_versions = Channel.empty()

    SAMTOOLS_VIEW(
        ch_samplesheet,
        ch_reference,
        [],
        [],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam,
        ch_reference,
        ""
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_FASTQ(
        SAMTOOLS_SORT.out.bam,
        false,
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    PEAR(
        SAMTOOLS_FASTQ.out.fastq
    )
    ch_versions = ch_versions.mix(PEAR.out.versions.first())
    ch_merge_input = PEAR.out.assembled.join(SAMTOOLS_FASTQ.out.singleton, failOnDuplicate: true, failOnMismatch: true)

    MERGE_READS(
        ch_merge_input
    )
    ch_versions = ch_versions.mix(MERGE_READS.out.versions.first())

    def query_list = file("${queries}/*.txt")

    ch_queries = ch_samplesheet.map { meta, _cram, _crai ->
        def query = query_list.find { file -> file.name.startsWith(meta.design) }
        if (!query) {
            error("Could not find a query file for design ${meta.design} in the query directory (${queries})")
        }
        tuple(meta, query)
    }
    ch_hotcount_input = MERGE_READS.out.merged
        .join(ch_queries, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, fastq, query -> tuple(meta, query, fastq) }

    HOTCOUNT(
        ch_hotcount_input
    )
    ch_versions = ch_versions.mix(HOTCOUNT.out.versions.first())

    emit:
    bam        = SAMTOOLS_VIEW.out.bam
    sorted_bam = SAMTOOLS_SORT.out.bam
    versions   = ch_versions
}
