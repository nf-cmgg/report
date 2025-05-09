include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { PEAR } from '../../../modules/nf-core/pear/main.nf'
include { MERGE_READS } from '../../../modules/local/mergereads/main.nf'
include { HOTCOUNT } from '../../../modules/local/hotcount/main.nf'

workflow COUNT_READS_AT_TARGET {
    take:
    ch_samplesheet

    main:
    ch_reference = Channel.value( [ [:], file(params.fasta) ] )

    SAMTOOLS_VIEW(
        ch_samplesheet.map { meta, cram, crai, _query -> tuple(meta, cram, crai)},
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

    ch_all_queries = Channel
    .fromPath("${params.queries_dir}/*.txt")
    .map { query_file ->
        def fname = query_file.getName()
        def design = fname.replaceFirst(/_[^_]+\.txt$/, '')
        tuple(design, query_file)
    }
    ch_samplesheet_design = ch_samplesheet.map { meta, _cram, _crai, design ->
        tuple(design, meta)
    }
    ch_queries = ch_samplesheet_design
        .join(ch_all_queries, failOnDuplicate:true, failOnMismatch:true)
        .map { _design, meta, query_file -> tuple(meta, query_file) }
    ch_hotcount_input = MERGE_READS.out.merged.join(ch_queries, failOnDuplicate:true, failOnMismatch:true).map { meta, fastq, query -> tuple(meta, query, fastq)}

    HOTCOUNT(
        ch_hotcount_input
    )

    emit:
    bam = SAMTOOLS_VIEW.out.bam
    sorted_bam = SAMTOOLS_SORT.out.bam
}
