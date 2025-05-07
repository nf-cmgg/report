include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT } from '../../../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main.nf'
include { PEAR } from '../../../modules/nf-core/pear/main.nf'

workflow COUNT_READS_AT_TARGET {
    take:
    ch_samplesheet

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

    emit:
    bam = SAMTOOLS_VIEW.out.bam
    sorted_bam = SAMTOOLS_SORT.out.bam
}
