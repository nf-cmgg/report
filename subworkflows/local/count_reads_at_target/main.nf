include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main.nf'

workflow COUNT_READS_AT_TARGET {
    take:
    ch_samplesheet

    main:
    ch_reference = Channel.value( [ [:], file(params.fasta) ] )

    SAMTOOLS_VIEW(
        ch_samplesheet,
        ch_reference,
        [],    // qname
        "bai"   // index_format
    )

    emit:
    bam = SAMTOOLS_VIEW.out.bam
}
