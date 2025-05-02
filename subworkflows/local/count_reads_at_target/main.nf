include { SAMTOOLS_VIEW } from '../../../modules/nf-core/samtools/view/main.nf'


params.reference_fasta = '/home/guest/Documents/Internship/Data/genome/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna'

workflow COUNT_READS_AT_TARGET {
    take:
    ch_samplesheet
    
    main:
    ch_reference = Channel.value( [ [:], file(params.reference_fasta) ] )
    
    SAMTOOLS_VIEW(
        ch_samplesheet,
        ch_reference,
        Channel.value(file('/dev/null')),    // qname
        Channel.value('')   // index_format
    )

    emit:
    bam = SAMTOOLS_VIEW.out.bam
}