include { SAMTOOLS_VIEW  } from '../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT  } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ } from '../modules/nf-core/samtools/fastq/main.nf'
include { PEAR           } from '../modules/nf-core/pear/main.nf'
include { CAT_FASTQ      } from '../modules/nf-core/cat/fastq/main.nf'
include { HOTCOUNT       } from '../modules/local/hotcount/main.nf'

workflow TARGETED {
    take:
    ch_samplesheet // channel: sampleinfo => meta, cram, crai
    fasta_fai      // value: meta, fasta, fai
    queries        // string: path to the directory containing the query files
    gene           // string: gene name (used to select the query file)

    main:

    SAMTOOLS_VIEW(
        ch_samplesheet,
        fasta_fai,
        [],
        [],
    )

    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam,
        fasta_fai,
        "",
    )

    SAMTOOLS_FASTQ(
        SAMTOOLS_SORT.out.bam,
        false,
    )

    // Combine fastq and singleton before branching
    def ch_fastq_and_singleton = SAMTOOLS_FASTQ.out.fastq.join(SAMTOOLS_FASTQ.out.singleton)

    // Branch based on fastq content (keeping both fastq and singleton together)
    def ch_branched = ch_fastq_and_singleton.branch { _meta, fastq, _singleton ->
        non_empty: fastq.any { f -> f.size() > 31 } // empty gzipped files are 31 bytes
        empty: true
    }

    // Run PEAR only on non-empty fastq files
    PEAR(
        ch_branched.non_empty.map { meta, fastq, _singleton -> tuple(meta, fastq) }
    )

    // For non-empty fastq: merge PEAR assembled with singleton from branched output
    // CAT_FASTQ expects tuple(meta, reads) where reads is a path or list of paths.
    def ch_pear_with_singleton = PEAR.out.assembled
        .join(ch_branched.non_empty.map { meta, _fastq, singleton -> tuple(meta, singleton) })
        // PEAR assembled and singleton are single-end reads.
        .map { meta, assembled, singleton -> tuple(meta + [single_end: true], [assembled, singleton]) }

    CAT_FASTQ(ch_pear_with_singleton)

    // For empty fastq: use singleton files directly (skip PEAR and CAT_FASTQ)
    def ch_singleton_only = ch_branched.empty.map { meta, _fastq, singleton -> tuple(meta, singleton) }

    // Combine CAT_FASTQ output with singleton-only samples
    def ch_final_fastq = CAT_FASTQ.out.reads.mix(ch_singleton_only)

    def query_list = files("${queries}/${gene}/*.txt")

    def ch_hotcount_input = ch_final_fastq.map { meta, fastq ->
        def query = query_list.find { file -> file.name.startsWith(meta.design) }
        if (!query) {
            error("Could not find a query file for design ${meta.design} in the query directory (${queries})")
        }
        tuple(meta, query, fastq)
    }

    HOTCOUNT(
        ch_hotcount_input
    )

    emit:
    HOTCOUNT.out.counts
}
