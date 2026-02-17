include { SAMTOOLS_VIEW  } from '../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT  } from '..//modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ } from '..//modules/nf-core/samtools/fastq/main.nf'
include { PEAR           } from '..//modules/nf-core/pear/main.nf'
include { MERGE_READS    } from '..//modules/local/mergereads/main.nf'
include { HOTCOUNT       } from '..//modules/local/hotcount/main.nf'

workflow TARGETED {
    take:
    ch_samplesheet
    fasta
    queries
    gene

    main:
    SAMTOOLS_VIEW(
        ch_samplesheet,
        fasta.map { meta, fa -> tuple(meta, fa, [])},
        [],
        [],
    )

    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam,
        fasta,
        ""
    )

    SAMTOOLS_FASTQ(
        SAMTOOLS_SORT.out.bam,
        false,
    )

    // Combine fastq and singleton before branching
    ch_fastq_and_singleton = SAMTOOLS_FASTQ.out.fastq
        .join(SAMTOOLS_FASTQ.out.singleton)

    // Branch based on fastq content (keeping both fastq and singleton together)
    ch_branched = ch_fastq_and_singleton.branch {
        non_empty: it[1].any { f -> f.countLines() > 0 }  // it[1] is the fastq files
        empty: true
    }

    // Run PEAR only on non-empty fastq files
    PEAR(
        ch_branched.non_empty.map { meta, fastq, singleton -> tuple(meta, fastq) }
    )

    // For non-empty fastq: merge PEAR assembled with singleton from branched output
    ch_pear_with_singleton = PEAR.out.assembled
        .join(ch_branched.non_empty.map { meta, fastq, singleton -> tuple(meta, singleton) })

    MERGE_READS(
        ch_pear_with_singleton
    )

    // For empty fastq: use singleton files directly (skip PEAR and MERGE_READS)
    ch_singleton_only = ch_branched.empty
        .map { meta, fastq, singleton -> tuple(meta, singleton) }

    // Combine MERGE_READS output with singleton-only samples
    ch_final_fastq = MERGE_READS.out.merged.mix(ch_singleton_only)

    def query_list = file("${queries}/${gene}/*.txt")

    ch_queries = ch_samplesheet.map { meta, _cram, _crai ->
        def query = query_list.find { file -> file.name.startsWith(meta.design) }
        if (!query) {
            error("Could not find a query file for design ${meta.design} in the query directory (${queries})")
        }
        tuple(meta, query)
    }
    ch_hotcount_input = ch_final_fastq
        .join(ch_queries, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, fastq, query -> tuple(meta, query, fastq) }

    HOTCOUNT(
        ch_hotcount_input
    )

    emit:
    hotcount = HOTCOUNT.out.counts
}
