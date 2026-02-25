include { SAMTOOLS_VIEW             } from '../modules/nf-core/samtools/view/main.nf'
include { SAMTOOLS_SORT             } from '../modules/nf-core/samtools/sort/main.nf'
include { SAMTOOLS_FASTQ            } from '../modules/nf-core/samtools/fastq/main.nf'
include { PEAR                      } from '../modules/nf-core/pear/main.nf'
include { MERGE_READS               } from '../modules/local/mergereads/main.nf'
include { HOTCOUNT                  } from '../modules/local/hotcount/main.nf'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'

// Modules
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'

workflow TARGETED {
    take:
    ch_samplesheet
    fasta
    queries
    gene
    multiqc_config              // string: path to the multiqc config file
    multiqc_logo                // string: path to the multiqc logo
    multiqc_methods_description // string: path to the multiqc methods description

    main:
    def ch_multiqc_files = channel.empty()

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
    ch_branched = ch_fastq_and_singleton.branch {  _meta, fastq, _singleton ->
        non_empty: fastq.any { f -> f.countLines() > 0 }
        empty: true
    }

    // Run PEAR only on non-empty fastq files
    PEAR(
        ch_branched.non_empty.map { meta, fastq, _singleton -> tuple(meta, fastq) }
    )

    // For non-empty fastq: merge PEAR assembled with singleton from branched output
    ch_pear_with_singleton = PEAR.out.assembled
        .join(ch_branched.non_empty.map { meta, _fastq, singleton -> tuple(meta, singleton) })

    MERGE_READS(
        ch_pear_with_singleton
    )

    // For empty fastq: use singleton files directly (skip PEAR and MERGE_READS)
    ch_singleton_only = ch_branched.empty
        .map { meta, _fastq, singleton -> tuple(meta, singleton) }

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

    //
    // Collate and save software versions
    //
    def ch_collated_versions = softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${outdir}/${params.unique_out}",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // Perform multiQC on all QC data
    //

    def ch_multiqc_config                     = channel.fromPath(
                                                "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    def ch_multiqc_custom_config              = multiqc_config ?
                                                channel.fromPath(multiqc_config, checkIfExists: true) :
                                                channel.empty()
    def ch_multiqc_logo                       = multiqc_logo ?
                                                channel.fromPath(multiqc_logo, checkIfExists: true) :
                                                channel.empty()

    def summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary                   = channel.value(paramsSummaryMultiqc(summary_params))
    def ch_multiqc_custom_methods_description = multiqc_methods_description ?
                                                file(multiqc_methods_description, checkIfExists: true) :
                                                file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description                = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files                          = ch_multiqc_files.mix(
                                                    ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
                                                    ch_collated_versions,
                                                    ch_methods_description.collectFile(
                                                        name: 'methods_description_mqc.yaml',
                                                        sort: false
                                                    ),
                                                    ch_reports
                                                )
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    hotcount            = HOTCOUNT.out.counts
    multiqc_report      = MULTIQC.out.report            // channel: /path/to/multiqc_report.html
    multiqc_data        = MULTIQC.out.data              // channel: /path/to/multiqc_data
}
