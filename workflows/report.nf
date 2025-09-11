/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNTAR                  } from '../modules/nf-core/untar/main.nf'
include { VARCOV                 } from '../modules/local/varcov/main.nf'
include { COUNT_READS_AT_TARGET  } from '../subworkflows/local/count_reads_at_target/main.nf'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_report_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REPORT {
    take:
    ch_samplesheet           // channel: samplesheet read in from --input
    queries
    ch_rnafusion_samplesheet // channel: samplesheet read in from --rnafusion_input
    ch_genes
    ch_fusions
    ch_mane

    main:
    def ch_versions = Channel.empty()

    //
    // Read counts
    //

    COUNT_READS_AT_TARGET(ch_samplesheet, queries)
    ch_versions = ch_versions.mix(COUNT_READS_AT_TARGET.out.versions)

    //
    // Rnafusion report flow
    //

    def ch_rnafusion_branch = ch_rnafusion_samplesheet.branch { _meta, dir ->
        tarzipped: dir.extension == "gz"
        dir: true
    }

    UNTAR(ch_rnafusion_branch.tarzipped)
    ch_versions = ch_versions.mix(UNTAR.out.versions.first())

    def ch_varcov_input = ch_rnafusion_branch.dir
        .mix(UNTAR.out.untar)
        .map { meta, dir ->
            [
                meta,
                getFilesAndCheck(dir, "vcf/*.vcf"),
                getFilesAndCheck(dir, "stringtie/*.gene.abundance.txt"),
                getFilesAndCheck(dir, "fusionreport/*/*.fusions.csv"),
                getFilesAndCheck(dir, "ctatsplicing/*.cancer.introns"),
                getFilesAndCheck(dir, "multiqc/multiqc_data/multiqc_general_stats.txt"),
                getFilesAndCheck(dir, "star/*.Aligned.sortedByCoord.out.bam"),
                getFilesAndCheck(dir, "star/*.Aligned.sortedByCoord.out.bam.bai"),
                meta.run,
            ]
        }

    VARCOV(
        ch_varcov_input,
        ch_genes,
        ch_fusions,
        ch_mane,
        workflow.manifest.version,
    )
    ch_versions = ch_versions.mix(VARCOV.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'report_software_' + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getFilesAndCheck(dir, glob) {
    def files = file("${dir.toUri()}/${glob}")
    if (!files) {
        error("No files found matching pattern '${glob}' in directory '${dir}'")
    }
    return files
}
