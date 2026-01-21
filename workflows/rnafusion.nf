/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNTAR     } from '../modules/nf-core/untar/main.nf'
include { VARCOV    } from '../modules/local/varcov/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNAFUSION {
    take:
    ch_samplesheet // channel: [ meta, path(directory) ]
    ch_genes
    ch_fusions
    ch_mane
    workflow_version

    main:

    def ch_versions = Channel.empty()

    def ch_input_branch = ch_samplesheet.branch { _meta, dir ->
        tarzipped: dir.extension == "gz"
        dir: true
    }

    UNTAR(ch_input_branch.tarzipped)
    ch_versions = ch_versions.mix(UNTAR.out.versions.first())

    def ch_varcov_input = ch_input_branch.dir
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
                getFilesAndCheck(dir, "arriba/*.fusions.tsv"),
                meta.run,
                meta.design
            ]
        }

    VARCOV(
        ch_varcov_input,
        ch_genes,
        ch_fusions,
        ch_mane,
        workflow_version,
    )

    emit:
    excels   = VARCOV.out.output // channel: [ val(meta), path(excel) ]
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
