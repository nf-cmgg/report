#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/report
----------------------------------------------------------------------------------------
*/

nextflow.preview.output = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_report_pipeline'

params.targeted.fasta       = params.targeted.fasta ?: getGenomeAttribute('fasta')
params.targeted.queries_dir = "${projectDir}/assets/queries/"

include { TARGETED                } from './workflows/targeted'
include { RNAFUSION               } from './workflows/rnafusion'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_report_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_report_pipeline'
include { softwareVersionsToYAML  } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { samplesheetToList       } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION(
        params.version,
        params.validate_params,
        args,
        params.outdir
    )

    //
    // Run reporting flows
    //

    def ch_versions = Channel.empty()

    def out_targeted_hotcount = Channel.empty()
    if(params.targeted.input) {
        def required_parameters = ['fasta', 'queries_dir']
        check_required_params('targeted', required_parameters)
        def targeted_params = params.targeted

        def ch_samplesheet = Channel.fromList(samplesheetToList(targeted_params.input, "${projectDir}/assets/schema_targeted_input.json"))
        def fasta = Channel.value([[id:'reference'], file(targeted_params.fasta)])

        TARGETED(
            ch_samplesheet,
            fasta,
            targeted_params.queries_dir
        )
        out_targeted_hotcount = TARGETED.out.hotcount
        ch_versions = ch_versions.mix(TARGETED.out.versions)
    }

    def out_rnafusion_excels = Channel.empty()
    if(params.rnafusion.input) {
        def required_parameters = ['genes', 'fusions', 'mane']
        check_required_params('rnafusion', required_parameters)
        def rnafusion_params = params.rnafusion

        def ch_samplesheet = Channel.fromList(samplesheetToList(rnafusion_params.input, "${projectDir}/assets/schema_rnafusion_input.json"))
        def genes = Channel.value(file(rnafusion_params.genes))
        def fusions = Channel.value(file(rnafusion_params.fusions))
        def mane = Channel.value(file(rnafusion_params.mane))

        RNAFUSION(
            ch_samplesheet,
            genes,
            fusions,
            mane,
            workflow.manifest.version
        )
        out_rnafusion_excels = RNAFUSION.out.excels
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'report_software_' + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION(
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )

    publish:
    targeted_hotcount = out_targeted_hotcount
    rnafusion_excels  = out_rnafusion_excels
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    targeted_hotcount {
        path { meta, counts ->
            counts >> "targeted/hotcount/${meta.id}.counts.txt"
        }
    }
    rnafusion_excels {
        path { meta, excel ->
            excel >> "rnafusion/varcov/${meta.run}/"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def check_required_params(scope_name, required_params) {
    def scope_params = params.get(scope_name, null)
    if (scope_params == null) {
        error "Could not find a parameters scope called '${scope_name}'"
    }

    def missing_params = []
    required_params.each { param ->
        if (scope_params.get(param, null) == null) {
            missing_params << param
        }
    }

    if (missing_params) {
        error "The following required parameters are missing from the '${scope_name}' parameters scope: ${missing_params.join(', ')}"
    }
}
