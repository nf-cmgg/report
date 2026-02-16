#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-cmgg/report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-cmgg/report
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_report_pipeline'

params.targeted.fasta       = params.targeted.fasta ?: getGenomeAttribute('fasta')
params.targeted.gene        = params.targeted.gene
params.targeted.queries_dir = "${projectDir}/assets/targeted/"

include { TARGETED                } from './workflows/targeted'
include { RNAFUSION               } from './workflows/rnafusion'
include { PACVAR_REPEAT           } from './modules/local/pacvarrepeat'
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

    def out_targeted_hotcount = channel.empty()
    if(params.targeted.input) {
        def required_parameters = ['fasta', 'queries_dir', 'gene']
        check_required_params(params.get('targeted'), 'targeted', required_parameters)
        def targeted_params = params.targeted

        def ch_samplesheet = channel.fromList(samplesheetToList(targeted_params.input, "${projectDir}/assets/schema_targeted_input.json"))
        def fasta = channel.value([[id:'reference'], file(targeted_params.fasta)])

        TARGETED(
            ch_samplesheet,
            fasta,
            targeted_params.queries_dir,
            targeted_params.gene
        )
        out_targeted_hotcount = TARGETED.out.hotcount
    }

    def out_rnafusion_excels = channel.empty()
    if(params.rnafusion.input) {
        def required_parameters = ['genes', 'fusions', 'mane']
        check_required_params(params.get('rnafusion'), 'rnafusion', required_parameters)
        def rnafusion_params = params.rnafusion

        def ch_samplesheet = channel.fromList(samplesheetToList(rnafusion_params.input, "${projectDir}/assets/schema_rnafusion_input.json"))
        def genes = channel.value(file(rnafusion_params.genes))
        def fusions = channel.value(file(rnafusion_params.fusions))
        def mane = channel.value(file(rnafusion_params.mane))

        RNAFUSION(
            ch_samplesheet,
            genes,
            fusions,
            mane,
            workflow.manifest.version
        )
        out_rnafusion_excels = RNAFUSION.out.excels
    }

    def out_pacvar_repeat_excels = channel.empty()
    if(params.pacvar_repeat.input) {
        def required_parameters = ['input']
        check_required_params(params.get('pacvar_repeat'), 'pacvar_repeat', required_parameters)
        def pacvar_repeat_params = params.pacvar_repeat
        def ch_samplesheet = channel.fromList(samplesheetToList(file(pacvar_repeat_params.input), "${projectDir}/assets/schema_pacvar_repeat_input.json"))

        PACVAR_REPEAT(ch_samplesheet)
        out_pacvar_repeat_excels = PACVAR_REPEAT.out.excels
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
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
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
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
    pacvar_repeat_excels  = out_pacvar_repeat_excels
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    OUTPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

output {
    targeted_hotcount {
        path { meta, counts ->
            counts >> "targeted/${params.targeted.gene}/${meta.id}.counts.txt"
        }
    }
    rnafusion_excels {
        path { meta, excel ->
            excel >> "rnafusion/varcov/${meta.run}/"
        }
    }
    pacvar_repeat_excels {
        path { _meta, excel ->
            excel >> "pacvar_repeat/reports/"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def check_required_params(Map scope_params, String scope_name, List<String> required_params) {
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
