include { SMALLVARIANTS_TO_EXCEL } from '../modules/local/smallvariantstoexcel/main.nf'

workflow SEQCAP_SMALLVARIANTS {
    take:
    samplesheet        // value: path to the run samplesheet (validated with schema_seqcap_smallvariants_input.json)
    ch_rows            // channel: [meta, panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage]
    api_data           // value: path to the combined api_data.json bundle
    run_name           // value: string used as the run id
    threshold_coverage // value: integer coverage threshold for low coverage flagging
    build              // value: string genome build label
    variant_caller     // value: string variant caller name
    runtype            // value: string runtype label

    main:

    // The python script processes the whole run (a single samplesheet) at once, so it needs
    // every file referenced by the samplesheet staged in, but the samplesheet itself is used
    // unmodified (it already contains valid, schema-checked absolute paths).
    def ch_input_files = ch_rows
        .flatMap { meta, panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage ->
            [panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage].findAll { file -> file }
        }
        .unique { file -> file.name }
        .collect()

    SMALLVARIANTS_TO_EXCEL(
        channel.value([[id: run_name], samplesheet]),
        api_data,
        ch_input_files,
        threshold_coverage,
        build,
        variant_caller,
        runtype,
    )

    emit:
    excels  = SMALLVARIANTS_TO_EXCEL.out.excels
    reports = SMALLVARIANTS_TO_EXCEL.out.reports
}
