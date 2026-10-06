include { SMALLVARIANTS_TO_EXCEL } from '../modules/local/smallvariantstoexcel/main.nf'

workflow SEQCAP_SMALLVARIANTS {
    take:
    ch_rows            // channel: [meta, panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage] (parsed from the schema-validated input samplesheet)
    api_data           // value: path to the combined api_data.json bundle
    run_name           // value: string used as the run id
    threshold_coverage // value: integer coverage threshold for low coverage flagging
    build              // value: string genome build label
    variant_caller     // value: string variant caller name
    runtype            // value: string runtype label

    main:

    // The python script processes the whole run at once. The input samplesheet can be in any
    // format supported by samplesheetToList (csv, tsv, yaml, json), so a standardized csv is
    // generated from the parsed rows. It references the staged input files by file name only,
    // since remote files (e.g. URLs) cannot be opened by the script directly.
    def ch_input_files = ch_rows
        .flatMap { meta, panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage ->
            [panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage].findAll { file -> file }
        }
        .unique { file -> file.name }
        .collect()

    def ch_samplesheet = ch_rows
        .collectFile(
            name: 'samplesheet_smallvariants_to_excel.csv',
            seed: 'sample,panel,design,panel_bed,design_bed,panel_genelist,design_genelist,transcript_file,vcf,coverage',
            sort: true,
            newLine: true,
        ) { meta, panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage ->
            [meta.id, meta.panel, meta.design, panel_bed, design_bed, panel_genelist, design_genelist, transcript_file, vcf, coverage]
                .collect { value -> value ? (value instanceof Path ? value.name : value.toString()) : '' }
                .join(',')
        }
        .map { samplesheet -> [[id: run_name], samplesheet] }
        .first()

    SMALLVARIANTS_TO_EXCEL(
        ch_samplesheet,
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
