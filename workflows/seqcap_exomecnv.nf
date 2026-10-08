include { EXOMECNV_TO_CSV } from '../modules/local/exomecnvtocsv/main.nf'

workflow SEQCAP_EXOMECNV {
    take:
    ch_rows   // channel: [meta, design_genelist, panel_genelist, vcf, cnv_database] (parsed from the schema-validated input samplesheet)
    run_name  // value: string used as the run id

    main:

    // The python script processes the whole run at once. The input samplesheet can be in any
    // format supported by samplesheetToList (csv, tsv, yaml, json), so a standardized csv is
    // generated from the parsed rows. It references the staged input files by file name only,
    // since remote files (e.g. URLs) cannot be opened by the script directly.
    def ch_input_files = ch_rows
        .flatMap { meta, design_genelist, panel_genelist, vcf, cnv_database ->
            [design_genelist, panel_genelist, vcf, cnv_database].findAll { file -> file }
        }
        .unique { file -> file.name }
        .collect()

    def ch_samplesheet = ch_rows
        .collectFile(
            name: 'samplesheet_exomecnv_to_csv.csv',
            seed: 'sample,design,panel,design_genelist,panel_genelist,vcf,cnv_database',
            sort: true,
            newLine: true,
        ) { meta, design_genelist, panel_genelist, vcf, cnv_database ->
            [meta.id, meta.design, meta.panel, design_genelist, panel_genelist, vcf, cnv_database]
                .collect { value -> value ? (value instanceof Path ? value.name : value.toString()) : '' }
                .join(',')
        }
        .map { samplesheet -> [[id: run_name], samplesheet] }

    EXOMECNV_TO_CSV(
        ch_samplesheet,
        ch_input_files,
    )

    emit:
    summaries = EXOMECNV_TO_CSV.out.summaries
}
