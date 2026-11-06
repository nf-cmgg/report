#!/usr/bin/env python

# This script was originally written by Karl Vandepoele
# It has been modified by Nicolas Vannieuwkerke to make Nextflow implementation easier

############################
### new in version 0.3.0 ###
############################

# The merging of the .vcf and fusion report .csv files was changed to a new approach
# removed MANE transcript worksheet



###########################################################

# importing modules
import argparse
import csv
import json
import re
import os
from pathlib import Path
from typing import Any
import pandas as pd
import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill, Border, Side
from openpyxl.utils.cell import get_column_letter
from openpyxl.utils import get_column_letter

###############
## constants ##
###############

SCORE_THRESHOLD: float = 0.18
TOOL_HITS_THRESHOLD: int = 1

CSV_POSITION_COLUMNS = ("arriba", "fusioncatcher", "starfusion")
CSV_CALLER_PREFIXES = {
    "arriba": "ar_",
    "fusioncatcher": "fc_",
    "starfusion": "sf_",
}
KEY_COLUMNS = ("GENEA", "GENEB", "CHRA", "POSA", "CHRB", "POSB")
POSITION_PATTERN = re.compile(r"position:\s*([^,\s]+)#([^,\s]+)")
BREAKPOINT_PATTERN = re.compile(r"([^:]+):(\d+)(?::[+-])?")

########################
## defining functions ##
########################

def normalize_chromosome(chromosome: str) -> str:
    chromosome = chromosome.strip()
    if chromosome.lower().startswith("chr"):
        chromosome = chromosome[3:]
    return chromosome


def parse_breakpoint(value: str) -> tuple[str, str] | None:
    match = BREAKPOINT_PATTERN.fullmatch(value.strip())
    if match is None:
        return None
    return normalize_chromosome(match.group(1)), str(int(match.group(2)))


def parse_csv_positions(value: Any) -> tuple[str, str, str, str] | None:
    if not isinstance(value, str):
        return None
    match = POSITION_PATTERN.search(value)
    if match is None:
        return None
    first = parse_breakpoint(match.group(1))
    second = parse_breakpoint(match.group(2))
    if first is None or second is None:
        return None
    return first[0], first[1], second[0], second[1]


def parse_caller_info(value: Any, prefix: str) -> dict[str, str]:
    if not isinstance(value, str) or not value:
        return {}

    parsed: dict[str, str] = {}
    for item in value.split(","):
        key, separator, item_value = item.partition(":")
        if separator and key.strip():
            parsed[f"{prefix}{key.strip()}"] = item_value.strip()
    return parsed


def read_fusion_csv(csv_path: Path) -> tuple[pd.DataFrame, list[str]]:
    df = pd.read_csv(csv_path, dtype=str, keep_default_na=False)
    if "Fusion" not in df.columns:
        raise ValueError(f"{csv_path} must contain a 'Fusion' column.")

    expanded_rows = []
    for _, row in df.iterrows():
        expanded_row: dict[str, str] = {}
        for column, prefix in CSV_CALLER_PREFIXES.items():
            if column in df.columns:
                expanded_row.update(parse_caller_info(row[column], prefix))
        expanded_rows.append(expanded_row)
    expanded_df = pd.DataFrame(expanded_rows, index=df.index)
    duplicates = set(df.columns).intersection(expanded_df.columns)
    if duplicates:
        raise ValueError(
            "Expanded caller fields would overwrite CSV columns: "
            + ", ".join(sorted(duplicates))
        )
    df = pd.concat([df, expanded_df], axis=1)

    position_columns = [
        column for column in CSV_POSITION_COLUMNS if column in df.columns
    ]
    keys: list[dict[str, Any]] = []
    for row_number, row in df.iterrows():
        genes = row["Fusion"].split("--", maxsplit=1)
        if len(genes) != 2 or not all(genes):
            raise ValueError(
                f"Invalid fusion value on CSV row {row_number + 2}: "
                f"{row['Fusion']!r}"
            )

        positions = {
            parsed
            for column in position_columns
            if (parsed := parse_csv_positions(row[column])) is not None
        }
        if not positions:
            positions.add((None, None, None, None))

        for chra, posa, chrb, posb in sorted(
            positions,
            key=lambda item: tuple("" if part is None else str(part) for part in item),
        ):
            keys.append(
                {
                    "_CSV_ROW": row_number,
                    "GENEA": genes[0],
                    "GENEB": genes[1],
                    "CHRA": chra,
                    "POSA": posa,
                    "CHRB": chrb,
                    "POSB": posb,
                }
            )

    key_df = pd.DataFrame(keys)
    csv_df = df.reset_index(drop=True).copy()
    csv_df["_CSV_ROW"] = csv_df.index
    return key_df.merge(csv_df, on="_CSV_ROW", how="left"), position_columns


def read_fusion_vcf(vcf_path: Path) -> pd.DataFrame:
    header: list[str] | None = None
    records: list[dict[str, Any]] = []
    info_columns: list[str] = []

    with vcf_path.open(newline="") as vcf_file:
        for fields in csv.reader(vcf_file, delimiter="\t"):
            if not fields:
                continue
            if fields[0] == "#CHROM":
                header = [field.lstrip("#") for field in fields]
                continue
            if fields[0].startswith("#"):
                continue
            if header is None:
                raise ValueError(f"{vcf_path} has data before its #CHROM header.")
            if len(fields) != len(header):
                raise ValueError(
                    f"Malformed VCF row {len(records) + 1} in {vcf_path}: "
                    f"expected {len(header)} fields, found {len(fields)}."
                )

            record = dict(zip(header, fields))
            info = {
                item.split("=", maxsplit=1)[0]: item.split("=", maxsplit=1)[1]
                for item in record.get("INFO", "").split(";")
                if "=" in item
            }
            for name in info:
                if name not in info_columns:
                    info_columns.append(name)
            key: dict[str, Any] = {
                "GENEA": info.get("GENEA"),
                "GENEB": info.get("GENEB"),
                "CHRA": info.get("CHRA"),
                "POSA": info.get("POSA"),
                "CHRB": info.get("CHRB"),
                "POSB": info.get("POSB"),
            }
            for chromosome_column in ("CHRA", "CHRB"):
                if key[chromosome_column] not in (None, "", "nan"):
                    key[chromosome_column] = normalize_chromosome(
                        key[chromosome_column]
                    )
            for position_column in ("POSA", "POSB"):
                value = key[position_column]
                key[position_column] = (
                    str(int(value)) if value not in (None, "", "nan") else None
                )
            records.append(
                {
                    **info,
                    **key,
                    **{f"VCF_{name}": value for name, value in record.items()},
                }
            )

    if header is None:
        raise ValueError(f"{vcf_path} is missing its #CHROM header.")
    required_info = {"GENEA", "GENEB", "CHRA", "CHRB", "POSA", "POSB"}
    for index, record in enumerate(records, start=1):
        if any(record[column] is None for column in required_info):
            raise ValueError(
                f"VCF record {index} in {vcf_path} is missing a fusion gene "
                "or breakpoint INFO field."
            )
    return pd.DataFrame(
        records,
        columns=[
            *KEY_COLUMNS,
            *(name for name in info_columns if name not in KEY_COLUMNS),
            *(f"VCF_{name}" for name in header),
        ],
    )


def collapse_duplicate_keys(df: pd.DataFrame) -> pd.DataFrame:
    if not df.duplicated(list(KEY_COLUMNS), keep=False).any():
        return df

    def combine(values: pd.Series) -> Any:
        unique_values: list[str] = []
        for value in values:
            if pd.isna(value):
                continue
            text = str(value)
            if text not in unique_values:
                unique_values.append(text)
        if not unique_values:
            return pd.NA
        if len(unique_values) == 1:
            return unique_values[0]
        return json.dumps(unique_values, ensure_ascii=False)

    value_columns = [column for column in df.columns if column not in KEY_COLUMNS]
    return (
        df.groupby(list(KEY_COLUMNS), dropna=False, sort=False, as_index=False)
        .agg({column: combine for column in value_columns})
    )


def remove_redundant_source_columns(df: pd.DataFrame) -> pd.DataFrame:
    raw_columns = [
        *CSV_CALLER_PREFIXES,
        "VCF_INFO",
    ]
    df = df.drop(columns=raw_columns, errors="ignore")

    score_columns = [
        column for column in df.columns if column.casefold() == "score"
    ]
    if len(score_columns) > 1:
        def non_empty_count(column: str) -> int:
            values = df[column]
            return int(
                (
                    values.notna()
                    & values.astype("string").str.strip().ne("")
                    & values.astype("string").str.casefold().ne("nan")
                ).sum()
            )

        keep_column = max(
            score_columns,
            key=lambda column: (
                non_empty_count(column),
                column == "Score",
            ),
        )
        df = df.drop(
            columns=[column for column in score_columns if column != keep_column]
        )
    return df


def combine_fusion_data(
    csv_path: Path, vcf_path: Path
) -> tuple[pd.DataFrame, list[str]]:
    csv_df, position_columns = read_fusion_csv(csv_path)
    vcf_df = read_fusion_vcf(vcf_path)

    csv_df = collapse_duplicate_keys(csv_df)
    vcf_df = collapse_duplicate_keys(vcf_df)
    combined = csv_df.merge(
        vcf_df,
        on=list(KEY_COLUMNS),
        how="outer",
        indicator=True,
        sort=False,
    )
    combined["Match_status"] = combined.pop("_merge").map(
        {"both": "matched", "left_only": "csv_only", "right_only": "vcf_only"}
    )
    combined = combined.drop(columns=["_CSV_ROW"], errors="ignore")
    combined = remove_redundant_source_columns(combined)
    ordered_columns = [
        *KEY_COLUMNS,
        "Match_status",
        *[
            column
            for column in combined.columns
            if column not in KEY_COLUMNS + ("Match_status",)
        ],
    ]
    return combined[ordered_columns], position_columns


# function to create dir if it does not exist yet
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")

# Adds chr to chromosomes missing it
def update_chroms(chrom:any) -> str:
    if pd.isna(chrom): return chrom
    str_chrom = str(chrom)
    return str_chrom if str_chrom.startswith("chr") else f"chr{str_chrom}"

# create a button in the summary sheet to link to the relevant worksheet
def create_cell_button(ws, cell_ref, text, link=None,
                        fill='00E8F4F8', font_color='000000',
                        border_color='00B3E5FC', border_style='medium', bold=True):
    c = ws[cell_ref]
    c.value = text
    c.font = Font(color=font_color, bold=bold)
    c.alignment = Alignment(horizontal='center', vertical='center')
    c.fill = PatternFill('solid', fgColor=fill)

    # Medium border for “button” frame
    side = Side(style=border_style, color=border_color)
    c.border = Border(left=side, right=side, top=side, bottom=side)

    # Optional hyperlink
    if link:
        c.hyperlink = link

    # Tidy sizing for a button feel
    ws.row_dimensions[c.row].height = 20
    col_letter = get_column_letter(c.column)
    ws.column_dimensions[col_letter].width = max(len(text) + 6, 18)

# read in the arriba output file to correct if necessary
def read_arriba_file(arriba_file:str) -> pd.DataFrame:
    df = pd.read_csv(arriba_file, sep='\t')
    df['Fusion'] = df['#gene1'] + "--" + df['gene2']
    df[['str1gene', 'str1fusion']] = df['strand1(gene/fusion)'].str.split('/', n=1,expand=True)
    df[['str2gene', 'str2fusion']] = df['strand2(gene/fusion)'].str.split('/', n=1,expand=True)
    df['ORIENTATION'] = df['str1fusion'] + '/' + df['str2fusion']
    df['TRANSCRIPT_A'] = df['transcript_id1']
    df['TRANSCRIPT_B'] = df['transcript_id2']
    df['FRAME_STATUS'] = df['reading_frame']
    df[['CHRA', 'POSA']] = df['breakpoint1'].str.split(':', n=1, expand=True)
    df[['CHRB', 'POSB']] = df['breakpoint2'].str.split(':', n=1, expand=True)
    df['CHRA'] = df['CHRA'].apply(update_chroms)
    df['CHRB'] = df['CHRB'].apply(update_chroms)
    df['POSA'] = pd.to_numeric(df['POSA'])
    df['POSB'] = pd.to_numeric(df['POSB'])
    df = df[['Fusion', 'CHRA', 'CHRB', 'POSA', 'POSB', 'ORIENTATION', 'TRANSCRIPT_A', 'TRANSCRIPT_B', 'FRAME_STATUS']]
    return df

#############################################
## read in arguments and store as variable ##
#############################################

parser = argparse.ArgumentParser(description="Create a report from output files of the nf-core/rnafusion pipeline")
parser.add_argument('--input', metavar='DIR', type=str, help="Path to a directory containing the input VCFs", required=True)
parser.add_argument('--stringtie', metavar='DIR', type=str, help="Path to a directory containing StringTie gene abundance output files", required=True)
parser.add_argument('--fusionreport', metavar='DIR', type=str, help="Path to a directory containing fusion report files with detected fusion and the amount of reads analysed by each fusion caller", required=True)
parser.add_argument('--ctat', metavar='DIR', type=str, help="Path to a directory containing CTAT-SPLICING cancer intron files", required=True)
parser.add_argument('--multiqc', metavar='DIR', type=str, help="Path to a directory containing the multiQC general stats file", required=True)
parser.add_argument('--output', metavar='DIR', type=str, help="Path to the output directory", required=True)
parser.add_argument('--bams', metavar='DIR', type=str, help="Path to a directory containing BAM files produced by STAR", required=True)
parser.add_argument('--genes', metavar='FILE', type=str, help="Path to a file containing the targeted genes list for the analysis.", required=True)
parser.add_argument('--fusion_whitelist', metavar='FILE', type=str, help="Path to a file containing a fusion whitelist for specific filtering.", required=True)
# parser.add_argument('--mane', metavar='FILE', type=str, help="Path to a file containing the MANE transcripts for filtering.", required=True)
parser.add_argument('--run', metavar='STR', type=str, help="Name of the run analysed", required=True)
parser.add_argument('--pipeline_version', metavar='STR', type=str, help="Version of the reporting pipeline used", required=True)
parser.add_argument('--arriba', metavar='STR', type=str, help="Path to arriba output directory", required=True)
parser.add_argument('--design', metavar='STR', type=str, help="TWIST design used in the analysis")

args = parser.parse_args()
input_path:str = args.input + "/"
cov_path:str = args.stringtie + "/"
reads_path:str = args.fusionreport + "/"
splicing_path:str = args.ctat + "/"
qc_path:str = args.multiqc + "/"
output_dir:str = args.output
bam_path:str = args.bams + "/"
genes:str = pd.read_csv(args.genes, sep='\t')
fusions:str = pd.read_csv(args.fusion_whitelist, sep='\t')
# mane:str = pd.read_csv(args.mane)
run_nr:str = args.run
pipeline_version:str = args.pipeline_version
arriba_path:str = args.arriba + "/"
design:str = args.design

###############################################################
### Loop over .vcf files, extract info and export to report ###
###############################################################

# Loop over all files in the folder
for filename in os.listdir(input_path):
    # Check if the file is a VCF file
    if filename.endswith(".vcf"):
        # Full path to the VCF file
        vcf_path = os.path.join(input_path, filename)
        basename = filename.split('_')[0]
        reads_file = os.path.join(reads_path + basename + ".fusions.csv")
        df, _ = combine_fusion_data(Path(reads_file), Path(vcf_path))

        # Get sample names
        columns = df.columns.tolist()
        format_column_index = columns.index("VCF_FORMAT")
        samples: list[str] = [
            column.removeprefix("VCF_")
            for column in columns[format_column_index + 1:]
            if column.startswith("VCF_")
        ]
        if len(samples) != 1:
            raise ValueError(f"Expected exactly one sample in VCF file {filename}, but found {len(samples)} samples.")
        sample_name: str = samples[0]

        print(f"Creating report for sample: {sample_name}")

        df_clean = df
        df_clean["CHRA"] = df_clean["CHRA"].apply(update_chroms)
        df_clean["CHRB"] = df_clean["CHRB"].apply(update_chroms)
        df_clean["Fusion"] = df_clean["GENEA"] + "--" + df_clean["GENEB"]

        # change annotation of reference sequence to one column
        # output is seen as float, not string (is sometimes empty so pandas consideres it as float)
        df_clean['TRANSCRIPT_A'] = df_clean['TRANSCRIPT_ID_A'].astype(str)
        df_clean = df_clean.drop('TRANSCRIPT_VERSION_A', axis = 1)
        df_clean = df_clean.drop('TRANSCRIPT_ID_A', axis = 1)

        df_clean['TRANSCRIPT_B'] = df_clean['TRANSCRIPT_ID_B'].astype(str)
        df_clean = df_clean.drop('TRANSCRIPT_VERSION_B', axis = 1)
        df_clean = df_clean.drop('TRANSCRIPT_ID_B', axis = 1)

        merged_df = df_clean
        merged_df = merged_df.drop(
            columns=["ar_position", "fc_position", "sf_position"],
            errors="ignore",
        )

        # set certain columns to numerical
        present_columns_numerical: list[str] = list(merged_df.columns.intersection(["POSA", "POSB", 'fc_longest_anchor',
            'fc_common_mapping_reads', 'fc_spanning_pairs','fc_spanning_unique_reads', 'ar_coverage1',
            'ar_coverage2', 'ar_discordant_mates', 'ar_split_reads1', 'ar_split_reads2', 'sf_ffmp',
            'sf_junction_reads', 'sf_spanning_reads', 'SCORE', 'TOOL_HITS']))
        for column in present_columns_numerical:
            merged_df[column] = pd.to_numeric(merged_df[column])

        # limit sf_ffmp to one decimal
        if 'sf_ffmp' in merged_df.columns:
            merged_df['sf_ffmp'] = round(merged_df['sf_ffmp'], 1)

        # add missing annotation data
        ## arriba check
        arriba_nan_df: pd.DataFrame = merged_df[merged_df['FOUND_IN'].str.contains('arriba')]
        arriba_nan_df = arriba_nan_df[(arriba_nan_df["TRANSCRIPT_A"] == 'nan') | (arriba_nan_df["TRANSCRIPT_B"] == 'nan')]
        if not arriba_nan_df.empty:
            arriba_df = read_arriba_file(os.path.join(arriba_path + basename + ".arriba.fusions.tsv"))
            arriba_indexed = arriba_df.set_index(['Fusion', 'CHRA', 'CHRB', 'POSA', 'POSB'])
            merged_indexed = merged_df.set_index(['Fusion', 'CHRA', 'CHRB', 'POSA', 'POSB'])
            merged_df = arriba_indexed.combine_first(merged_indexed).reset_index()

        # subset the dataframe to only contain relevant columns
        relevant_columns: list[str] = [
            "Fusion",
            "CHRA",
            "CHRB",
            "GENEA",
            "GENEB",
            "POSA",
            "POSB",
            "ORIENTATION",
            "FOUND_IN",
            "TOOL_HITS",
            "SCORE",
            "FRAME_STATUS",
            "EXON_NUMBER_A",
            "EXON_NUMBER_B",
            "TRANSCRIPT_A",
            "TRANSCRIPT_B"
        ]

        def starts_with_list(string:str, prefixes:list[str]) -> bool:
            return any(string.startswith(prefix) for prefix in prefixes)

        tool_specific_columns: list[str] = [column for column in merged_df.columns if starts_with_list(column, ['fc_', 'ar_', 'sf_'])]

        relevant_columns.extend(tool_specific_columns)

        df_filt = merged_df[relevant_columns]

        # sort on score
        df_filt = df_filt.sort_values('SCORE', ascending = False)

        # move some columns
        column_to_move = df_filt.pop('EXON_NUMBER_A')
        df_filt.insert(5, 'EXONA', column_to_move)

        column_to_move = df_filt.pop('EXON_NUMBER_B')
        df_filt.insert(6, 'EXONB', column_to_move)

        column_to_move = df_filt.pop('Fusion')
        df_filt.insert(0, 'Fusion', column_to_move)

        df_filt.dropna(subset=tool_specific_columns, how='all', inplace=True)

        # filter base on score and number of callers, both conditions have to be true
        df_final = df_filt[(df_filt['SCORE'] > SCORE_THRESHOLD) & (df_filt['TOOL_HITS'] > TOOL_HITS_THRESHOLD)]

        # limit filtered to versioned MANE transcript only
        # df_MANE = mane['Transcript stable ID version'].astype(str)
        # df_final_MANE = df_final[df_final['TRANSCRIPT_A'].apply(lambda x: any(i in x for i in df_MANE)) & df_final['TRANSCRIPT_B'].apply(lambda x: any(i in x for i in df_MANE))]


        #########################
        ### Coverage analysis ###
        #########################

        # read coverage into a new df
        cov_file = os.path.join(cov_path + basename + ".gene.abundance.txt")
        expression = pd.read_csv(cov_file, sep = '\t', dtype={2: 'str'})

        # Round the columns with numerical values to one digit
        columns_to_round = ['Coverage', 'FPKM', 'TPM']
        expression[columns_to_round] = expression[columns_to_round].round(1)

        # filter expression values to contain only targeted genes
        genes_in_df = genes['genes']
        genes_expr = expression[expression['Gene Name'].isin(genes_in_df)]

        # rename column and sort based on gene name
        genes_expr = genes_expr.rename(columns={'Reference': 'Chromosome'})
        genes_expr = genes_expr.sort_values('Gene Name')

        # make separate df for control genes
        ref_genes = ['CHMP2A', 'GPI', 'RAB7A', 'VCP']
        ref_genes_expr = genes_expr[genes_expr['Gene Name'].isin(ref_genes)]

        # make df for whitelisted fusions  (should be placed more upwards in fusion section)
        fusions_in_df = fusions['gene A']
        df_fusions = df_filt[df_filt['GENEA'].isin(fusions_in_df) | df_filt['GENEB'].isin(fusions_in_df)]

        # calculate coverage for DUX4 (mapping quality 0)
        bashCommand = "samtools coverage -r chr4:190173000-190175500 " + bam_path + basename + "*.bam" + " -q 0 -A -H > " + basename + "_DUX4.txt"
        os.system(bashCommand)

        DUX4_file = basename + "_DUX4.txt"

        DUX4_reads = "0"
        if(os.path.exists(DUX4_file)):
            with open(DUX4_file, 'r') as file:
                data = file.read()
                # alternatively look for 'Number of reads' and print these (see separate script for this code)
                result = re.search(r'\((.*?) filtered\)', data)
                if result:
                    DUX4_reads = result.group(1)

        #####################
        ### Splicing data ###
        #####################

        # read in splicing parameters
        splicing_file = os.path.join(splicing_path + basename + ".cancer.introns")
        splicing = pd.read_csv(splicing_file, sep = '\t', dtype={2: 'str'})

        # move variant_name column up front for easy identification
        column_to_move = splicing.pop('variant_name')
        splicing.insert(0, 'variant_name', column_to_move)

        ctat_variants_found: bool = splicing['variant_name'].count()

        ###############
        ### QC data ###
        ###############

        # read in QC parameters
        qc_file = os.path.join(qc_path + "multiqc_general_stats.txt")
        qc = pd.read_csv(qc_file, sep = '\t', dtype={2: 'str'})

        # filter so only the data from this sample are retained
        qc_filt = qc[qc['Sample'].str.contains(basename)]

        def clean_column_name(col_name: str) -> str:
            col_name = col_name.strip().replace('-', ' ').replace('_', ' ')
            return col_name.title()

        qc_filt.columns = [clean_column_name(col) for col in qc_filt.columns]
        # Set Sample as a temporary index to convert all other columns to a float
        qc_filt = qc_filt.set_index(['Sample'], append=True).astype(float).round(2).reset_index(['Sample'])

        ####################################
        ### Save all df to an Excel file ###
        ####################################

        excel_path: str = f"{output_dir}/{sample_name}.xlsx"

        # Sheet names
        fusions_filt_sheet = "fusions_filt"
        #fusions_filt_MANE_sheet = "fusions_filt_MANE"
        fusions_specific_sheet = "fusions_specific"
        fusions_all_sheet = "fusions_all"
        splicing_sheet = "CTAT_splicing"
        coverage_ref_sheet = "coverage_ref"
        coverage_all_sheet = "coverage_all"
        qc_sheet = "QC"
        summary_sheet = "Summary"

        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            df_final.to_excel(writer, index=False, sheet_name=fusions_filt_sheet)
            #df_final_MANE.to_excel(writer, index = False, sheet_name = fusions_filt_MANE_sheet)
            df_fusions.to_excel(writer, index= False, sheet_name=fusions_specific_sheet)
            merged_df.to_excel(writer, index=False, sheet_name=fusions_all_sheet)
            splicing.to_excel(writer, index=False, sheet_name = splicing_sheet)
            ref_genes_expr.to_excel(writer, index= False, sheet_name=coverage_ref_sheet)
            genes_expr.to_excel(writer, index=False, sheet_name=coverage_all_sheet)
            qc_filt.to_excel(writer, index=False, sheet_name=qc_sheet)


        ###############################
        ### Change layout of report ###
        ###############################

        # open report file
        workbook = openpyxl.load_workbook(excel_path)

        # create a summary sheet
        workbook.create_sheet(title=summary_sheet, index=0)
        ws = workbook[summary_sheet]
        ws['A1'] = f"Summary for {basename}"
        ws['A1'].font = Font(size = 14, bold = True)
        ws['A1'].fill = PatternFill("solid", fgColor='00FFCC99')

        ws['A3'] = "Filtered fusions found:"
        ws['B3'] = df_final.shape[0]
        create_cell_button(ws, 'C3', 'Go to data', link=f"#{fusions_filt_sheet}!A1")

        # ws['A4'] = "Filtered fusion present in MANE transcripts:"
        # ws['B4'] = df_final_MANE.shape[0]
        # create_cell_button(ws, 'C4', 'Go to data', link=f"#{fusions_filt_MANE_sheet}!A1")

        ws['A4'] = "Fusions in whitelist found:"
        ws['B4'] = df_fusions.shape[0]
        create_cell_button(ws, 'C4', 'Go to data', link=f"#{fusions_specific_sheet}!A1")

        ws['A5'] = "Total fusions found:"
        ws['B5'] = merged_df.shape[0]
        create_cell_button(ws, 'C5', 'Go to data', link=f"#{fusions_all_sheet}!A1")

        ws['A6'] = 'Splice sites found:'
        ws['B6'] = splicing.shape[0]
        create_cell_button(ws, 'C7', 'Go to data', link=f"#{splicing_sheet}!A1")
        ws['D6'] = f' {ctat_variants_found} CTAT variant{"s" if ctat_variants_found != 1 else ""} found'
        if ctat_variants_found > 0:
            ws['D6'].font = Font(color='00FF0000', bold=True)

        ws['A8'] = "Coverage analysis of control genes:"
        create_cell_button(ws, 'B8', 'Go to data', link=f"#{coverage_ref_sheet}!A1")

        ws['A9'] = "Coverage analysis of all targeted genes:"
        create_cell_button(ws, 'B9', 'Go to data', link=f"#{coverage_all_sheet}!A1")

        ws['A10'] = "QC metrics:"
        create_cell_button(ws, 'B10', 'Go to data', link=f"#{qc_sheet}!A1")

        ws['A12'] = "Patient:"
        ws['A13'] = "RNA-nr:"
        ws['B13'] = basename
        ws['A14'] = "Date of birth:"
        ws['A15'] = "TWIST design:"
        ws['B15'] = design if design != None else "?"
        ws['A16'] = "Sequencing run:"
        ws['B16'] = run_nr

        ws.column_dimensions["A"].width = 40
        ws.column_dimensions["B"].width = len(design) + 5 if design != None and len(design) > 20 else 20
        ws.column_dimensions["D"].width = 25

        # apply styling to columns in summary
        for idx, cell in enumerate(ws['A']):
            if idx == 0: continue  # skip header
            cell.alignment = Alignment(horizontal='right', vertical='center')
            cell.font = Font(bold=True)

        for cell in ws['B']:
            cell.alignment = Alignment(horizontal='center', vertical='center')

        for cell in ws['D']:
            cell.alignment = Alignment(horizontal='left', vertical='center')

        # modify fusion worksheets
        for worksheet in [fusions_all_sheet, fusions_filt_sheet, fusions_specific_sheet]:
            ws = workbook[worksheet]

            # insert sample name in new empty row and add hyperlink to visualisation files
            ws.insert_rows(1)
            ws['A1'] = "fusions " + basename

            # change layout
            c = ws['A1']
            c.font = Font(size = 14, bold = True)
            c.fill = PatternFill("solid", fgColor='00FFCC99')

            # change column width
            for column_index in range(1, 38):
                column_letter = get_column_letter(column_index)
                ws.column_dimensions[column_letter].width = 12

            for column_index in range(15, 21):
                column_letter = get_column_letter(column_index)
                ws.column_dimensions[column_letter].width = 25

            ws.column_dimensions["A"].width = 20
            ws.column_dimensions["K"].width = 32
            ws.column_dimensions["N"].width = 20

            # change alignment/cell height for header
            for cell in ws["2:2"]:
                cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)

            ws.row_dimensions[2].height = 50

            '''
            # add ensembl link to transcript columns
            def link_transcript(column:str) -> None:
                column_letter: str = openpyxl.utils.cell.get_column_letter(column)
                for cell in ws[f'{column_letter}3:{column_letter}{ws.max_row}']:
                    cell0 = cell[0]
                    if cell0.value not in ["", "nan"]:
                        cell0.hyperlink = f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t={cell0.value}"

            for cell in ws['2']:
                if cell.value == 'TRANSCRIPT_A':
                    link_transcript(cell.column)
                if cell.value == 'TRANSCRIPT_B':
                    link_transcript(cell.column)
            '''

        #loop over coverage worksheets to change layout
        for worksheet in [coverage_ref_sheet, coverage_all_sheet]:
            ws = workbook[worksheet]

            # add new line and add sample name
            ws.insert_rows(1)
            ws['A1']= "coverage " + basename

            # change layout
            c = ws['A1']
            c.font = Font(size = 14, bold = True)
            c.fill = PatternFill("solid", fgColor='00FFCC99')

            for column_index in range(1, 10):
                column_letter = get_column_letter(column_index)
                ws.column_dimensions[column_letter].width = 14

            ws.column_dimensions["A"].width = 20

        # add DUX4 data to coverage file
        ws = workbook[coverage_all_sheet]

        ws['J26'] = DUX4_reads + " filtered reads"

        # change layout of splicing worksheet
        for worksheet in [splicing_sheet]:
            ws = workbook[worksheet]

            # add new line and add sample name
            ws.insert_rows(1)
            ws['A1']= "CTAT splicing " + basename

            # change layout
            c = ws['A1']
            c.font = Font(size = 14, bold = True)
            c.fill = PatternFill("solid", fgColor='00FFCC99')

            ws.column_dimensions["A"].width = 20
            ws.column_dimensions["B"].width = 30
            ws.column_dimensions["D"].width = 30
            ws.column_dimensions["E"].width = 20
            ws.column_dimensions["F"].width = 20
            ws.column_dimensions["G"].width = 40
            ws.column_dimensions["H"].width = 40

        # change layout of QC worksheet
        for worksheet in [qc_sheet]:
            ws = workbook[worksheet]

            # add new line and add sample name
            ws.insert_rows(1)
            ws['A1']= "QC " + basename

            # add script version for logging
            ws['E1'] = "nf-cmgg/report version: " + pipeline_version

            # change layout
            c = ws['A1']
            c.font = Font(size = 14, bold = True)
            c.fill = PatternFill("solid", fgColor='00FFCC99')

            ws.column_dimensions["A"].width = 25

            for column_index in range(2, 28):
                column_letter = get_column_letter(column_index)
                ws.column_dimensions[column_letter].width = 15

            for cell in ws["2:2"]:
                cell.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)

            ws.row_dimensions[2].height = 100

        # add back to summary button on all sheets except for summary
        for sheet_name in workbook.sheetnames:
            if sheet_name != summary_sheet:
                ws = workbook[sheet_name]
                create_cell_button(ws, 'C1', 'Back to summary', link=f"#{summary_sheet}!A1", border_style='thin', border_color='000000')

        # save file
        workbook.save(excel_path)
