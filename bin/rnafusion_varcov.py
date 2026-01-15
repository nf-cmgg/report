#!/usr/bin/env python

# This script was originally written by Karl Vandepoele
# It has been modified by Nicolas Vannieuwkerke to make Nextflow implementation easier

############################
### new in version 0.2.8 ###
############################

# deal with empty columns for arriba and starfusion

# TO DO: make script compatible with .gz .vcf files
# always check first if file exists and if not, display error message in the worksheet instead of crashing the script



###########################################################

# importing modules
import argparse
import re
import os
import pandas as pd
import openpyxl
import cyvcf2
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils.cell import get_column_letter

###############
## constants ##
###############

VCF_INFO_COLUMNS = [
    'CHRA',
    'CHRB',
    'GENEA',
    'GENEB',
    'POSA',
    'POSB',
    'SCORE',
    'TOOL_HITS',
    'SCORE',
    'FRAME_STATUS',
    'TRANSCRIPT_ID_A',
    'TRANSCRIPT_ID_B',
    'EXON_NUMBER_A',
    'EXON_NUMBER_B',
    'TRANSCRIPT_VERSION_A',
    'TRANSCRIPT_VERSION_B',
    'HGNC_ID_A',
    'HGNC_ID_B',
    'ANNOTATIONS',
    'ORIENTATION',
    'FOUND_DB',
    'FOUND_IN',
]

VCF_COLUMNS = [
    "CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "FORMAT"
]

########################
## defining functions ##
########################

# function to create dir if it does not exist yet
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")

# Function to fetch needed VCF fields and convert them to a DataFrame
def vcf_to_df(vcf):
    processed_variants: list[dict[str,any]] = []
    for variant in vcf:
        var_dict: dict[str,any] = {}
        var_dict["CHROM"] = update_chroms(variant.CHROM)
        var_dict["POS"] = variant.POS
        var_dict["ID"] = variant.ID
        var_dict["REF"] = variant.REF
        var_dict["ALT"] = update_alt_chroms(variant.ALT)
        var_dict["QUAL"] = variant.QUAL
        var_dict["FILTER"] = variant.FILTER
        var_dict["FORMAT"] = str(variant.FORMAT)
        for column in VCF_INFO_COLUMNS:
            var_dict[column] = variant.INFO.get(column, None)
        processed_variants.append(var_dict)
    df = pd.DataFrame(processed_variants)
    if len(df.index) == 0:
        # Create an empty DataFrame with the appropriate headers if the VCF file is empty
        return pd.DataFrame(columns=VCF_INFO_COLUMNS)
    return df

# Adds chr to chromosomes missing it
def update_chroms(chrom:str) -> str:
    str_chrom = str(chrom)
    return str_chrom if str_chrom.startswith("chr") else f"chr{str_chrom}"

# Adds chr to alt allele chromosomes missing it
def update_alt_chroms(alts: list[str]) -> (list[str]|str):
    outputs = []
    for alt in alts:
        if "chr" not in alt:
            chr_match = re.match(r"^.*[\[\]]([\dXY]{1,2}):.*$", alt)
            if chr_match:
                chrom:str = chr_match.group(1)
                alt = alt.replace(f"{chrom}:", f"chr{chrom}:")
        outputs.append(alt)
    return outputs[0] if len(outputs) == 1 else outputs

# Function to remove brackets and single quotes
def remove_brackets(df, cols):
    for col in cols:
        df[col] = df[col].apply(lambda x: x[0] if isinstance(x, list) and len(x) == 1 else x)
    return df


def split_string(row, column_name, prefix):
    # check if the cell is a string
    if isinstance(row[column_name], str):
        # split the string into a dictionary
        data = dict(item.split(': ') for item in row[column_name].split(','))
        # add the prefix to the keys
        data = {prefix + key: value for key, value in data.items()}
    else:
        # the cell is not a string, return an empty dictionary (empty values are seen as float by pd)
        data = {}
    # convert the dictionary to a Series and add it to the row
    return pd.concat([row, pd.Series(data, dtype = 'object')])


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
parser.add_argument('--mane', metavar='FILE', type=str, help="Path to a file containing the MANE transcripts for filtering.", required=True)
parser.add_argument('--run', metavar='STR', type=str, help="Name of the run analysed", required=True)
parser.add_argument('--pipeline_version', metavar='STR', type=str, help="Version of the reporting pipeline used", required=True)

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
mane:str = pd.read_csv(args.mane)
run_nr:str = args.run
pipeline_version:str = args.pipeline_version

###############################################################
### Loop over .vcf files, extract info and export to report ###
###############################################################

# extract zipped vcf files or read in zipped format
'''
import gzip

with gzip.open('file.gz', 'rb') as f:
    file_content = f.read()

OR if you want to extract all files:

import os
import gzip
import shutil

for filename in os.listdir(input_path):
    if filename.endswith(".gz"):
        with gzip.open(filename, 'rb') as f_in:
            with open(filename[:-3], 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

'''

# Loop over all files in the folder
for filename in os.listdir(input_path):
    # Check if the file is a VCF file
    if filename.endswith(".vcf"):
        # Full path to the VCF file
        vcf_path = os.path.join(input_path, filename)
        vcf = cyvcf2.VCF(vcf_path)

        # Convert VCF to DataFrame
        df = vcf_to_df(vcf)

        # Get sample names
        samples: list[str] = vcf.samples
        if len(samples) != 1:
            raise ValueError(f"Expected exactly one sample in VCF file {filename}, but found {len(samples)} samples.")
        sample_name: str = samples[0]

        # clean up data extracted from info
        # Specify the columns to apply the function to
        # Apply the function to your DataFrame
        df_clean = remove_brackets(df, VCF_INFO_COLUMNS)

        # change annotation of reference sequence to one column
        # output is seen as float, not string (is sometimes empty so pandas consideres it as float)
        df_clean['TRANSCRIPT_A'] = df_clean['TRANSCRIPT_ID_A'].astype(str)
        df_clean = df_clean.drop('TRANSCRIPT_VERSION_A', axis = 1)
        df_clean = df_clean.drop('TRANSCRIPT_ID_A', axis = 1)

        df_clean['TRANSCRIPT_B'] = df_clean['TRANSCRIPT_ID_B'].astype(str)
        df_clean = df_clean.drop('TRANSCRIPT_VERSION_B', axis = 1)
        df_clean = df_clean.drop('TRANSCRIPT_ID_B', axis = 1)

        # Extract the basename and use it as the output Excel file name and worksheet name
        basename = filename.split('_')[0]

        # read info about reads into a new df
        reads_file = os.path.join(reads_path + basename + ".fusions.csv")
        df_reads = pd.read_csv(reads_file)

        # split into different df, depending on algorithm
        df_reads_fc = df_reads.apply(split_string, args=('fusioncatcher', 'fc_'), axis=1)
        df_reads_ar = df_reads.apply(split_string, args=('arriba', 'ar_'), axis=1)
        df_reads_sf = df_reads.apply(split_string, args=('starfusion', 'sf_'), axis=1)

        # split into distinct columns so the positions can be reads
        df_reads_fc['fc_position'] = df_reads_fc.get('fc_position', default=pd.Series(['nan'] * len(df_reads_fc))).astype(str)
        df_reads_fc[['CHRA', 'POSA', 'junk', 'POSB', 'junk2']] = df_reads_fc['fc_position'].apply(lambda x: pd.Series(x.split(':')))
        df_reads_fc = df_reads_fc.drop(['CHRA', 'junk', 'junk2'], axis = 1)

        df_reads_ar['ar_position'] = df_reads_ar.get('ar_position', default=pd.Series(['nan'] * len(df_reads_ar))).astype(str)
        df_reads_ar[['CHRA', 'POSA', 'CHRB', 'POSB']] = df_reads_ar['ar_position'].apply(lambda x: pd.Series(re.split('[:#]', x)))
        df_reads_ar = df_reads_ar.drop(['CHRA', 'CHRB'], axis = 1)

        df_reads_sf['sf_position'] = df_reads_sf.get('sf_position', default=pd.Series(['nan'] * len(df_reads_sf))).astype(str)
        df_reads_sf[['CHRA', 'POSA', 'junk', 'POSB', 'junk2']] = df_reads_sf['sf_position'].apply(lambda x: pd.Series(x.split(':')))
        df_reads_sf = df_reads_sf.drop(['CHRA', 'junk', 'junk2'], axis = 1)

        # combine fusion df with reads df on 'Fusion' column
        df_clean['Fusion'] = df_clean['GENEA'] + "--" + df_clean['GENEB']

        # merge read info one by one into the overview
        merged_df1 = pd.merge(df_clean, df_reads_fc, on = ['Fusion', 'POSA', 'POSB'], how = 'left')
        merged_df2 = pd.merge(df_clean, df_reads_ar, on = ['Fusion', 'POSA', 'POSB'], how = 'left')
        merged_df3 = pd.merge(df_clean, df_reads_sf, on = ['Fusion', 'POSA', 'POSB'], how = 'left')

        # concatenate the merged df
        merged_df = pd.concat([merged_df1, merged_df2, merged_df3], axis = 1)

        # drop duplicate columns
        merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]

        # update chromsome names to have chr in front
        merged_df["CHRA"] = merged_df["CHRA"].apply(update_chroms)
        merged_df["CHRB"] = merged_df["CHRB"].apply(update_chroms)

        # drop original fusioncatcher, arribe and starfusion columns
        merged_df = merged_df.drop(['fusioncatcher', 'arriba', 'starfusion'], axis = 1)

        # set certain columns to numerical
        list_numerical = ["POSA", "POSB", 'fc_longest_anchor', 'fc_common_mapping_reads', 'fc_spanning_pairs',
                          'fc_spanning_unique_reads', 'ar_coverage1', 'ar_coverage2', 'ar_discordant_mates', 'ar_split_reads1',
                          'ar_split_reads2', 'sf_ffmp', 'sf_junction_reads', 'sf_spanning_reads']

        for column in list_numerical:
            merged_df[column] = pd.to_numeric(merged_df[column], errors = 'coerce')

        # limit sf_ffmp to one decimal
        merged_df['sf_ffmp'] = round(merged_df['sf_ffmp'], 1)

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
            "TRANSCRIPT_B",
            "fc_common_mapping_reads",
            "fc_fusion_type",
            "fc_longest_anchor",
            "fc_position",
            "fc_spanning_pairs",
            "fc_spanning_unique_reads",
            "ar_confidence",
            "ar_coverage1",
            "ar_coverage2",
            "ar_discordant_mates",
            "ar_position",
            "ar_reading-frame",
            "ar_split_reads1",
            "ar_split_reads2",
            "ar_type",
            "sf_ffmp",
            "sf_junction_reads",
            "sf_position",
            "sf_spanning_reads"
        ]

        tool_specific_columns: list[str] = [
            "fc_common_mapping_reads",
            "fc_fusion_type",
            "fc_longest_anchor",
            "fc_position",
            "fc_spanning_pairs",
            "fc_spanning_unique_reads",
            "ar_confidence",
            "ar_coverage1",
            "ar_coverage2",
            "ar_discordant_mates",
            "ar_position",
            "ar_reading-frame",
            "ar_split_reads1",
            "ar_split_reads2",
            "ar_type",
            "sf_ffmp",
            "sf_junction_reads",
            "sf_position",
            "sf_spanning_reads"
        ]

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
        df_final = df_filt[(df_filt['SCORE'] > 0.18) & (df_filt['TOOL_HITS'] > 1)]

        # ALTERNATIVE: either one of the two conditions is true: gives too many hits
        # df_final = df_filt[(df_filt['SCORE'] > 0.2) | (df_filt['TOOL_HITS'] > 2)]

        # limit filtered to versioned MANE transcript only
        df_MANE = mane['MANE'].astype(str) + "." + mane['Version'].astype(str)
        df_final_MANE = df_final[df_final['TRANSCRIPT_A'].apply(lambda x: any(i in x for i in df_MANE)) & df_final['TRANSCRIPT_B'].apply(lambda x: any(i in x for i in df_MANE))]

        #df_final_MANE = df_final[df_final['TRANSCRIPT_A'].isin(df_MANE) | df_final['TRANSCRIPT_B'].isin(df_MANE)]

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
        # finaal path nog eens controleren van bashCommand!
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

        with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
            df_final.to_excel(writer, index=False, sheet_name="fusions_filt")
            df_final_MANE.to_excel(writer, index = False, sheet_name = "fusions_filt_MANE")
            df_fusions.to_excel(writer, index= False, sheet_name="fusions_specific")
            merged_df.to_excel(writer, index=False, sheet_name="fusions_all")
            splicing.to_excel(writer, index=False, sheet_name = "CTAT splicing")
            ref_genes_expr.to_excel(writer, index= False, sheet_name="coverage_ref")
            genes_expr.to_excel(writer, index=False, sheet_name="coverage_all")
            qc_filt.to_excel(writer, index=False, sheet_name="QC")


        ###############################
        ### Change layout of report ###
        ###############################

        # open report file
        workbook = openpyxl.load_workbook(excel_path)

        # modify fusion worksheets
        for worksheet in ["fusions_all", "fusions_filt", "fusions_filt_MANE", "fusions_specific"]:
            ws = workbook[worksheet]

            # insert sample name in new empty row and add hyperlink to visualisation files
            ws.insert_rows(1)
            ws['A1'] = "fusions " + basename

            # TODO: Nextflow -> Find another way to make this work
            # # format links
            # link_arriba = "https://login.hpc.ugent.be/pun/sys/dashboard/files/fs//data/gent/vo/000/gvo00082/research/LabMDG/RNASeq/" + \
            #     run_nr + "/output/arriba_visualisation/" + basename + "_combined_fusions_arriba_visualisation.pdf"

            # link_fusioninspector = "https://login.hpc.ugent.be/pun/sys/dashboard/files/fs//data/gent/vo/000/gvo00082/research/LabMDG/RNASeq/" + \
            #     run_nr + "/output/fusioninspector/" + basename + ".fusion_inspector_web.html"

            # # check full path to .bam once on server!
            # link_igv = "http://localhost:60151/load?file=" + bam_path + "/" + basename + ".bam&genome=hg38"

            # # change layout so it is recognized as hyperlink
            # # links not functional so far (Excel removes a / after fs, to change once on final server)
            # ws['D1'].hyperlink = link_arriba
            # ws['D1'].value = "Arriba visualisation"
            # ws['D1'].font = Font(color="0000EE", underline="single")

            # ws['F1'].hyperlink = link_fusioninspector
            # ws['F1'].value = "FusionInspector"
            # ws['F1'].font = Font(color="0000EE", underline="single")

            # ws['H1'].hyperlink = link_igv
            # ws['H1'].value = "IGV visualisation"
            # ws['H1'].font = Font(color="0000EE", underline="single")

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

        # TODO: Nextflow -> Find another way to make this work
        # # add extra hyperlinks to the worksheets that are useful (if needed for fusions_all, columns have to be changed)
        # for worksheet in ["fusions_filt", "fusions_filt_MANE", "fusions_specific"]:
        #     ws = workbook[worksheet]

        #     # individual links to positions of fusions (added as links in columns H and I)
        #     for row in ws.iter_rows(min_row=3):
        #         cell_B = row[1]
        #         cell_H = row[7]

        #         bp_link = f"http://localhost:60151/goto?locus=chr{cell_B.value}:{cell_H.value}"
        #         cell_H.hyperlink = bp_link

        #         cell_C = row[2]
        #         cell_I = row[8]
        #         bp_link = f"http://localhost:60151/goto?locus=chr{cell_C.value}:{cell_I.value}"
        #         cell_I.hyperlink = bp_link


        #loop over coverage worksheets to change layout
        for worksheet in ["coverage_ref", "coverage_all"]:
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
        ws = workbook["coverage_all"]

        ws['J26'] = DUX4_reads + " filtered reads"

        # change layout of splicing worksheet
        for worksheet in ["CTAT splicing"]:
            ws = workbook[worksheet]

            # add new line and add sample name
            ws.insert_rows(1)
            ws['A1']= "CTAT splicing " + basename

            # TODO: Nextflow -> Find another way to make this work
            # link_splicing = "https://login.hpc.ugent.be/pun/sys/dashboard/files/fs//data/gent/vo/000/gvo00082/research/LabMDG/RNASeq/" + \
            #     run_nr + "/output/ctat_splicing_arriba/" + basename + ".ctat-splicing.igv.html"

            # # links not functional so far (Excel removes a / after fs, to change once on final server)
            # ws['D1'].hyperlink = link_splicing
            # ws['D1'].value = "CTAT IGV"
            # ws['D1'].font = Font(color="0000EE", underline="single")

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
        for worksheet in ["QC"]:
            ws = workbook[worksheet]

            # add new line and add sample name
            ws.insert_rows(1)
            ws['A1']= "QC " + basename

            # add script version for logging, to be replaced in the future
            ws['D1'] = "nf-cmgg/report version: " + pipeline_version

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

        # save file
        workbook.save(excel_path)

# remove unzipped .vcf files again / not necessary if zipped files have been read
