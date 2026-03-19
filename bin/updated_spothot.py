#!/usr/bin/env python3

import sys
import argparse
import configparser
import os
import grp
import shutil
import subprocess
import xlsxwriter
import re
import glob
import gzip

#########################
### read config file  ###
#########################
config_file = "./config_seqcap_v2.ini"
config = configparser.ConfigParser()
config.read(config_file)

#########################
###     functions     ###
#########################

def config_section_map(section, config_obj):
    """
    Reads a section of the config file and returns it as a dictionary.
    Raise KeyError if the section does not exist.
    """
    if section not in config_obj:
        raise KeyError(f"Section '{section}' not found in config file.")
    return dict(config_obj[section])

def read_panel_file(file, paneldictionary):
    """
    Extract panel and design info from the panellist.
    """
    with open(file, "r") as f1:
        for line in f1:
            items = line.rstrip("\n").split("\t")
            panel = items[0]
            design = items[1]
            paneldictionary[panel] = design
    return paneldictionary

def read_runinfo_file(file, paneldictionary, patients, design_query):
    """
    Extract all patients from a specific design from HyperCap runinfo file.
    """
    with open(file, mode="r") as f2:
        for line in f2:
            items = line.rstrip("\n").split("\t")
            patient = items[0]
            first_panel = items[1].split(',')[0]
            if paneldictionary.get(first_panel) == design_query:
                patients.append(patient)
    return patients

### user and group asked but nowhere used, can be removed? ###
def count_reads_at_target(crams, outputdir, pos, user, group, samtools, hotcount, pear, log, output, query):
    """
    Count reads with certain variant/wildtype sequence overlapping a target position.
    """
    sample_coverage_too_low = []

    for sample, cram in crams.items():
        ### Log sample that is being processed ###
        print(f"\t{sample}")
        with open(log, "a") as log_file:
            log_file.write(f"\t{sample}\n")
        
        ### Defining file paths for intermediate files ###
        bam = os.path.join(outputdir, f"{sample}_var.bam")
        bam_sorted = os.path.join(outputdir, f"{sample}_var_sorted.bam")
        fq_r1 = os.path.join(outputdir, f"{sample}_R1.fastq")
        fq_r2 = os.path.join(outputdir, f"{sample}_R2.fastq")
        fq_p = os.path.join(outputdir, f"{sample}_paired")
        fq_s = os.path.join(outputdir, f"{sample}_singletons.fastq")
        fq_all = os.path.join(outputdir, f"{sample}_all.fastq.gz")

        ### extract reads from BAM/CRAM ###
        subprocess.run([samtools, "view", "-T", reference_genome, cram, pos, "-F", "12", "-b", "-o", bam], check=True, stderr=subprocess.STDOUT)
        subprocess.run([samtools, "sort", "-n", "-o", bam_sorted, bam], check=True, stderr=subprocess.STDOUT)
        
        with open(log, "a") as log_file:
            subprocess.run([samtools, "fastq", "-n", "-1", fq_r1, "-2", fq_r2, "-s", fq_s, bam_sorted], check=True, stderr=log_file)

        ### check if file is not empty (otherwise pear goes in ERROR) ###
        if os.path.exists(fq_r1) and os.path.getsize(fq_r1) != 0:
            ### merge overlapping read pairs using pear ###
            with open(log, "a") as log_file_out:
                subprocess.run(
                    [pear, "--forward-fastq", fq_r1, "--reverse-fastq", fq_r2, "--output", fq_p], 
                    check=True, 
                    stdout=log_file_out, 
                    stderr=subprocess.STDOUT
                )
            
            ### count target sequences in merged/singleton reads ###
            peared_files = glob.glob(f"{fq_p}.*") + [fq_s]
            with gzip.open(fq_all, 'wb') as f_out:
                for file_to_cat in peared_files:
                    if os.path.exists(file_to_cat):
                        with open(file_to_cat, 'rb') as f_in:
                            shutil.copyfileobj(f_in, f_out)
            
            with open(output, "a") as out_file:
                subprocess.run([hotcount, query, fq_all], check=True, stdout=out_file)
            
            files_to_remove = [bam, bam_sorted, fq_r1, fq_r2, fq_s] + glob.glob(f"{fq_p}*")
            for file_to_remove in files_to_remove:
                if os.path.exists(file_to_remove) and file_to_remove != fq_all:
                    os.remove(file_to_remove)
        else:
            sample_coverage_too_low.append(sample)
            for file_to_remove in [bam, bam_sorted, fq_r1, fq_r2, fq_s]:
                if os.path.exists(file_to_remove):
                    os.remove(file_to_remove)

    return sample_coverage_too_low

def parse_count_output(outp, var_dict):
    """
    Extract counts from hotcount output.
    """
    if not os.path.exists(outp):
        return var_dict
        
    with open(outp, "r") as f3:
        for line in f3:
            if line.startswith('/'):
                dnanumber = line.split(' ')[0].split('/')[-1].split('_all.fastq')[0]
                var_dict[dnanumber] = tuple(line.split(' ')[1:])
    return var_dict

################################
###     input parameters     ###
################################

parser = argparse.ArgumentParser(usage='Usage: python SpotHot.py --run NVQ_828')
parser.add_argument('--run', type=str, required=True, help='run name')
args = parser.parse_args()
run = args.run

### config sections ###
targeted_variants = config_section_map("targeted_variants", config)
directories = config_section_map("dir", config)
files = config_section_map("files", config)
analysis = config_section_map("analysis", config)

### config keys ###
designs = targeted_variants["designs"]
designs_list = designs.split(',')
panelfile = files["seqcap_panels_file"]
results_dir = f"{directories['documents']}{directories['results']}_{analysis['variant_caller']}"
runinfofile = files["runinfofile"]
bin_dir = directories["bin"]
samtools_path = os.path.expanduser(targeted_variants["samtools"])
hotcount_path = targeted_variants["hotcount"]
pear_path = os.path.expanduser(targeted_variants["pear"])
msh2_hotspot_var = targeted_variants["msh2"].split(',')
vaf_detection_threshold = float(targeted_variants["vaf_detection_threshold"])
reference_genome = files["reference_genome"]

### error messages ###
if not os.path.exists(results_dir):
    print(f'Directory "{results_dir}" does not exist, creating it.')
    os.makedirs(results_dir)

if not os.path.exists(runinfofile):
    print(f'File "{runinfofile}" does not exist.')

### get user and group id ###
user_id = os.getuid()
try:
    group_id = grp.getgrnam("seqcap").gr_gid
except KeyError:
    group_id = os.getgid()

### read panel file ###
panel_vs_design = read_panel_file(panelfile, {})
samples_coverage_too_low = []

for designs_plus_gene in designs_list:
    ### variables ###
    design, rest = designs_plus_gene.split('~')
    gene, position = rest.split('/')
    hotspot_count = {}
    query_list = f"{design}_{gene}.txt"

    ### read runinfo ###
    patients_design = read_runinfo_file(runinfofile, panel_vs_design, [], design)
    print(f"\n*** Start analysis for {len(patients_design)} samples ***")
    
    if patients_design:
        ### check output ###
        excel_output_sub_dir = f"{results_dir}/{design}_targeted"
        if os.path.exists(excel_output_sub_dir):
            print(f"-> !!! Warning: directory '{excel_output_sub_dir}' already existed and will be overwritten !!!")
            shutil.rmtree(excel_output_sub_dir)
            
        os.makedirs(excel_output_sub_dir)
        print(f"-> Directory \"{excel_output_sub_dir}\" created")
        os.chmod(excel_output_sub_dir, 0o775)
        os.chown(excel_output_sub_dir, user_id, group_id)
        
        print(f"-> {len(patients_design)} {design} patients to analyze")
        
        ### read variant query file ###
        if not os.path.isfile(query_list):
            print(f"Error: {gene} variant file \"{query_list}\" could not be opened\n")
            sys.exit(1)

        ### retrieving CRAM files ###
        print("-> Retrieving CRAM files")
        MAPPING_files = {}
        for patient in patients_design:
            cram_path = os.path.abspath(f"{patient}-ready.cram")
            bam_path = os.path.abspath(f"{patient}-ready.bam")
            if os.path.exists(cram_path):
                print(f"Found {cram_path}")
                MAPPING_files[patient] = cram_path
            elif os.path.exists(bam_path):
                print(f"Found {bam_path}")
                MAPPING_files[patient] = bam_path
            else:
                print(f"Missing {cram_path} or {bam_path}")
                
        if len(MAPPING_files) < len(patients_design):
            missing = set(patients_design) - set(MAPPING_files)
            print(f"Error: CRAM/BAM files missing for: {', '.join(missing)}")
            sys.exit(1)
            
        ### counting reads with specific target sequences ###
        print(f"-> Counting reads at position {gene} {position}:\n")
        logfile = f"{excel_output_sub_dir}/{run}_{design}_{gene}_varcount.log"
        outputfile = f"{excel_output_sub_dir}/{run}_{design}_{gene}_varcount.txt"
        output_excel = f"{results_dir}/{run}_{design}_targeted_{gene}.xlsx"
        output_txt = f"{excel_output_sub_dir}/{run}_{design}_targeted_{gene}.txt"
        
        samples_coverage_too_low = count_reads_at_target(
            crams=MAPPING_files,
            outputdir=excel_output_sub_dir,
            pos=position,
            user=user_id,
            group=group_id,
            samtools=samtools_path,
            hotcount=hotcount_path,
            pear=pear_path,
            log=logfile,
            output=outputfile,
            query=query_list
        )
        
        if samples_coverage_too_low:
            print(f"\tWARNING: the following samples had no coverage in the ROI of the variant: {samples_coverage_too_low}")
            
        hotspot_count = parse_count_output(outputfile, hotspot_count)

        ### Write to excel and txt output file ###
        with open(output_txt, "w+") as txt:
            workbook = xlsxwriter.Workbook(output_excel)
            worksheet = workbook.add_worksheet(f"{design} patients")
            worksheet.set_column('A:B', 16)
            worksheet.set_column('C:F', 25)

            format_overview_title = workbook.add_format({'bold': 1, 'font_size': 18, 'font_color': 'black', 'align': 'left'})
            run_overview_format = workbook.add_format({'bold': 1, 'font_size': 14, 'font_color': 'black', 'bg_color': '#D8D8D8', 'align': 'center'})
            header_overview_format = workbook.add_format({'bold': 1, 'font_size': 11, 'font_color': 'black', 'bg_color': '#D8D8D8', 'align': 'left'})
            body_overview_format = workbook.add_format({'bold': 0, 'font_size': 11, 'font_color': 'black', 'align': 'left'})
            body_overview_format_bg = workbook.add_format({'bold': 1, 'font_size': 11, 'font_color': 'black', 'bg_color': '#4db8ff', 'align': 'left'})

            if gene == 'MSH2':
                worksheet.write(0, 0, f"Check of {gene} hotspot variants {','.join(msh2_hotspot_var)}", format_overview_title)
            worksheet.write(1, 0, f"Run {run}", run_overview_format)
            worksheet.write(4, 0, "DNA-number", header_overview_format)
            worksheet.write(4, 1, "WT count", header_overview_format)
            my_header = "DNA-number\tWT count"
            
            for column, hotspot_var in enumerate(msh2_hotspot_var, start=2):
                label = hotspot_var.split(':')[1] if ':' in hotspot_var else hotspot_var
                worksheet.write(4, column, f"{label} frequency (count)", header_overview_format)
                my_header += f"\t{label} frequency (count)"
            txt.write(my_header + "\n")
            rownum = 5

            ### sort on 2nd key of tuple, reverse order ###
            for patient, my_tuple in sorted(hotspot_count.items(), key=lambda item: int(item[1][1]), reverse=True):
                wt, *alt_list = my_tuple
                worksheet.write(rownum, 0, patient, body_overview_format)
                worksheet.write(rownum, 1, wt, body_overview_format)
                my_row = f"{patient}\t{wt}"
                
                for i, hotspot_var in enumerate(msh2_hotspot_var):
                    colnum = i + 2
                    try:
                        alt = alt_list[i]
                        # Safe division against ZeroDivisionError
                        total_reads = float(alt) + float(wt)
                        frequency = round(float(alt) / total_reads * 100, 2) if total_reads > 0 else 0.0
                        
                        if frequency >= vaf_detection_threshold:
                            freq_str = f"{frequency}({alt}) - Perform Sanger!"
                            worksheet.write(rownum, colnum, freq_str, body_overview_format_bg)
                        else:
                            worksheet.write(rownum, colnum, f"{frequency}")
                        my_row += f"\t{frequency}({alt})"
                    except IndexError:
                        # If alt_list is too short, write zero
                        worksheet.write(rownum, colnum, "0", body_overview_format)
                        my_row += "\t0(0)"
                txt.write(my_row + "\n")
                rownum += 1

            for patient in samples_coverage_too_low:
                worksheet.write(rownum, 0, patient, body_overview_format)
                worksheet.write(rownum, 1, "0", body_overview_format)
                my_row = f"{patient}\t0"
                for colnum in range(2, 2 + len(msh2_hotspot_var)):
                    worksheet.write(rownum, colnum, "0", body_overview_format)
                    my_row += "\t0(0)"
                txt.write(my_row + "\n")
                rownum += 1
            
            workbook.close()
            
        os.chmod(output_excel, 0o664)
        os.chown(output_excel, user_id, group_id)
        os.chmod(output_txt, 0o664)
        os.chown(output_txt, user_id, group_id)
        print(f"\nThe analysis is finished ({output_excel})\n")
    else:
        print(f"-> No patients found for design {design}")
