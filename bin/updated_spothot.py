#!/usr/bin/python3

import sys
import argparse
import configparser
import os.path
import grp
import shutil
import subprocess
import xlsxwriter
import re
import glob

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
    reads a section of the config file and returns it as a dictionary.
    raise keyerror if the section does not exist
    """
    if section not in config_obj:
        raise KeyError(f"Section '{section}' not found in config file.")
    return dict(config_obj[section])

def read_panel_file(file, paneldictionary):
    """
    extract panel and design info from the panellist
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
    Extract all patients from a specific design from HyperCap runinfo file
    """
    with open(file, mode = "r") as f2:
        for line in f2:
            items = line.rstrip("\n").split("\t")
            patient = items[0]
            first_panel = items[1].split(',')[0]
            if paneldictionary[first_panel] == design_query:
                patients.append(patient)
    return patients

### user and group asked but nowhere used, can be removed? ###
def count_reads_at_target(crams, outputdir, pos, user, group, samtools, hotcount, pear, log, output, query):
    """
    count reads with certain variant/wildtype sequence overlapping a target position
    """
    sample_coverage_to_low = []

    for sample, cram in crams.items():
        ### Log sample that is being processed ###
        print(f"\t{sample}")
        subprocess.check_output(f"echo\t{sample} >> {log}", stderr=subprocess.STDOUT, shell=True)
        
        ### Defining file paths for intermediate files ###
        bam = f"{outputdir}/{sample}_var.bam"
        bam_sorted = f"{outputdir}/{sample}_var_sorted.bam"
        fq_r1 = f"{outputdir}/{sample}_R1.fastq"
        fq_r2 = f"{outputdir}/{sample}_R2.fastq"
        fq_p = f"{outputdir}/{sample}_paired"
        fq_s = f"{outputdir}/{sample}_singletons.fastq"
        fq_all = f"{outputdir}/{sample}_all.fastq.gz"

        ### extract reads from BAM/CRAM ###
        subprocess.check_output(f"{samtools} view -T {reference_genome} {cram} {pos} -F 12 -b -o {bam}", stderr=subprocess.STDOUT, shell=True)
        subprocess.check_output(f"{samtools} sort -n -o {bam_sorted} {bam}", stderr=subprocess.STDOUT, shell=True)
        subprocess.check_output(f"{samtools} fastq -n -1 {fq_r1} -2 {fq_r2} -s {fq_s} {bam_sorted} 2>> {log}", shell=True)

        ### check if file is not empty (otherwise pear goes in ERROR) ###
        if os.path.getsize(fq_r1) != 0:
            ### merge overlapping read pairs using pear ###
            subprocess.check_output(f"{pear} --forward-fastq {fq_r1} --reverse-fastq {fq_r2} --output {fq_p} >> {log}", stderr=subprocess.STDOUT, shell=True)
            ### count target sequences in merged/singleton reads ###
            subprocess.check_output(f"cat {fq_p}.* {fq_s} | gzip -c > {fq_all}", stderr=subprocess.STDOUT, shell=True)
            subprocess.check_output(f"{hotcount} {query} {fq_all} >> {output}", shell=True)
            subprocess.check_output(f"rm {bam} {fq_r1} {fq_r2} {fq_p}* {fq_s}", shell=True)
        else:
            sample_coverage_to_low.append(sample)
            subprocess.check_output(f"rm {bam} {fq_r1} {fq_r2} {fq_s}", shell=True)
    return sample_coverage_to_low

def parse_count_output(outp, var_dict):
    """
    Extract counts from hotcount output
    """
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
results_dir = f"{directories["documents"]}{directories["results"]}_{analysis["variant_caller"]}" 
runinfofile = files["runinfofile"]
bin = directories["bin"]
samtools_path = os.path.expanduser(targeted_variants["samtools"])
hotcount_path = targeted_variants["hotcount"]
pear_path = os.path.expanduser(targeted_variants["pear"])
msh2_hotspot_var = targeted_variants["msh2"].split(',')
vaf_detection_treshold = float(targeted_variants["vaf_detection_threshold"])
reference_genome = files["reference_genome"]

### error messages ###
if not os.path.exists(results_dir):
    print(f'Map "{results_dir}" bestaat niet')
    os.makedirs(results_dir)
    
if not os.path.exists(runinfofile):
    print(f'Bestand "{runinfofile}" bestaat niet')

### get user and group id ###
user_id = os.getuid()
group_id = grp.getgrnam("seqcap").gr_gid

### read panel file ###
panel_vs_design = read_panel_file(panelfile, {})
samples_coverage_to_low = []

for designs_plus_gene in designs_list:
    ### variables ###
    design, rest = designs_plus_gene.split('~')
    gene, position = rest.split('/')  
    hotspot_count = {}
    query_list = f"{design}_{gene}.txt"

    ### read runinfo ###
    patients_design = read_runinfo_file(runinfofile, panel_vs_design, [], design)
    print(f"\n*** Start analyse voor {len(patients_design)} stalen ***")
    if patients_design:
        ### check output ###
        excel_output_sub_dir = f"{results_dir}/{design}_targeted"
        if os.path.exists(excel_output_sub_dir):
            print(f"-> !!! Opgelet: de map {excel_output_sub_dir} bestond al en wordt overschreven !!!")
            shutil.rmtree(excel_output_sub_dir)
        if not os.path.exists(excel_output_sub_dir):
            os.makedirs(excel_output_sub_dir)
            print(f"-> Map \"{excel_output_sub_dir}\" gemaakt")
            os.chmod(excel_output_sub_dir, 0o775)
            os.chown(excel_output_sub_dir, user_id, group_id)
        print(f"-> {len(patients_design)} {design} patienten te analyseren")
        
        ### read variant query file ###
        if not os.path.isfile(query_list):
            print(f"Error: {gene} variant bestand \"{query_list}\" kon niet worden geopend\n")
            sys.exit()

        ### retrieving CRAM files ###
        print("-> Ophalen CRAM files")
        MAPPING_files = {
            patient: cram
            for patient in patients_design
            for cram in [os.path.abspath(f"{patient}-ready.cram"),
                         os.path.abspath(f"{patient}-ready.bam")]
            if (print(f"Checking {cram}... Exists: {os.path.exists(cram)}") or os.path.exists(cram))
        }
        print(MAPPING_files)
        if len(MAPPING_files) < len(patients_design):
            missing = set(patients_design) - set(MAPPING_files)
            print(f"Error: CRAM/BAM files missing for: {', '.join(missing)}")
            sys.exit()
        
        ### counting reads with specific target sequences ###
        print(f"-> Tellen van reads in positie {gene} {position}:\n")
        logfile = f"{excel_output_sub_dir}/{run}_{design}_{gene}_varcount.log"
        outputfile = f"{excel_output_sub_dir}/{run}_{design}_{gene}_varcount.txt"
        output_excel = f"{results_dir}/{run}_{design}_targeted_{gene}.xlsx"
        output_txt = f"{excel_output_sub_dir}/{run}_{design}_targeted_{gene}.txt"
        samples_coverage_to_low = count_reads_at_target(
            MAPPING_files, excel_output_sub_dir, position, user_id, group_id, 
            samtools_path, hotcount_path, pear_path, logfile, outputfile, query_list
        )
        if samples_coverage_to_low:
            print(f"\tWARNING: volgende stalen hadden geen coverage in de ROI van de variant: {samples_coverage_to_low}")
        hotspot_count = parse_count_output(outputfile, hotspot_count)

        ### Write to excel and txt output file ###
        with open(output_txt, "w+") as txt:
            workbook = xlsxwriter.Workbook(output_excel)
            worksheet = workbook.add_worksheet(f"{design} patienten")
            worksheet.set_column('A:B', 16)
            worksheet.set_column('C:F', 25)

            format_overview_title = workbook.add_format({'bold': 1, 'font_size': 18, 'font_color': 'black', 'align': 'left'})
            run_overview_format = workbook.add_format({'bold': 1, 'font_size': 14, 'font_color': 'black','bg_color': '#D8D8D8', 'align': 'centre'})
            header_overview_format = workbook.add_format({'bold': 1, 'font_size': 11, 'font_color': 'black', 'bg_color': '#D8D8D8', 'align': 'left'})
            body_overview_format = workbook.add_format({'bold': 0, 'font_size': 11, 'font_color': 'black', 'align': 'left'})
            body_overview_format_bg = workbook.add_format({'bold': 1, 'font_size': 11, 'font_color': 'black', 'bg_color': '#4db8ff', 'align': 'left'})

            if gene == 'MSH2':
                worksheet.write(0, 0, f"Controle van {gene} hotspot varianten {','.join(msh2_hotspot_var)}", format_overview_title)
            worksheet.write(1, 0, f"Run {run}", run_overview_format)
            worksheet.write(4, 0, "DNA-nummer", header_overview_format)
            worksheet.write(4, 1, "WT count", header_overview_format)
            my_header = "DNA-nummer\tWT count"
            for column, hotspot_var in enumerate(msh2_hotspot_var, start=2):
                label = hotspot_var.split(':')[1] if ':' in hotspot_var else hotspot_var
                worksheet.write(4, column, f"{label} frequentie (count)", header_overview_format)
                my_header += f"\t{label} frequentie (count)"
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
                        frequency = round(float(alt) / (float(alt) + float(wt)) * 100, 2)
                        if frequency >= vaf_detection_treshold:
                            freq_str = f"{frequency}({alt}) - Sangeren!"
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

            for patient in samples_coverage_to_low:
                worksheet.write(rownum, 0, patient, body_overview_format)
                worksheet.write(rownum, 1, "0", body_overview_format)
                my_row = f"{patient}\t0"
                for colnum in range(2, 2 + len(msh2_hotspot_var)):
                    worksheet.write(rownum, colnum, "0", body_overview_format)
                    my_row += "\t0(0)"
                txt.write(my_row + "\n")
                rownum += 1
            workbook.close()
            txt.close()
        os.chmod(output_excel, 0o664)
        os.chown(output_excel, user_id, group_id)
        os.chmod(output_txt, 0o664)
        os.chown(output_txt, user_id, group_id)
        print(f"\nDe analyse is afgelopen ({output_excel})\n")
    else:
        print(f"->Geen patienten aanwezig voor design {design}")
