#!/usr/bin/env python3
"""
exomecnv_to_csv.py

Summarizes ExomeDepth CNV calls (VEP-annotated VCFs) for a SeqCap run into
per-design and per-panel TSV reports, annotated with recurrency information
from a run-level sample set and a historical CMGG high-confidence CNV
database.

Designed for standalone execution and integration into Nextflow pipelines.

Inputs:
  1. samplesheet.csv: Comma-separated file with columns sample, design, panel,
     design_genelist, panel_genelist, vcf, cnv_database.

Outputs:
  - cnvs_exomedepth_designs_summary.tsv: One row per patient/design/CNV.
  - cnvs_exomedepth_panels_summary.tsv: One row per patient/panel/design/CNV.
"""

import argparse
import csv
import gzip
import os
import re
import shutil

from collections import defaultdict


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Summarize ExomeDepth CNVs from a samplesheet."
    )
    parser.add_argument(
        "--samplesheet",
        required=True,
        help="CSV samplesheet with columns: sample, design, panel, design_genelist, panel_genelist, vcf, cnv_database.",
    )
    parser.add_argument(
        "--output-dir",
        default="results_exomedepth",
        help="Directory for TSV output files (default: ./results_exomedepth).",
    )
    return parser.parse_args()


def build_output_header(include_panel=False):
    header = ["Patient"]
    if include_panel:
        header.append("Panel")
    header.extend([
        "Design",
        "Type",
        "Region",
        "Recurrency_run",
        "Recurrency_HC-CNV",
        "Exon count",
        "BF",
        "Reads.expected",
        "Reads.observed",
        "Reads.ratio",
        "Exons",
        "HC-CNV.reads_ratio.[0-0.2[",
        "HC-CNV.reads_ratio.[0.2-1[",
        "HC-CNV.reads_ratio.[1-1.83[",
        "HC-CNV.reads_ratio.>=1.83",
        "HC-CNV.patient_included",
        "HC-CNV.patients_subsample",
    ])
    return header

def _resolve_samplesheet_path(path, samplesheet_dir):
    path = os.path.expanduser(path.strip())
    if not os.path.isabs(path):
        path = os.path.join(samplesheet_dir, path)
    return os.path.abspath(path)


def _read_gene_list(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Gene list file not found: {path}")
    with open(path, encoding="utf-8-sig") as handle:
        genes = {line.strip() for line in handle if line.strip()}
    return {gene for gene in genes if not gene.startswith("#")}


def read_samplesheet(samplesheet_path):
    """Load samples and per-design/panel resources from a Nextflow samplesheet."""
    if not os.path.isfile(samplesheet_path):
        raise FileNotFoundError(f"Samplesheet not found: {samplesheet_path}")

    required_columns = {
        "sample", "design", "panel", "design_genelist", "panel_genelist",
        "vcf", "cnv_database",
    }
    samples = {}
    design_gene_paths = {}
    panel_gene_paths = {}
    database_paths = {}
    samplesheet_dir = os.path.dirname(os.path.abspath(samplesheet_path))

    print(f"Reading samplesheet: {samplesheet_path}")
    with open(samplesheet_path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Empty or invalid samplesheet: {samplesheet_path}")
        field_map = {name.strip().lower(): name for name in reader.fieldnames if name}
        missing = required_columns - field_map.keys()
        if missing:
            raise ValueError(
                "Samplesheet is missing required column(s): " + ", ".join(sorted(missing))
            )

        for row_idx, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(f"Samplesheet row {row_idx} has extra fields.")
            values = {
                name: (row.get(field_map[name]) or "").strip()
                for name in required_columns
            }
            if not any(values.values()):
                continue
            for name in required_columns:
                if not values[name]:
                    raise ValueError(
                        f"Samplesheet row {row_idx} is missing a value for '{name}'."
                    )

            sample = values["sample"]
            design = values["design"]
            panel = values["panel"]
            vcf = _resolve_samplesheet_path(values["vcf"], samplesheet_dir)
            design_genelist = _resolve_samplesheet_path(
                values["design_genelist"], samplesheet_dir
            )
            panel_genelist = _resolve_samplesheet_path(
                values["panel_genelist"], samplesheet_dir
            )
            cnv_database = _resolve_samplesheet_path(
                values["cnv_database"], samplesheet_dir
            )
            for path in (vcf, design_genelist, panel_genelist, cnv_database):
                if not os.path.isfile(path):
                    raise FileNotFoundError(f"Samplesheet row {row_idx}: file not found: {path}")

            if design in design_gene_paths and design_gene_paths[design] != design_genelist:
                raise ValueError(f"Design '{design}' uses multiple design_genelist files.")
            if design in database_paths and database_paths[design] != cnv_database:
                raise ValueError(f"Design '{design}' uses multiple cnv_database files.")
            panel_key = (design, panel)
            if panel_key in panel_gene_paths and panel_gene_paths[panel_key] != panel_genelist:
                raise ValueError(f"Panel '{panel}' in design '{design}' uses multiple gene lists.")
            design_gene_paths[design] = design_genelist
            database_paths[design] = cnv_database
            panel_gene_paths[panel_key] = panel_genelist

            if sample not in samples:
                samples[sample] = {
                    "design": design,
                    "vcf": vcf,
                    "panels": set(),
                }
            elif samples[sample]["design"] != design:
                raise ValueError(f"Sample '{sample}' uses more than one design.")
            elif samples[sample]["vcf"] != vcf:
                raise ValueError(f"Sample '{sample}' uses more than one VCF file.")
            samples[sample]["panels"].add(panel)

    design_gene_sets = {}
    for design, path in design_gene_paths.items():
        print(f"Reading gene list for design '{design}': {path}")
        design_gene_sets[design] = _read_gene_list(path)
    panel_gene_sets = {}
    for (design, panel), path in panel_gene_paths.items():
        print(f"Reading gene list for panel '{panel}' in design '{design}': {path}")
        panel_gene_sets[(design, panel)] = _read_gene_list(path)
    cnvs_databases = {}
    cnvs_databases_count_per_design = {}
    for design, path in database_paths.items():
        print(f"Reading CNV database for design '{design}': {path}")
        cnvs_databases[design], cnvs_databases_count_per_design[design] = (
            read_cmgg_hc_cnv_database(path)
        )

    return samples, design_gene_sets, panel_gene_sets, cnvs_databases, cnvs_databases_count_per_design


def _parse_patients_reads_ratio(field_value):
    """
    Parse value like: D2309982:1.45;D2311428:1.42
    into: {"D2309982": "1.45", "D2311428": "1.42"}
    """
    parsed = {}
    if not field_value:
        return parsed

    for token in field_value.split(";"):
        token = token.strip()
        if not token or ":" not in token:
            continue
        patient, ratio = token.split(":", 1)
        patient = patient.strip()
        ratio = ratio.strip()
        if patient:
            parsed[patient] = ratio

    return parsed


def read_cmgg_hc_cnv_database(database_path):
    """
    Equivalent of Perl read_cmgg_hc_cnv_database():
    returns (cnv_map, total_patient_count)

    cnv_map[cnv_id] contains keys used later in reporting:
    - patient_count
    - total_patient_count
    - [0-0.2[
    - [0.2-1[
    - [1-1.83[
    - >=1.83
    - patients_subsample
    - patients_reads_ratio (dict patient -> ratio)
    """
    cnv_map = {}
    total_patient_count = 0

    with open(database_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cnv_id = (row.get("id") or "").strip()
            if not cnv_id:
                continue

            total_str = (row.get("total_patient_count") or "").strip()
            if total_str.isdigit():
                total_patient_count = int(total_str)

            patients_reads_ratio_raw = (row.get("patients_reads_ratio") or "").strip()
            patients_reads_ratio = _parse_patients_reads_ratio(patients_reads_ratio_raw)
            patients_subsample = ";".join(list(patients_reads_ratio.keys())[:15])

            cnv_map[cnv_id] = {
                "patient_count": (row.get("patient_count") or "0").strip(),
                "total_patient_count": total_str if total_str else str(total_patient_count),
                "[0-0.2[": (row.get("[0-0.2[") or "0").strip(),
                "[0.2-1[": (row.get("[0.2-1[") or "0").strip(),
                "[1-1.83[": (row.get("[1-1.83[") or "0").strip(),
                ">=1.83": (row.get(">=1.83") or "0").strip(),
                "patients_subsample": patients_subsample,
                "patients_reads_ratio": patients_reads_ratio,
            }

    return cnv_map, total_patient_count


def _parse_info_field(info_str):
    info = {}
    for token in info_str.split(";"):
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            info[key] = value
        else:
            info[token] = True
    return info


def _extract_genes_from_csq(info, gene_field_index=0):
    genes = set()
    csq = info.get("CSQ", "")
    if not csq:
        return genes

    for csq_item in csq.split(","):
        fields = csq_item.split("|")
        if gene_field_index < len(fields):
            gene = fields[gene_field_index].strip()
            if gene and gene not in {".", "-"}:
                genes.add(gene)
    return genes


def iter_vcf_records(vcf):
    svtype_map = {
        "DUP": "duplication",
        "DEL": "deletion",
    }
    csq_gene_index = 0

    with gzip.open(vcf, "rt") as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ,"):
                match = re.search(r'Format:\s*([^">]+)', line)
                if match:
                    csq_fields = match.group(1).split("|")
                    for field in ("SYMBOL", "Gene"):
                        if field in csq_fields:
                            csq_gene_index = csq_fields.index(field)
                            break
                continue
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _id, _ref, alt, _qual, _flt, info_str = parts[:8]
            info = _parse_info_field(info_str)
            svtype = info.get("SVTYPE", alt.strip("<>"))
            type_label = svtype_map.get(svtype, svtype)
            end = info.get("END", pos)
            region = f"{chrom}:{pos}-{end}"

            yield {
                "type": type_label,
                "region": region,
                "nexons": info.get("NEXONS", "NA"),
                "bf": info.get("BF", "NA"),
                "reads_expected": info.get("READSEXP", "NA"),
                "reads_observed": info.get("READSOBS", "NA"),
                "reads_ratio": info.get("READRATIO", "NA"),
                "exons": info.get("EXONSENS", "NA"),
                "genes": _extract_genes_from_csq(info, csq_gene_index),
            }


# Sentinel used to mark rows that only indicate the absence of CNVs for a
# patient. These rows must not take part in recurrency calculations and keep
# "NA" in every column besides Patient/Panel/Design.
NA_VALUE = "NA"


def _build_na_row(patient, design, panel=None):
    row = {
        "Patient": patient,
        "Design": design,
        "Type": NA_VALUE,
        "Region": NA_VALUE,
        "Exon count": NA_VALUE,
        "BF": NA_VALUE,
        "Reads.expected": NA_VALUE,
        "Reads.observed": NA_VALUE,
        "Reads.ratio": NA_VALUE,
        "Exons": NA_VALUE,
        "Recurrency_run": NA_VALUE,
        "Recurrency_HC-CNV": NA_VALUE,
        "HC-CNV.reads_ratio.[0-0.2[": NA_VALUE,
        "HC-CNV.reads_ratio.[0.2-1[": NA_VALUE,
        "HC-CNV.reads_ratio.[1-1.83[": NA_VALUE,
        "HC-CNV.reads_ratio.>=1.83": NA_VALUE,
        "HC-CNV.patient_included": NA_VALUE,
        "HC-CNV.patients_subsample": NA_VALUE,
    }
    if panel is not None:
        row["Panel"] = panel
    return row


def build_output_rows(samples, design_gene_sets, panel_gene_sets):
    design_rows = []
    panel_rows = []

    # Prevent duplicate rows when the same patient/design/cnv appears multiple times.
    seen_design_rows = set()
    seen_panel_rows = set()

    # Track which patients/panels actually got a CNV row, so a placeholder
    # "NA" row can be added for the ones that did not.
    patients_with_design_cnv = set()
    patients_with_panel_cnv = set()

    for patient, patient_info in samples.items():
        vcf = patient_info.get("vcf", "")
        if not vcf:
            continue

        panels = patient_info.get("panels", {})
        design = patient_info["design"]

        print(f"Reading VCF for sample '{patient}': {vcf}")
        for record in iter_vcf_records(vcf):
            record_genes = record["genes"]
            if not record_genes:
                continue

            design_genes = design_gene_sets.get(design, set())
            if design_genes and not record_genes.isdisjoint(design_genes):
                row_key = (patient, design, record["region"])
                if row_key not in seen_design_rows:
                    seen_design_rows.add(row_key)
                    patients_with_design_cnv.add(patient)
                    design_rows.append(
                        {
                            "Patient": patient,
                            "Design": design,
                            "Type": record["type"],
                            "Region": record["region"],
                            "Exon count": record["nexons"],
                            "BF": record["bf"],
                            "Reads.expected": record["reads_expected"],
                            "Reads.observed": record["reads_observed"],
                            "Reads.ratio": record["reads_ratio"],
                            "Exons": record["exons"],
                        }
                    )

            for panel in panels:
                panel_genes = panel_gene_sets.get((design, panel), set())
                if not panel_genes:
                    continue
                if record_genes.isdisjoint(panel_genes):
                    continue

                # Panel output includes panel column: keep one row per patient/panel/cnv.
                row_key = (patient, panel, design, record["region"])
                if row_key in seen_panel_rows:
                    continue
                seen_panel_rows.add(row_key)
                patients_with_panel_cnv.add((patient, panel))

                panel_rows.append(
                    {
                        "Patient": patient,
                        "Panel": panel,
                        "Design": design,
                        "Type": record["type"],
                        "Region": record["region"],
                        "Exon count": record["nexons"],
                        "BF": record["bf"],
                        "Reads.expected": record["reads_expected"],
                        "Reads.observed": record["reads_observed"],
                        "Reads.ratio": record["reads_ratio"],
                        "Exons": record["exons"],
                    }
                )

    # Add a placeholder "NA" row for patients (and patient/panel combinations)
    # for which no CNVs were found at all, so they still appear in the report.
    for patient, patient_info in samples.items():
        if not patient_info.get("vcf"):
            continue
        design = patient_info["design"]
        if patient not in patients_with_design_cnv:
            design_rows.append(_build_na_row(patient, design))
        for panel in patient_info.get("panels", {}):
            if (patient, panel) not in patients_with_panel_cnv:
                panel_rows.append(_build_na_row(patient, design, panel=panel))

    return design_rows, panel_rows


def build_design_patient_counts(samples):
    design_to_patients = defaultdict(set)
    for patient, patient_info in samples.items():
        design_to_patients[patient_info["design"]].add(patient)
    return {design: len(patients) for design, patients in design_to_patients.items()}


def annotate_rows(rows, cnvs_databases, cnvs_databases_count_per_design, design_patient_counts):
    patient_sets_per_cnv = defaultdict(set)
    for row in rows:
        # Placeholder "no CNV found" rows keep their preset "NA" values and
        # must not be counted towards recurrency.
        if row["Region"] == NA_VALUE:
            continue
        # Recurrency is per design and genomic region (chr:start-stop).
        key = (row["Design"], row["Region"])
        patient_sets_per_cnv[key].add(row["Patient"])

    for row in rows:
        if row["Region"] == NA_VALUE:
            continue

        design = row["Design"]
        region = row["Region"]
        patient = row["Patient"]

        recurrency_count = len(patient_sets_per_cnv[(design, region)])
        total_design_patients = design_patient_counts.get(design, 0)
        row["Recurrency_run"] = f"{recurrency_count}/{total_design_patients}"

        design_db = cnvs_databases.get(design, {})
        hc = design_db.get(region)

        if hc:
            row["Recurrency_HC-CNV"] = (
                f"{hc.get('patient_count', '0')}/{hc.get('total_patient_count', '0')}"
            )
            row["HC-CNV.reads_ratio.[0-0.2["] = hc.get("[0-0.2[", "NA")
            row["HC-CNV.reads_ratio.[0.2-1["] = hc.get("[0.2-1[", "NA")
            row["HC-CNV.reads_ratio.[1-1.83["] = hc.get("[1-1.83[", "NA")
            row["HC-CNV.reads_ratio.>=1.83"] = hc.get(">=1.83", "NA")
            row["HC-CNV.patients_subsample"] = hc.get("patients_subsample", "NA")
            row["HC-CNV.patient_included"] = (
                "True" if patient in hc.get("patients_reads_ratio", {}) else "False"
            )
        else:
            total_count = cnvs_databases_count_per_design.get(design, 0)
            row["Recurrency_HC-CNV"] = f"0/{total_count}"
            row["HC-CNV.reads_ratio.[0-0.2["] = "NA"
            row["HC-CNV.reads_ratio.[0.2-1["] = "NA"
            row["HC-CNV.reads_ratio.[1-1.83["] = "NA"
            row["HC-CNV.reads_ratio.>=1.83"] = "NA"
            row["HC-CNV.patients_subsample"] = "NA"
            row["HC-CNV.patient_included"] = "NA"


def write_output(rows, output_path, include_panel=False):
    header = build_output_header(include_panel=include_panel)

    with open(output_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow([row.get(col, "") for col in header])


def prepare_output_dir(path):
    output_dir = os.path.abspath(path)
    if output_dir == os.path.abspath(os.sep):
        raise ValueError(f"Refusing to remove protected output directory: {output_dir}")

    # Writing into the current directory (e.g. a Nextflow work dir): keep it, outputs are overwritten
    if output_dir == os.path.abspath(os.getcwd()):
        return output_dir

    if os.path.lexists(output_dir):
        print(f"WARNING: Output directory '{output_dir}' exists; removing it before writing outputs.")
        if os.path.isdir(output_dir) and not os.path.islink(output_dir):
            shutil.rmtree(output_dir)
        else:
            os.remove(output_dir)

    os.makedirs(output_dir)
    return output_dir


def main():
    args = parse_arguments()
    (
        samples,
        design_gene_sets,
        panel_gene_sets,
        cnvs_databases,
        cnvs_databases_count_per_design,
    ) = read_samplesheet(args.samplesheet)
    if not samples:
        raise ValueError(f"No samples found in samplesheet: {args.samplesheet}")

    design_rows, panel_rows = build_output_rows(
        samples, design_gene_sets, panel_gene_sets
    )
    design_patient_counts = build_design_patient_counts(samples)
    for row in sorted(design_rows, key=lambda r: (r["Design"], r["Patient"])):
        if row["Region"] == NA_VALUE:
            print(f"No CNVs found for patient '{row['Patient']}' in design '{row['Design']}'.")
    for row in sorted(panel_rows, key=lambda r: (r["Design"], r["Panel"], r["Patient"])):
        if row["Region"] == NA_VALUE:
            print(
                f"No CNVs found for patient '{row['Patient']}' in panel '{row['Panel']}' "
                f"of design '{row['Design']}'."
            )
    annotate_rows(
        design_rows,
        cnvs_databases,
        cnvs_databases_count_per_design,
        design_patient_counts,
    )
    annotate_rows(
        panel_rows,
        cnvs_databases,
        cnvs_databases_count_per_design,
        design_patient_counts,
    )

    output_dir = prepare_output_dir(args.output_dir)
    write_output(
        design_rows,
        os.path.join(output_dir, "cnvs_exomedepth_designs_summary.tsv"),
    )
    write_output(
        panel_rows,
        os.path.join(output_dir, "cnvs_exomedepth_panels_summary.tsv"),
        include_panel=True,
    )


if __name__ == "__main__":
    main()
