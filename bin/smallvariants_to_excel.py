#!/usr/bin/env python3
"""
smallvariants_to_excel.py

Processes variant and coverage data for SeqCap runs and outputs text reports
and Excel workbooks in the exact format of smallvariantsToExcel_no_job.pl.

Designed for standalone execution and integration into Nextflow pipelines.
Does NOT make any API calls; all API data is passed as input files.

Inputs:
  1. API data file (api_data_YYYYMMDD_HHMMSS.json) or directory containing the API JSON files.
  2. samplesheet.csv: Comma-separated runinfo file with one panel, panel BED,
      design BED, and transcript mapping file per row.
  3. Directory containing coverage and variant files (with the standard structure:
     <sample>/output*/*.vcf.gz and <sample>/<sample>*/*.per-base.bed.gz).

Outputs:
  - Text files per patient & panel in <output_dir>/<patient>/
  - Consolidated raw data files in <output_dir>/_RAWdata/
  - Excel files in <output_dir>/ with exact styling, colors, and freeze panes:
    - <patient>_<panel>.xlsx
    - <patient>/_<patient>_<design>.xlsx and <patient>_<design>.xlsx when the design is a samplesheet panel
    - _RAWdata_<runName>.xlsx
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import logging
import math
import os
import re
import shutil
import sys
import xml.etree.ElementTree as ET
import zipfile
from collections import defaultdict
from decimal import Decimal, InvalidOperation
from typing import Any

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process nf-cmgg/smallvariants data to produce exact Excel and text reports."
    )
    parser.add_argument(
        "--api-data",
        dest="api_data",
        required=True,
        help="Path to combined api_data.json OR directory containing individual API JSON files.",
    )
    parser.add_argument(
        "--samplesheet",
        dest="samplesheet",
        required=True,
        help="Path to samplesheet.csv linking samples, panels, designs, and BED files.",
    )
    parser.add_argument(
        "--input-dir",
        dest="input_dir",
        required=False,
        default=".",
        help="nf-cmgg/smallvariants directory containing per-sample coverage and VCF files. "
        "Not needed when vcf/coverage paths are provided explicitly in the samplesheet (default: '.').",
    )
    parser.add_argument(
        "--output-dir",
        dest="output_dir",
        default="results",
        help="Directory to write output text files and Excel workbooks (default: ./results).",
    )
    parser.add_argument(
        "--run-name",
        dest="run_name",
        default=None,
        help="SMAPLE_RUN_NAME (e.g. NVQ_671). Inferred from input-dir or samplesheet if omitted.",
    )
    parser.add_argument(
        "--threshold-coverage",
        dest="threshold_coverage",
        type=int,
        default=31,
        help="Threshold for low coverage flag (default: 31).",
    )
    parser.add_argument(
        "--build",
        dest="build",
        default="GRCh38/hg38",
        help="Genome build label for reports (default: GRCh38/hg38).",
    )
    parser.add_argument(
        "--variant-caller",
        dest="variant_caller",
        default="vardict",
        help="Variant caller name to match VCF filename (default: vardict).",
    )
    parser.add_argument(
        "--runtype",
        dest="runtype",
        default="SeqCap",
        help="Runtype label for Excel title styling (default: SeqCap).",
    )
    return parser.parse_args()


# ==============================================================================
# API Data Loader
# ==============================================================================

def load_api_data(api_path: str) -> dict[str, Any]:
    """
    Loads API resources from a single JSON bundle or directory of JSON files.
    Returns dictionary with keys:
      - 'ford_approved_assays': { chrom: { assay: {start, stop} } }
      - 'mdg_approved_assays': { chrom: { assay: {start, stop} } }
      - 'cmgg_variants': { ensg: { gnot: {'class': ..., 'tags': ...} } }
      - 'cmgg_variants_mdg': { ensg: { gnot: class_str } }
    """
    api_data: dict[str, Any] = {
        "ford_approved_assays": {},
        "mdg_approved_assays": {},
        "cmgg_variants": {},
        "cmgg_variants_mdg": {},
    }

    if os.path.isfile(api_path):
        logger.info(f"Loading API data from bundle: {api_path}")
        with open(api_path, "r", encoding="utf-8") as f:
            bundle = json.load(f)
            for k in api_data:
                if k in bundle:
                    api_data[k] = bundle[k]
    elif os.path.isdir(api_path):
        logger.info(f"Loading API data from directory: {api_path}")
        mapping = {
            "ford_approved_assays": "ford_approved_assays.json",
            "mdg_approved_assays": "mdg_approved_assays.json",
            "cmgg_variants": "cmgg_variants.json",
            "cmgg_variants_mdg": "cmgg_variants_mdg.json",
        }
        for key, fname in mapping.items():
            fpath = os.path.join(api_path, fname)
            if os.path.isfile(fpath):
                with open(fpath, "r", encoding="utf-8") as f:
                    api_data[key] = json.load(f)
            else:
                logger.warning(f"API data file {fpath} not found.")
    else:
        raise FileNotFoundError(f"API data path '{api_path}' does not exist.")

    return api_data


# ==============================================================================
# Samplesheet & Input Discovery
# ==============================================================================

class SampleInfo:
    def __init__(self, sample: str, design: str):
        self.sample: str = sample
        self.design: str = design
        self.panel: list[str] = []
        self.design_bed: str = ""
        self.design_genelist: str = ""
        self.panel_bed: dict[str, str] = {}
        self.panel_genelist: dict[str, str] = {}
        self.transcripts_file: dict[str, str] = {}
        self.vcf_path: str = ""
        self.coverage_bed_path: str = ""


def load_samplesheet(samplesheet_path: str) -> dict[str, SampleInfo]:
    """
    Parses comma-separated samplesheet.csv.
    Supported column names:
      - sample / patient / samplename / dna_nr
      - panel (exactly one panel per row)
      - design
      - design_bed
      - design_genelist (path to a gene list file, one gene symbol per line, for the design)
      - panel_bed
      - panel_genelist (path to a gene list file, one gene symbol per line, per panel)
      - transcript_file (required, per design)
      - vcf / vcf_file (optional)
      - coverage / coverage_bed (optional)
    """
    if not os.path.isfile(samplesheet_path):
        raise FileNotFoundError(f"Samplesheet '{samplesheet_path}' not found.")

    samples: dict[str, SampleInfo] = {}

    with open(samplesheet_path, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"Empty or invalid samplesheet: {samplesheet_path}")

        # Normalise column headers
        field_map = {col.strip().lower(): col for col in reader.fieldnames if col}

        sample_col = None
        for candidate in ["sample", "id"]:
            if candidate in field_map:
                sample_col = field_map[candidate]
                break

        if not sample_col:
            raise ValueError(
                f"Could not identify sample column in samplesheet. Found: {list(field_map.keys())}"
            )

        panel_col = None
        for candidate in ["panel"]:
            if candidate in field_map:
                panel_col = field_map[candidate]
                break

        design_col = None
        for candidate in ["design"]:
            if candidate in field_map:
                design_col = field_map[candidate]
                break

        design_bed_col = field_map.get("design_bed")
        design_genelist_col = field_map.get("design_genelist")
        panel_bed_col = field_map.get("panel_bed")
        panel_genelist_col = field_map.get("panel_genelist")
        transcripts_file_col = field_map.get("transcript_file")

        vcf_col = field_map.get("vcf") or field_map.get("vcf_file") or field_map.get("vcf_path")
        cov_col = field_map.get("coverage") or field_map.get("coverage_bed") or field_map.get("mosdepth_bed")

        missing_columns = [
            name
            for name, column in {
                "panel": panel_col,
                "design": design_col,
                "panel_bed": panel_bed_col,
                "design_bed": design_bed_col,
                "panel_genelist": panel_genelist_col,
                "design_genelist": design_genelist_col,
                "transcript_file": transcripts_file_col,
            }.items()
            if column is None
        ]
        if missing_columns:
            raise ValueError(
                f"Samplesheet is missing required column(s): {', '.join(missing_columns)}."
            )

        for row_idx, row in enumerate(reader, start=2):
            if None in row:
                raise ValueError(
                    f"Samplesheet row {row_idx} contains more comma-separated fields than its header."
                )

            sample = (row.get(sample_col) or "").strip()
            if not sample:
                continue

            panel_val = (row.get(panel_col) or "").strip() if panel_col else ""
            design_val = (row.get(design_col) or "").strip() if design_col else ""
            panel_bed_val = (row.get(panel_bed_col) or "").strip() if panel_bed_col else ""
            panel_genelist_val = (
                (row.get(panel_genelist_col) or "").strip() if panel_genelist_col else ""
            )
            design_bed_val = (row.get(design_bed_col) or "").strip() if design_bed_col else ""
            design_genelist_val = (
                (row.get(design_genelist_col) or "").strip() if design_genelist_col else ""
            )
            transcripts_file_val = (
                (row.get(transcripts_file_col) or "").strip() if transcripts_file_col else ""
            )
            vcf_val = (row.get(vcf_col) or "").strip() if vcf_col else ""
            cov_val = (row.get(cov_col) or "").strip() if cov_col else ""

            if (
                not panel_val
                or not design_val
                or not panel_bed_val
                or not design_bed_val
                or not panel_genelist_val
                or not design_genelist_val
                or not transcripts_file_val
            ):
                raise ValueError(
                    f"Samplesheet row {row_idx} must contain sample, panel, design, "
                    "panel_bed, design_bed, panel_genelist, design_genelist, and transcript_file values."
                )
            if re.search(r"[,;]", panel_val):
                raise ValueError(
                    f"Samplesheet row {row_idx} contains multiple panels; provide one panel per row."
                )

            if sample not in samples:
                sample_obj = SampleInfo(sample=sample, design=design_val)
                samples[sample] = sample_obj
            else:
                sample_obj = samples[sample]
                if sample_obj.design != design_val:
                    raise ValueError(
                        f"Sample '{sample}' uses multiple designs ({sample_obj.design} and {design_val}); "
                        "split the analysis into separate samplesheets."
                    )

            if sample_obj.design_bed and sample_obj.design_bed != design_bed_val:
                raise ValueError(
                    f"Sample '{sample}' uses multiple design_bed files for design '{design_val}'."
                )
            sample_obj.design_bed = design_bed_val
            if sample_obj.design_genelist and sample_obj.design_genelist != design_genelist_val:
                raise ValueError(
                    f"Sample '{sample}' uses multiple design_genelist files for design '{design_val}'."
                )
            sample_obj.design_genelist = design_genelist_val

            if vcf_val and not sample_obj.vcf_path:
                sample_obj.vcf_path = vcf_val

            if cov_val and not sample_obj.coverage_bed_path:
                sample_obj.coverage_bed_path = cov_val

            if panel_val not in sample_obj.panel:
                sample_obj.panel.append(panel_val)
            if panel_val in sample_obj.panel_bed and sample_obj.panel_bed[panel_val] != panel_bed_val:
                raise ValueError(
                    f"Sample '{sample}' and panel '{panel_val}' have multiple panel_bed files."
                )
            sample_obj.panel_bed[panel_val] = panel_bed_val
            if (
                panel_val in sample_obj.panel_genelist
                and sample_obj.panel_genelist[panel_val] != panel_genelist_val
            ):
                raise ValueError(
                    f"Sample '{sample}' and panel '{panel_val}' have multiple panel_genelist files."
                )
            sample_obj.panel_genelist[panel_val] = panel_genelist_val

            existing_file = sample_obj.transcripts_file.get(design_val)
            if existing_file and existing_file != transcripts_file_val:
                raise ValueError(
                    f"Design '{design_val}' uses multiple transcript mapping files."
                )
            sample_obj.transcripts_file[design_val] = transcripts_file_val

    logger.info(f"Loaded {len(samples)} samples from samplesheet.")
    return samples


def load_transcripts_mapping(mapping_path: str) -> dict[str, dict[str, str]]:
    """
    Loads one design-specific transcript mapping file (gene<TAB>gene_id<TAB>ENST<TAB>NM).
    Returns ENST -> {"gene": gene_name, "gene_id": ensembl_gene_id, "nm": NM}.
    """
    if not os.path.isfile(mapping_path):
        raise FileNotFoundError(f"Transcript mapping file '{mapping_path}' not found.")

    mapping: dict[str, dict[str, str]] = {}
    with open(mapping_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                mapping[parts[2]] = {"gene": parts[0], "gene_id": parts[1], "nm": parts[3]}
    return mapping


def discover_sample_files(
    sample_obj: SampleInfo, input_dir: str, variant_caller: str = "vardict"
) -> None:
    """
    Locates VCF and Mosdepth per-base bed.gz files for a sample if not explicitly provided in samplesheet.
    Standard directory structure:
      <input_dir>/<sample>/output*/*.<variant_caller>.vcf.gz (or *.vcf.gz)
      <input_dir>/<sample>/<sample>*/*.per-base.bed.gz (or *.per-base.bed.gz)
    """
    sample = sample_obj.sample
    sample_dir = os.path.join(input_dir, sample)

    # 1. Discover VCF
    if not sample_obj.vcf_path:
        found_vcf = None
        # Check standard output subdirectory
        if os.path.isdir(sample_dir):
            for entry in os.listdir(sample_dir):
                sub_path = os.path.join(sample_dir, entry)
                if os.path.isdir(sub_path) and entry.startswith("output"):
                    for fname in os.listdir(sub_path):
                        if fname.endswith(f".{variant_caller}.vcf.gz"):
                            found_vcf = os.path.join(sub_path, fname)
                            break
                        elif fname.endswith(".vcf.gz") and not found_vcf:
                            found_vcf = os.path.join(sub_path, fname)

            # Fallback to direct search in sample_dir
            if not found_vcf:
                for fname in os.listdir(sample_dir):
                    if fname.endswith(f".{variant_caller}.vcf.gz"):
                        found_vcf = os.path.join(sample_dir, fname)
                        break
                    elif fname.endswith(".vcf.gz") and not found_vcf:
                        found_vcf = os.path.join(sample_dir, fname)

        # Fallback to search in input_dir root
        if not found_vcf:
            for fname in os.listdir(input_dir):
                if fname.startswith(sample) and fname.endswith(f".{variant_caller}.vcf.gz"):
                    found_vcf = os.path.join(input_dir, fname)
                    break

        if found_vcf:
            sample_obj.vcf_path = found_vcf
        else:
            logger.warning(f"No VCF file found for sample '{sample}' in {input_dir}.")

    # 2. Discover Mosdepth Coverage file
    if not sample_obj.coverage_bed_path:
        found_cov = None
        if os.path.isdir(sample_dir):
            for entry in os.listdir(sample_dir):
                sub_path = os.path.join(sample_dir, entry)
                if os.path.isdir(sub_path) and entry.startswith(sample):
                    for fname in os.listdir(sub_path):
                        if fname.endswith(".per-base.bed.gz"):
                            found_cov = os.path.join(sub_path, fname)
                            break
                elif os.path.isdir(sub_path):
                    for fname in os.listdir(sub_path):
                        if fname.endswith(".per-base.bed.gz"):
                            found_cov = os.path.join(sub_path, fname)
                            break

            # Fallback direct search in sample_dir
            if not found_cov:
                for fname in os.listdir(sample_dir):
                    if fname.endswith(".per-base.bed.gz"):
                        found_cov = os.path.join(sample_dir, fname)
                        break

        # Fallback to search in input_dir root
        if not found_cov:
            for fname in os.listdir(input_dir):
                if fname.startswith(sample) and fname.endswith(".per-base.bed.gz"):
                    found_cov = os.path.join(input_dir, fname)
                    break

        if found_cov:
            sample_obj.coverage_bed_path = found_cov
        else:
            raise FileNotFoundError(
                f"ERROR: No mosdepth per-base.bed.gz coverage file found for sample '{sample}' in {input_dir}."
            )


# ==============================================================================
# BED Parsing & Coverage Statistics Calculation
# ==============================================================================

class BedRegion:
    def __init__(self, chrom: str, start: int, end: int, attribute: str):
        self.chrom: str = chrom
        self.start: int = start
        self.end: int = end
        self.attribute: str = attribute
        self.gene: str = ""
        self.exon_number: int = 0
        self.gene_or_assay: str = ""
        self._parse_attribute()

    def _parse_attribute(self) -> None:
        parts = self.attribute.split(";")
        self.gene_or_assay = parts[0] if len(parts) > 0 else ""
        gene = self.gene_or_assay
        if "---" in gene:
            gene = gene.split("---")[0]
        if "_Exon" in gene:
            gene = gene.split("_Exon")[0]
        self.gene = gene

        exon_str = parts[4] if len(parts) > 4 else "0"
        try:
            self.exon_number = int(exon_str)
        except ValueError:
            self.exon_number = 0


def load_bed_regions(bed_path: str) -> list[BedRegion]:
    """Reads target BED file (chrom, start, end, attribute)."""
    if not os.path.isfile(bed_path):
        raise FileNotFoundError(f"Target BED file '{bed_path}' not found.")

    regions: list[BedRegion] = []
    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            attr = cols[3] if len(cols) > 3 else f"Region_{start}_{end}"
            regions.append(BedRegion(chrom, start, end, attr))
    return regions


def load_gene_list(genelist_path: str) -> set[str]:
    """Loads panel gene symbols from a samplesheet-supplied file, one gene name per line."""
    if not os.path.isfile(genelist_path):
        raise FileNotFoundError(f"Panel gene list file '{genelist_path}' not found.")
    genes: set[str] = set()
    with open(genelist_path, "r", encoding="utf-8-sig") as f:
        for line in f:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.add(gene)
    return genes


def calculate_coverage_statistics(
    coverage_bed_gz: str,
    target_regions: list[BedRegion],
    genome_build: str = "hg38",
    design: str = "",
) -> list[list[str]]:
    """
    Pure Python, highly efficient streaming intersection between mosdepth per-base.bed.gz
    and target BED regions. Calculates:
      build, chrom, start, end, attribute, length, min, max, mean, median, stdev,
      zero_coverage_bases, proportion_covered, length_above_10X, %_above_10X,
      length_above_20X, %_above_20X, length_above_30X, %_above_30X.

    Matches seqcap_calculate_coverage_statistics_per_bed_region.pl exactly.
    """
    # Group target regions by chromosome
    regions_by_chr: dict[str, list[BedRegion]] = defaultdict(list)
    for r in target_regions:
        chrom_key = r.chrom if r.chrom.startswith("chr") else f"chr{r.chrom}"
        regions_by_chr[chrom_key].append(r)

    # Sort target regions per chromosome by start position
    for chrom_key in regions_by_chr:
        regions_by_chr[chrom_key].sort(key=lambda x: (x.start, x.end))

    # Track depth intervals per region: index -> list of (depth, length)
    region_intervals: dict[BedRegion, list[tuple[int, int]]] = {r: [] for r in target_regions}

    # Stream through mosdepth per-base bed.gz
    with gzip.open(coverage_bed_gz, "rt", encoding="utf-8") as f:
        current_chrom: str | None = None
        active_regions: list[BedRegion] = []
        target_idx = 0
        chrom_targets: list[BedRegion] = []

        for line in f:
            if not line:
                continue
            cols = line.strip().split("\t")
            if len(cols) < 4:
                continue

            m_chr = cols[0]
            if not m_chr.startswith("chr"):
                m_chr = f"chr{m_chr}"

            if m_chr != current_chrom:
                current_chrom = m_chr
                chrom_targets = regions_by_chr.get(current_chrom, [])
                target_idx = 0
                active_regions = []

            if not chrom_targets:
                continue

            m_start = int(cols[1])
            m_end = int(cols[2])
            m_depth = int(cols[3])

            # Add newly reachable target regions to active list
            while target_idx < len(chrom_targets) and chrom_targets[target_idx].start < m_end:
                active_regions.append(chrom_targets[target_idx])
                target_idx += 1

            # Remove target regions that end before or at m_start
            active_regions = [r for r in active_regions if r.end > m_start]

            # Accumulate overlap for active target regions
            for r in active_regions:
                overlap_start = max(r.start, m_start)
                overlap_end = min(r.end, m_end)
                if overlap_end > overlap_start:
                    length = overlap_end - overlap_start
                    region_intervals[r].append((m_depth, length))

    # Compute statistics for each region
    output_rows: list[list[str]] = []

    for r in target_regions:
        target_len = r.end - r.start
        intervals = region_intervals[r]
        covered_len = sum(l for _, l in intervals)

        # Account for gaps in mosdepth coverage as depth 0
        if covered_len < target_len:
            intervals.append((0, target_len - covered_len))

        count_bases = target_len
        if count_bases == 0:
            continue

        min_cov = min(d for d, _ in intervals)
        max_cov = max(d for d, _ in intervals)
        sum_cov = sum(d * l for d, l in intervals)
        mean_cov = sum_cov / count_bases
        mean_str = f"{mean_cov:.2f}".replace(".", ",")

        variance = sum(l * ((d - mean_cov) ** 2) for d, l in intervals) / count_bases
        stdev = math.sqrt(variance)
        stdev_str = f"{stdev:.2f}".replace(".", ",")

        # Median calculation matching Perl
        depths: list[int] = []
        for d, l in intervals:
            depths.extend([d] * l)
        depths.sort()
        n = len(depths)
        if n % 2 == 0:
            med = (depths[n // 2 - 1] + depths[n // 2]) / 2
        else:
            med = depths[n // 2]
        median_str = str(int(med)) if med == int(med) else str(med).replace(".", ",")

        zero_count = sum(l for d, l in intervals if d == 0)
        prop_cov = (count_bases - zero_count) / count_bases * 100
        prop_cov_str = f"{prop_cov:.2f}".replace(".", ",")

        count_10 = sum(l for d, l in intervals if d >= 10)
        pct_10_str = f"{(count_10 / count_bases * 100):.2f}".replace(".", ",")

        count_20 = sum(l for d, l in intervals if d >= 20)
        pct_20_str = f"{(count_20 / count_bases * 100):.2f}".replace(".", ",")

        count_30 = sum(l for d, l in intervals if d >= 30)
        pct_30_str = f"{(count_30 / count_bases * 100):.2f}".replace(".", ",")

        row = [
            genome_build,
            r.chrom,
            str(r.start),
            str(r.end),
            r.attribute,
            str(count_bases),
            str(min_cov),
            str(max_cov),
            mean_str,
            median_str,
            stdev_str,
            str(zero_count),
            prop_cov_str,
            str(count_10),
            pct_10_str,
            str(count_20),
            pct_20_str,
            str(count_30),
            pct_30_str,
        ]
        output_rows.append(row)

    # Sort design coverage lines alphabetically by gene, numerically by exon number, alphabetically by assay
    def sort_key(row_cols: list[str]) -> tuple[str, int, str]:
        attr = row_cols[4]
        parts = attr.split(";")
        g_assay = parts[0] if len(parts) > 0 else ""
        g = g_assay
        if "---" in g:
            g = g.split("---")[0]
        if "_Exon" in g:
            g = g.split("_Exon")[0]
        exon_num = 0
        if len(parts) > 4:
            try:
                exon_num = int(parts[4])
            except ValueError:
                exon_num = 0
        return (g, exon_num, g_assay)

    if not design.startswith(("Solid", "Hemato")):
        output_rows.sort(key=sort_key)
    return output_rows


COVERAGE_HEADER = [
    "#build",
    "chromosome",
    "start",
    "end",
    "attribute",
    "length",
    "min",
    "max",
    "mean",
    "median",
    "stdev",
    "zero_coverage_bases",
    "proportion_covered",
    "length_above_10X",
    "%_above_10X",
    "length_above_20X",
    "%_above_20X",
    "length_above_30X",
    "%_above_30X",
]


# ==============================================================================
# Variant Processing & Recurrency
# ==============================================================================

VARIANT_ORDER_CLASSES = [
    "unknown",
    "CLASS 5",
    "CLASS 4",
    "CLASS 3",
    "CLASS 2",
    "CLASS 1",
    "KNOWN FALSE POSITIVE",
]

VARIANT_HEADERS_OUTPUT = [
    "Build",
    "Patient",
    "VCF_g_nomenclature",
    "HGVS_g_nomenclature",
    "HGVS_Coding_region_change",
    "HGVS_Amino_acid_change",
    "Recurrency",
    "Target",
    "Reads_ref",
    "Reads_alt",
    "Coverage",
    "Variant_allele_frequency",
    "Qual_score",
    "CMGG_result",
    "DNA_tags",
    "MDG_class",
    "Ford_assay",
    "Consequence",
    "Impact",
    "ClinVar",
    "dbSNP",
    "dbSNP_Alternate_IDs",
    "Strand",
    "VEP_Codons",
    "VEP_EXON",
    "VEP_INTRON",
    "VEP_cDNA_position",
    "MAX_AF (1000G, ESP and gnomAD)",
    "MAX_AF_POPS (1000G, ESP and gnomAD)",
    "AF_gnomAD",
    "gnomAD_Allele_Count",
    "gnomAD_Allele_Number",
    "gnomAD_Homozygous_Alleles",
    "EOG_AF",
    "REVEL_score",
    "CADD_score",
    "BayesDel_score (incl MaxAF)",
    "BayesDel_score (excl MaxAF)",
    "SpliceAI_score_AG (max distance 50bp)",
    "SpliceAI_score_AL (max distance 50bp)",
    "SpliceAI_score_DG (max distance 50bp)",
    "SpliceAI_score_DL (max distance 50bp)",
    "SpliceAI_position_AG (max distance 50bp)",
    "SpliceAI_position_AL (max distance 50bp)",
    "SpliceAI_position_DG (max distance 50bp)",
    "SpliceAI_position_DL (max distance 50bp)",
    "SpliceRegion",
    "MaxEntScan_alt",
    "MaxEntScan_diff",
    "MaxEntScan_ref",
    "rf_score",
    "ada_score",
    "Protein_id",
    "Interpro_domain",
    "Transcript_ENST",
    "Transcript_RefSeq",
    "Gene_id",
    "PubMed",
]


class VariantRecord:
    def __init__(self) -> None:
        self.data: dict[str, str] = {}
        self.chrom_num: int = 0
        self.pos: int = 0
        self.duplicate: int = 1
        self.class_category: str = "unknown"
        self.g_nomenclature: str = ""


def normalize_chrom_to_num(chrom: str) -> int:
    c = chrom.replace("chr", "").strip()
    if c in ("X", "x"):
        return 23
    if c in ("Y", "y"):
        return 24
    if c in ("M", "MT", "m", "mt"):
        return 25
    try:
        return int(c)
    except ValueError:
        return 99


def parse_vcf_variants(
    vcf_path: str,
    target_regions: list[BedRegion],
    patient: str,
    design: str,
    api_data: dict[str, Any],
    transcripts_map: dict[str, dict[str, dict[str, str]]],
    build: str = "GRCh38/hg38",
) -> list[VariantRecord]:
    """
    Parses VCF, filters variants against target regions, annotates against VEP CSQ,
    CMGG, and MDG databases, and calculates VAF.
    """
    if not vcf_path or not os.path.isfile(vcf_path):
        return []

    # Build target interval index per chromosome
    target_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for r in target_regions:
        c = r.chrom.replace("chr", "").strip()
        target_intervals[c].append((r.start, r.end))

    # Merge overlapping intervals for fast lookup
    merged_intervals: dict[str, list[tuple[int, int]]] = {}
    for c, intervals in target_intervals.items():
        intervals.sort(key=lambda x: x[0])
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append((start, end))
            else:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        merged_intervals[c] = merged

    def is_in_target(chrom_str: str, pos: int, ref_len: int = 1) -> bool:
        # BED start is 0-based; convert to 1-based to match VCF POS. BED end is
        # already the correct 1-based inclusive end. Overlap is checked against
        # the variant's full REF span (pos .. pos+ref_len-1), matching bcftools
        # view -R semantics, so indels overlapping only at their tail are kept.
        c = chrom_str.replace("chr", "").strip()
        intervals = merged_intervals.get(c, [])
        variant_end = pos + ref_len - 1
        for start, end in intervals:
            region_start = start + 1
            if region_start > variant_end:
                break
            if end >= pos:
                return True
        return False

    records: list[VariantRecord] = []
    csq_fields: list[str] = []
    feature_idx = -1

    pos_dup_counter: dict[tuple[int, int], int] = defaultdict(int)

    open_func = gzip.open if vcf_path.endswith(".gz") else open
    with open_func(vcf_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("##INFO=<ID=CSQ,"):
                match = re.search(r'Format:\s*([^">]+)', line)
                if match:
                    csq_fields = match.group(1).split("|")
                    if "Feature" in csq_fields:
                        feature_idx = csq_fields.index("Feature")
                continue

            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 10:
                continue

            raw_chr, raw_pos, rs_id, ref, alt, qual, filter_val, info, fmt, sample_val = cols[:10]
            if not alt:
                continue
            alt = re.sub(r",?<NON_REF>$", "", alt)
            if not alt:
                continue

            pos = int(raw_pos)
            chrom_num = normalize_chrom_to_num(raw_chr)
            clean_chr = raw_chr.replace("chr", "").strip()

            # Filter variant to target BED
            if not is_in_target(clean_chr, pos, len(ref)):
                continue

            pos_dup_counter[(chrom_num, pos)] += 1
            duplicate = pos_dup_counter[(chrom_num, pos)]

            # g_nomenclature
            chr_label = "X" if chrom_num == 23 else "Y" if chrom_num == 24 else str(chrom_num)
            if len(ref) == 1:
                g_nomenclature = f"{chr_label}_{pos}_{pos}_{ref}/{alt}"
            else:
                end_pos = pos + len(ref) - 1
                g_nomenclature = f"{chr_label}_{pos}_{end_pos}_{ref}/{alt}"

            # Overlapping assays
            chrom_key = str(chrom_num)
            ford_assays = []
            for assay, coords in api_data.get("ford_approved_assays", {}).get(chrom_key, {}).items():
                if coords["stop"] > pos and coords["start"] < pos:
                    ford_assays.append(assay)
            overlapping_ford = ",".join(sorted(ford_assays))

            mdg_assays = []
            for assay, coords in api_data.get("mdg_approved_assays", {}).get(chrom_key, {}).items():
                if coords["stop"] > pos and coords["start"] < pos:
                    mdg_assays.append(assay)
            overlapping_mdg = ",".join(sorted(mdg_assays))

            # Parse INFO tags
            info_dict: dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_dict[k] = v

            # Parse FORMAT & SAMPLE
            fmt_keys = fmt.split(":")
            sample_vals = sample_val.split(":")
            sample_dict = dict(zip(fmt_keys, sample_vals))

            # Allele depths
            ad_val = sample_dict.get("AD", "")
            ad_parts = ad_val.split(",") if ad_val else []
            ref_reads = ad_parts[0] if len(ad_parts) > 0 else "0"
            alt_reads = ",".join(ad_parts[1:]) if len(ad_parts) > 1 else "0"
            coverage_dp = sample_dict.get("DP", "0")

            # Calculate VAF
            vaf_str = "0"
            try:
                dp_num = float(coverage_dp)
                if dp_num > 0 and ad_parts:
                    if len(ad_parts) > 2:
                        freqs = [f"{(float(cnt) / dp_num * 100):.2f}" for cnt in ad_parts[1:]]
                        vaf_str = ",".join(freqs)
                    elif len(ad_parts) == 2:
                        vaf_str = f"{(float(ad_parts[1]) / dp_num * 100):.2f}"
            except (ValueError, ZeroDivisionError):
                vaf_str = "0"

            # Parse VEP CSQ
            target_gene = ""
            enst_id = ""
            hgvsc = ""
            hgvsp = ""
            hgvsg = ""
            consequence = ""
            impact = ""
            strand = ""
            codons = ""
            exon_str = ""
            intron_str = ""
            cdna_pos = ""
            max_af = ""
            max_af_pops = ""
            gnomad_af = ""
            eog_af = ""
            ensp = ""
            ensg_id = ""
            nm_id = "."
            pubmed = ""
            existing_var = ""
            splice_ai = {}

            # Only transcripts present in the design's transcripts_file are valid;
            # this disambiguates variants annotated against overlapping genes (matches modules/variants.pm).
            design_transcripts = transcripts_map.get(design, {})

            csq_str = info_dict.get("CSQ", "")
            if csq_str and csq_fields:
                transcripts_list = csq_str.split(",")
                for tr_entry in transcripts_list:
                    tr_fields = tr_entry.split("|")
                    tr_id = tr_fields[feature_idx] if feature_idx < len(tr_fields) else ""

                    transcript_info = design_transcripts.get(tr_id)
                    if not transcript_info:
                        continue

                    c_dict = dict(zip(csq_fields, tr_fields))
                    cand_hgvsc = c_dict.get("HGVSc", "")

                    target_gene = transcript_info["gene"]
                    ensg_id = transcript_info["gene_id"]
                    nm_id = transcript_info["nm"] or "."
                    enst_id = tr_id
                    hgvsc = cand_hgvsc
                    hgvsp = c_dict.get("HGVSp", "")
                    hgvsg = c_dict.get("HGVSg", "")
                    consequence = c_dict.get("Consequence", "")
                    impact = c_dict.get("IMPACT", "")
                    strand = c_dict.get("STRAND", "")
                    codons = c_dict.get("Codons", "")
                    exon_str = c_dict.get("EXON", "")
                    intron_str = c_dict.get("INTRON", "")
                    cdna_pos = c_dict.get("cDNA_position", "")
                    max_af = c_dict.get("MAX_AF", "")
                    max_af_pops = c_dict.get("MAX_AF_POPS", "")
                    gnomad_af = c_dict.get("gnomAD_AF", "")
                    eog_af = c_dict.get("EOG_AF", "")
                    ensp = c_dict.get("ENSP", "")
                    pubmed = c_dict.get("PUBMED", "").replace("&", ",")
                    existing_var = c_dict.get("Existing_variation", "").replace("&", ",")

                    for k in ["DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]:
                        splice_ai[k] = c_dict.get(f"SpliceAI_pred_{k}", "") or "."

                    # Exit on first transcript with a valid coding change
                    if cand_hgvsc:
                        break

            # Format HGVS strings
            if hgvsp:
                hgvsp = hgvsp.replace("%3D", "=")
                if ":" in hgvsp:
                    hgvsp = hgvsp.split(":", 1)[1]
                if hgvsp.startswith("p."):
                    hgvsp = f"p.({hgvsp[2:]})"
            elif hgvsc:
                hgvsp = "p.?"

            if ":" in hgvsc:
                hgvsc = hgvsc.split(":", 1)[1]

            # Look up variant classification in CMGG and MDG
            cmgg_class = ""
            dna_tags = ""
            mdg_class = ""

            if ensg_id and hgvsg:
                cmgg_info = api_data.get("cmgg_variants", {}).get(ensg_id, {}).get(hgvsg, {})
                cmgg_class = cmgg_info.get("class", "").rstrip(" - ")
                dna_tags = cmgg_info.get("tags", "")

                mdg_class = api_data.get("cmgg_variants_mdg", {}).get(ensg_id, {}).get(hgvsg, "").rstrip(" - ")

            # REVEL & Interpro
            revel = info_dict.get("REVEL_score", ".")
            if revel != ".":
                rev_unique = list(dict.fromkeys(revel.split(",")))
                revel = next((x for x in rev_unique if x != "."), ".")

            interpro = info_dict.get("Interpro_domain", "")
            if interpro:
                interpro = ",".join(dict.fromkeys(interpro.split(",")))

            # Clinvar
            clinvar = info_dict.get("clinvar_sig", ".").replace("&", ",")

            # Determine classification category for sorting
            category = "unknown"
            if cmgg_class:
                for cls_name in ["CLASS 5", "CLASS 4", "CLASS 3", "CLASS 2", "CLASS 1", "KNOWN FALSE POSITIVE"]:
                    if cmgg_class.startswith(cls_name):
                        category = cls_name
                        break

            # Build record dictionary
            rec = VariantRecord()
            rec.chrom_num = chrom_num
            rec.pos = pos
            rec.duplicate = duplicate
            rec.class_category = category
            rec.g_nomenclature = g_nomenclature

            # Populate fields
            d = rec.data
            d["Build"] = build
            d["Patient"] = patient
            d["VCF_g_nomenclature"] = g_nomenclature
            d["HGVS_g_nomenclature"] = hgvsg
            d["HGVS_Coding_region_change"] = hgvsc or "."
            d["HGVS_Amino_acid_change"] = hgvsp or "."
            d["Recurrency"] = ""  # Populated later
            d["Target"] = target_gene or "."
            d["Reads_ref"] = ref_reads
            d["Reads_alt"] = alt_reads
            d["Coverage"] = coverage_dp
            d["Variant_allele_frequency"] = vaf_str
            d["Qual_score"] = qual
            d["CMGG_result"] = cmgg_class or "."
            d["DNA_tags"] = dna_tags or "."
            d["MDG_class"] = mdg_class or "."
            d["Ford_assay"] = overlapping_ford or "."
            d["Consequence"] = consequence or "."
            d["Impact"] = impact or "."
            d["ClinVar"] = clinvar or "."
            d["dbSNP"] = rs_id if rs_id and rs_id != "." else "."
            d["dbSNP_Alternate_IDs"] = existing_var or "."
            d["Strand"] = strand or "."
            d["VEP_Codons"] = codons or "."
            d["VEP_EXON"] = exon_str or "."
            d["VEP_INTRON"] = intron_str or "."
            d["VEP_cDNA_position"] = cdna_pos or "."
            d["MAX_AF (1000G, ESP and gnomAD)"] = max_af or "."
            d["MAX_AF_POPS (1000G, ESP and gnomAD)"] = max_af_pops or "."
            d["AF_gnomAD"] = gnomad_af or "."
            d["gnomAD_Allele_Count"] = info_dict.get("gnomAD_AC", ".")
            d["gnomAD_Allele_Number"] = info_dict.get("gnomAD_AN", ".")
            d["gnomAD_Homozygous_Alleles"] = info_dict.get("gnomAD_Hom", ".")
            d["EOG_AF"] = eog_af or "."
            d["REVEL_score"] = revel
            d["CADD_score"] = info_dict.get("CADD_phred", ".")
            d["BayesDel_score (incl MaxAF)"] = info_dict.get("BayesDel_addAF_score", ".")
            d["BayesDel_score (excl MaxAF)"] = info_dict.get("BayesDel_noAF_score", ".")
            d["SpliceAI_score_AG (max distance 50bp)"] = splice_ai.get("DS_AG", ".")
            d["SpliceAI_score_AL (max distance 50bp)"] = splice_ai.get("DS_AL", ".")
            d["SpliceAI_score_DG (max distance 50bp)"] = splice_ai.get("DS_DG", ".")
            d["SpliceAI_score_DL (max distance 50bp)"] = splice_ai.get("DS_DL", ".")
            d["SpliceAI_position_AG (max distance 50bp)"] = splice_ai.get("DP_AG", ".")
            d["SpliceAI_position_AL (max distance 50bp)"] = splice_ai.get("DP_AL", ".")
            d["SpliceAI_position_DG (max distance 50bp)"] = splice_ai.get("DP_DG", ".")
            d["SpliceAI_position_DL (max distance 50bp)"] = splice_ai.get("DP_DL", ".")
            d["SpliceRegion"] = info_dict.get("SpliceRegion", ".")
            d["MaxEntScan_alt"] = info_dict.get("MaxEntScan_alt", ".")
            d["MaxEntScan_diff"] = info_dict.get("MaxEntScan_diff", ".")
            d["MaxEntScan_ref"] = info_dict.get("MaxEntScan_ref", ".")
            d["rf_score"] = info_dict.get("rf_score", ".")
            d["ada_score"] = info_dict.get("ada_score", ".")
            d["Protein_id"] = ensp or "."
            d["Interpro_domain"] = interpro or "."
            d["Transcript_ENST"] = enst_id or "."
            d["Transcript_RefSeq"] = nm_id
            d["Gene_id"] = ensg_id or "."
            d["PubMed"] = pubmed or "."

            # Format European decimal commas
            non_comma_cols = {
                "rs_number",
                "dbSNP",
                "Transcript_RefSeq",
                "HGVS_Coding_region_change",
                "HGVS_Amino_acid_change",
                "VCF_g_nomenclature",
                "HGVS_g_nomenclature",
                "CMGG_result",
            }
            for col_name, val in d.items():
                if val == "null":
                    d[col_name] = "0"
                elif col_name == "Qual_score" and val != ".":
                    try:
                        qual_value = Decimal(val)
                    except InvalidOperation:
                        d[col_name] = val.replace(".", ",")
                    else:
                        if qual_value.is_finite() and qual_value == qual_value.to_integral_value():
                            d[col_name] = str(int(qual_value))
                        else:
                            d[col_name] = val.replace(".", ",")
                elif col_name not in non_comma_cols and val != ".":
                    d[col_name] = val.replace(".", ",")

            records.append(rec)

    # Sort variants according to class hierarchy then genomic coordinates
    class_priority = {cls_name: idx for idx, cls_name in enumerate(VARIANT_ORDER_CLASSES)}

    def variant_sort_key(r: VariantRecord) -> tuple[int, int, int, int]:
        prio = class_priority.get(r.class_category, 0)
        return (prio, r.chrom_num, r.pos, r.duplicate)

    records.sort(key=variant_sort_key)
    return records


def calculate_recurrency(
    patient_variants_by_design: dict[str, dict[str, list[VariantRecord]]],
    total_patients_per_design: dict[str, int],
) -> dict[str, dict[str, str]]:
    """
    Computes recurrency for each variant per design across the entire run.
    Returns: recurrency_map[design][g_nomenclature] = f"{count}/{total_patients}"
    """
    recurrency_map: dict[str, dict[str, str]] = defaultdict(dict)
    for design, patient_dict in patient_variants_by_design.items():
        total_design_pts = total_patients_per_design.get(design, len(patient_dict))
        variant_patient_counts: dict[str, int] = defaultdict(int)

        for patient, records in patient_dict.items():
            seen_gnots: set[str] = set()
            for r in records:
                if r.g_nomenclature and r.g_nomenclature not in seen_gnots:
                    seen_gnots.add(r.g_nomenclature)
                    variant_patient_counts[r.g_nomenclature] += 1

        for gnot, cnt in variant_patient_counts.items():
            recurrency_map[design][gnot] = f"{cnt}/{total_design_pts}"

        for records in patient_dict.values():
            for r in records:
                r.data["Recurrency"] = recurrency_map[design].get(r.g_nomenclature, f"1/{total_design_pts}")

    return recurrency_map


# ==============================================================================
# Excel Writer (Native XML)
# ==============================================================================

class UniversalExcelWriter:
    """
    Writes Excel workbooks (.xlsx) with exact colors, fonts, freeze panes,
    and column widths as produced by Excel::Writer::XLSX in the Perl script.
    """

    def __init__(self, filename: str):
        self.filename: str = filename
        self.sheets: dict[str, list[list[Any]]] = {}
        self.sheet_col_widths: dict[str, dict[int, float]] = {}
        self.sheet_cell_formats: dict[str, dict[tuple[int, int], dict[str, Any]]] = {}
        self.sheet_freeze_panes: dict[str, int] = {}  # 0-indexed row below which to freeze

    def add_worksheet(self, name: str) -> None:
        if name not in self.sheets:
            self.sheets[name] = []
            self.sheet_col_widths[name] = {}
            self.sheet_cell_formats[name] = {}

    def set_column(self, sheet_name: str, col_range: str, width: float) -> None:
        """Sets width for column range like 'A:A' or 'B:D'."""
        self.add_worksheet(sheet_name)
        parts = col_range.split(":")
        start_letter = parts[0].strip().upper()
        end_letter = parts[1].strip().upper() if len(parts) > 1 else start_letter

        def letter_to_idx(col_str: str) -> int:
            idx = 0
            for char in col_str:
                idx = idx * 26 + (ord(char) - ord("A") + 1)
            return idx - 1

        s_idx = letter_to_idx(start_letter)
        e_idx = letter_to_idx(end_letter)
        for c in range(s_idx, e_idx + 1):
            self.sheet_col_widths[sheet_name][c] = width

    def freeze_panes(self, sheet_name: str, row_idx: int) -> None:
        """Freezes pane below row_idx (0-indexed)."""
        self.sheet_freeze_panes[sheet_name] = row_idx

    def write(
        self,
        sheet_name: str,
        row: int,
        col: int,
        value: Any,
        fmt: dict[str, Any] | None = None,
    ) -> None:
        self.add_worksheet(sheet_name)
        data = self.sheets[sheet_name]
        while len(data) <= row:
            data.append([])
        row_list = data[row]
        while len(row_list) <= col:
            row_list.append("")
        row_list[col] = value

        if fmt:
            self.sheet_cell_formats[sheet_name][(row, col)] = fmt

    def save(self) -> None:
        os.makedirs(os.path.dirname(os.path.abspath(self.filename)), exist_ok=True)
        self._save_native_xml()

    def _save_native_xml(self) -> None:
        """
        Pure Python fallback for creating standard .xlsx files using zipfile and xml,
        without requiring any external packages.
        """
        import html

        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
            # 1. [Content_Types].xml
            ct = [
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
                '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
                '<Default Extension="xml" ContentType="application/xml"/>',
                '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
                '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>',
            ]
            for i in range(1, len(self.sheets) + 1):
                ct.append(
                    f'<Override PartName="/xl/worksheets/sheet{i}.xml" '
                    'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                )
            ct.append("</Types>")
            zf.writestr("[Content_Types].xml", "\n".join(ct))

            # 2. _rels/.rels
            zf.writestr(
                "_rels/.rels",
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">\n'
                '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>\n'
                "</Relationships>",
            )

            # 3. xl/_rels/workbook.xml.rels
            wb_rels = [
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
                '<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>',
            ]
            for i in range(1, len(self.sheets) + 1):
                wb_rels.append(
                    f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i}.xml"/>'
                )
            wb_rels.append("</Relationships>")
            zf.writestr("xl/_rels/workbook.xml.rels", "\n".join(wb_rels))

            # 4. xl/styles.xml with custom colors (#66D9FF, #B3EDFF, #D8D8D8, #F78181)
            styles_xml = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
  <fonts count="4">
    <font><sz val="11"/><name val="Calibri"/></font>
    <font><b/><sz val="18"/><color rgb="FF000000"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><color rgb="FF000000"/><name val="Calibri"/></font>
    <font><b/><sz val="11"/><name val="Calibri"/></font>
  </fonts>
  <fills count="6">
    <fill><patternFill patternType="none"/></fill>
    <fill><patternFill patternType="gray125"/></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FF66D9FF"/><bgColor rgb="FF66D9FF"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFB3EDFF"/><bgColor rgb="FFB3EDFF"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFD8D8D8"/><bgColor rgb="FFD8D8D8"/></patternFill></fill>
    <fill><patternFill patternType="solid"><fgColor rgb="FFF78181"/><bgColor rgb="FFF78181"/></patternFill></fill>
  </fills>
  <borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>
  <cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
  <cellXfs count="6">
    <xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
    <xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>
    <xf numFmtId="0" fontId="2" fillId="3" borderId="0" xfId="0" applyFont="1" applyFill="1"/>
    <xf numFmtId="0" fontId="3" fillId="4" borderId="0" xfId="0" applyFont="1" applyFill="1"/>
    <xf numFmtId="0" fontId="0" fillId="5" borderId="0" xfId="0" applyFill="1"/>
    <xf numFmtId="0" fontId="3" fillId="0" borderId="0" xfId="0" applyFont="1"/>
  </cellXfs>
</styleSheet>"""
            zf.writestr("xl/styles.xml", styles_xml)

            # 5. xl/workbook.xml
            wb_xml = [
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
                "<sheets>",
            ]
            for idx, sheet_name in enumerate(self.sheets.keys(), start=1):
                clean_name = html.escape(sheet_name)
                wb_xml.append(f'<sheet name="{clean_name}" sheetId="{idx}" r:id="rId{idx}"/>')
            wb_xml.append("</sheets></workbook>")
            zf.writestr("xl/workbook.xml", "\n".join(wb_xml))

            # 6. xl/worksheets/sheetX.xml
            for s_idx, (sheet_name, rows) in enumerate(self.sheets.items(), start=1):
                ws_xml = [
                    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
                    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
                ]

                # Freeze panes
                if sheet_name in self.sheet_freeze_panes:
                    f_row = self.sheet_freeze_panes[sheet_name]
                    ws_xml.append(
                        f'<sheetViews><sheetView tabSelected="1" workbookViewId="0">'
                        f'<pane ySplit="{f_row}" topLeftCell="A{f_row+1}" activePane="bottomLeft" state="frozen"/>'
                        f"</sheetView></sheetViews>"
                    )

                # Column widths
                widths = self.sheet_col_widths.get(sheet_name, {})
                if widths:
                    ws_xml.append("<cols>")
                    for col_idx, width in sorted(widths.items()):
                        c_num = col_idx + 1
                        ws_xml.append(
                            f'<col min="{c_num}" max="{c_num}" width="{width}" customWidth="1"/>'
                        )
                    ws_xml.append("</cols>")

                ws_xml.append("<sheetData>")
                for r_idx, row in enumerate(rows):
                    r_num = r_idx + 1
                    ws_xml.append(f'<row r="{r_num}">')
                    for c_idx, val in enumerate(row):
                        # Calculate cell reference
                        c_temp = c_idx
                        c_letters = ""
                        while c_temp >= 0:
                            c_letters = chr(ord("A") + (c_temp % 26)) + c_letters
                            c_temp = c_temp // 26 - 1
                        cell_ref = f"{c_letters}{r_num}"

                        # Determine style index
                        fmt = self.sheet_cell_formats.get(sheet_name, {}).get((r_idx, c_idx), {})
                        s_style = 0
                        if fmt.get("bg_color") == "#66d9ff" or fmt.get("bg_color") == "#66D9FF":
                            s_style = 1
                        elif fmt.get("bg_color") == "#b3edff" or fmt.get("bg_color") == "#B3EDFF":
                            s_style = 2
                        elif fmt.get("bg_color") == "#D8D8D8":
                            s_style = 3
                        elif fmt.get("bg_color") == "#F78181":
                            s_style = 4
                        elif fmt.get("bold"):
                            s_style = 5

                        val_str = html.escape(str(val))
                        ws_xml.append(
                            f'<c r="{cell_ref}" s="{s_style}" t="inlineStr"><is><t>{val_str}</t></is></c>'
                        )
                    ws_xml.append("</row>")
                ws_xml.append("</sheetData></worksheet>")
                zf.writestr(f"xl/worksheets/sheet{s_idx}.xml", "\n".join(ws_xml))

        with open(self.filename, "wb") as f_out:
            f_out.write(buf.getvalue())


# ==============================================================================
# Excel Report Generation (writeExcel & writeExcel_rawdata)
# ==============================================================================

def write_patient_excel(
    workbook_path: str,
    patient: str,
    panel: str,
    design: str,
    run_name: str,
    design_gene_count: int,
    panel_gene_count: int,
    coverage_rows: list[list[str]],
    variant_records: list[VariantRecord],
    threshold_coverage: int,
    api_data: dict[str, Any],
    build: str = "GRCh38/hg38",
    runtype: str = "SeqCap",
) -> None:
    """Generates patient/panel Excel workbook with 'info', 'variants', and 'coverage' tabs."""
    writer = UniversalExcelWriter(workbook_path)

    # 1. INFO TAB
    writer.add_worksheet("info")
    writer.set_column("info", "A:A", 65)
    writer.set_column("info", "B:B", 15)
    writer.set_column("info", "C:C", 12)
    writer.set_column("info", "D:D", 15)
    writer.set_column("info", "E:F", 25)
    writer.set_column("info", "G:K", 15)
    writer.set_column("info", "H:I", 25)

    fmt_banner = {"bold": True, "size": 18, "color": "black", "bg_color": "#66d9ff"}
    fmt_header_info = {"bold": True, "size": 11, "color": "black", "bg_color": "#b3edff"}
    fmt_table_header = {"bold": True, "bg_color": "#D8D8D8"}

    writer.write("info", 0, 0, "Patient: .................", fmt_banner)
    writer.write("info", 1, 0, f"DNA-nr: {patient}", fmt_banner)
    writer.write("info", 2, 0, "Geboortedatum: .................", fmt_banner)
    writer.write("info", 3, 0, "", fmt_banner)
    writer.write("info", 4, 0, f"SeqCap design: {design} ({design_gene_count} genen)", fmt_banner)
    gen_text = "gen" if panel_gene_count == 1 else "genen"
    writer.write("info", 5, 0, f"Screening panel: {panel} ({panel_gene_count} {gen_text})", fmt_banner)
    writer.write("info", 6, 0, f"Sequencing run: {run_name}", fmt_banner)

    # Low coverage section
    curr_row = 9
    writer.write(
        "info",
        curr_row,
        0,
        f"Regio's met een te lage coverage (minimum coverage <{threshold_coverage} X)",
        fmt_header_info,
    )
    curr_row += 1

    low_cov_rows: list[tuple[str, int, tuple[str, str, str, str, str]]] = []
    for r in coverage_rows:
        # r[6] is min coverage
        try:
            min_c = int(float(r[6]))
        except ValueError:
            min_c = 0
        if min_c < threshold_coverage:
            attr = r[4]
            parts = attr.split(";")
            g = parts[0]
            exon_num = parts[4] if len(parts) > 4 else "0"
            exon_attr = f"{g}_exon{exon_num}"
            try:
                exon_sort_num = int(exon_num)
            except ValueError:
                exon_sort_num = 0
            mean_c = r[8]
            zero_c = r[11]

            # Overlapping Ford assays with at least 40bp overlap at ends
            chr_t = r[1].replace("chr", "").strip()
            chr_key = "23" if chr_t in ("X", "x") else "24" if chr_t in ("Y", "y") else chr_t
            start_t = int(r[2])
            stop_t = int(r[3])
            overlapping_ford = []
            for assay, coords in api_data.get("ford_approved_assays", {}).get(chr_key, {}).items():
                if coords["stop"] >= (start_t + 20) and coords["start"] <= (stop_t - 20):
                    overlapping_ford.append(assay)
            ford_str = ",".join(sorted(overlapping_ford)) if overlapping_ford else "NA"

            low_cov_rows.append(
                (g, exon_sort_num, (exon_attr, str(min_c), mean_c, zero_c, ford_str))
            )

    if low_cov_rows:
        low_cov_rows.sort(key=lambda row: (row[0], row[1]))
        low_cov_headers = ["Exon", "Minimum", "Mean", "Zero_coverage_bases", "Ford_assay"]
        for c_idx, h in enumerate(low_cov_headers):
            writer.write("info", curr_row, c_idx, h, fmt_table_header)
        curr_row += 1
        for _, _, row_data in low_cov_rows:
            for c_idx, val in enumerate(row_data):
                writer.write("info", curr_row, c_idx, val)
            curr_row += 1
    else:
        writer.write("info", curr_row, 0, "Alle regio's hebben voldoende coverage")
        curr_row += 1

    curr_row += 1

    # Variant summary tables per class
    variant_info_headers = [
        "Build",
        "DNA_tags",
        "Transcript",
        "HGVS_Coding_region_change",
        "HGVS_Amino_acid_change",
        "Coverage",
        "Variant_allele_frequency",
        "Ford_assay",
        "Recurrency",
        "Target",
    ]

    classes_to_report = [
        ("CLASS 5", "Klasse 5 varianten", "Geen gekende klasse 5 varianten gedetecteerd"),
        ("CLASS 4", "Klasse 4 varianten", "Geen gekende klasse 4 varianten gedetecteerd"),
        ("CLASS 3", "Klasse 3 varianten", "Geen gekende klasse 3 varianten gedetecteerd"),
        ("unknown", "Varianten zonder gekende klasse", "Geen varianten gedetecteerd zonder klasse"),
    ]

    for cls_key, section_title, empty_msg in classes_to_report:
        writer.write("info", curr_row, 0, section_title, fmt_header_info)
        curr_row += 1
        cls_records = [r for r in variant_records if r.class_category == cls_key]
        if cls_records:
            cls_records.sort(key=lambda record: (record.data.get("Target", "."), record.g_nomenclature))
            for c_idx, h in enumerate(variant_info_headers):
                writer.write("info", curr_row, c_idx, h, fmt_table_header)
            curr_row += 1
            for rec in cls_records:
                d = rec.data
                row_vals = [
                    d.get("Build", build),
                    d.get("DNA_tags", "."),
                    d.get("Transcript_RefSeq", "."),
                    d.get("HGVS_Coding_region_change", "."),
                    d.get("HGVS_Amino_acid_change", "."),
                    d.get("Coverage", "0"),
                    d.get("Variant_allele_frequency", "0"),
                    d.get("Ford_assay", "."),
                    d.get("Recurrency", "."),
                    d.get("Target", "."),
                ]
                for c_idx, val in enumerate(row_vals):
                    writer.write("info", curr_row, c_idx, val)
                curr_row += 1
        else:
            writer.write("info", curr_row, 0, empty_msg)
            curr_row += 1
        curr_row += 1

    # 2. VARIANTS TAB
    writer.add_worksheet("variants")
    writer.set_column("variants", "A:A", 30)
    writer.set_column("variants", "B:B", 15)
    writer.set_column("variants", "C:F", 30)
    writer.set_column("variants", "G:L", 12)
    writer.set_column("variants", "M:M", 20)
    writer.set_column("variants", "N:N", 10)
    writer.set_column("variants", "O:BH", 25)
    writer.set_column("variants", "Y:Y", 35)
    writer.set_column("variants", "AE:AF", 30)
    writer.set_column("variants", "AP:AW", 30)
    writer.set_column("variants", "BI:BI", 255)

    writer.write("variants", 0, 0, f"variants {patient}", fmt_banner)
    writer.freeze_panes("variants", 2)

    for c_idx, h in enumerate(VARIANT_HEADERS_OUTPUT):
        writer.write("variants", 1, c_idx, h, fmt_table_header)

    v_row = 2
    if variant_records:
        for rec in variant_records:
            d = rec.data
            for c_idx, h in enumerate(VARIANT_HEADERS_OUTPUT):
                writer.write("variants", v_row, c_idx, d.get(h, "."))
            v_row += 1
    else:
        writer.write("variants", v_row, 0, "Geen varianten gevonden")

    # 3. COVERAGE TAB
    writer.add_worksheet("coverage")
    writer.set_column("coverage", "A:A", 30)
    writer.set_column("coverage", "B:E", 20)
    writer.write("coverage", 0, 0, f"coverage {patient}", fmt_banner)
    writer.freeze_panes("coverage", 2)

    cov_headers = ["Exon", "Minimum", "Mean", "Zero_coverage_bases"]
    for c_idx, h in enumerate(cov_headers):
        writer.write("coverage", 1, c_idx, h, fmt_table_header)

    fmt_cov_low = {"bg_color": "#F78181"}
    c_row = 2
    for r in coverage_rows:
        attr = r[4]
        parts = attr.split(";")
        g = parts[0]
        exon_num = parts[4] if len(parts) > 4 else "0"
        exon_attr = f"{g}_exon{exon_num}"
        min_c_val = r[6]
        mean_c_val = r[8]
        zero_c_val = r[11]

        writer.write("coverage", c_row, 0, exon_attr)

        try:
            is_low = int(float(min_c_val)) < threshold_coverage
        except ValueError:
            is_low = False

        if is_low:
            writer.write("coverage", c_row, 1, min_c_val, fmt_cov_low)
        else:
            writer.write("coverage", c_row, 1, min_c_val)

        writer.write("coverage", c_row, 2, mean_c_val)
        writer.write("coverage", c_row, 3, zero_c_val)
        c_row += 1

    writer.save()


def write_rawdata_excel(
    workbook_path: str,
    run_name: str,
    variants_raw_path: str,
    coverage_raw_path: str,
) -> None:
    """Generates consolidated _RAWdata_<runName>.xlsx containing variants and coverage sheets."""
    writer = UniversalExcelWriter(workbook_path)
    fmt_table_header = {"bold": True, "bg_color": "#D8D8D8"}

    # Variants tab
    if os.path.isfile(variants_raw_path):
        writer.add_worksheet("variants")
        writer.set_column("variants", "A:A", 10)
        writer.set_column("variants", "B:D", 15)
        writer.set_column("variants", "E:H", 30)
        writer.set_column("variants", "I:N", 12)
        writer.set_column("variants", "O:O", 20)
        writer.set_column("variants", "P:P", 10)
        writer.set_column("variants", "Q:BJ", 25)
        writer.set_column("variants", "AA:AA", 35)
        writer.set_column("variants", "AR:AY", 30)
        writer.set_column("variants", "AG:AH", 30)
        writer.set_column("variants", "BK:BK", 255)
        writer.freeze_panes("variants", 1)

        with open(variants_raw_path, "r", encoding="utf-8") as f:
            for r_idx, line in enumerate(f):
                cols = line.strip("\r\n").split("\t")
                for c_idx, val in enumerate(cols):
                    if r_idx == 0:
                        writer.write("variants", r_idx, c_idx, val, fmt_table_header)
                    else:
                        writer.write("variants", r_idx, c_idx, val)

    # Coverage tab
    if os.path.isfile(coverage_raw_path):
        writer.add_worksheet("coverage")
        writer.set_column("coverage", "A:V", 18)
        writer.set_column("coverage", "H:H", 50)
        writer.freeze_panes("coverage", 1)

        with open(coverage_raw_path, "r", encoding="utf-8") as f:
            for r_idx, line in enumerate(f):
                cols = line.strip("\r\n").split("\t")
                for c_idx, val in enumerate(cols):
                    if r_idx == 0:
                        writer.write("coverage", r_idx, c_idx, val, fmt_table_header)
                    else:
                        writer.write("coverage", r_idx, c_idx, val)

    writer.save()


# ==============================================================================
# Main Pipeline Workflow
# ==============================================================================

def main() -> None:
    args = parse_arguments()

    run_name = args.run_name
    if not run_name:
        base_name = os.path.splitext(os.path.basename(args.samplesheet))[0]
        run_name = base_name.replace(".runinfo", "").replace(".samplesheet", "")
    logger.info(f"Starting smallvariants processing for run: {run_name}")

    # 1. Load API Data
    api_data = load_api_data(args.api_data)

    # 2. Load Samplesheet
    samples_map = load_samplesheet(args.samplesheet)
    if not samples_map:
        logger.error(f"No samples found in samplesheet {args.samplesheet}.")
        sys.exit(1)

    # 3. Discover Files for Each Sample
    for sample_obj in samples_map.values():
        discover_sample_files(sample_obj, args.input_dir, args.variant_caller)

    # 4. Load design-specific transcript mappings
    transcripts_map: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)
    transcript_files_by_design: dict[str, str] = {}
    for s_obj in samples_map.values():
        for design, mapping_path in s_obj.transcripts_file.items():
            existing_path = transcript_files_by_design.get(design)
            if existing_path and existing_path != mapping_path:
                raise ValueError(f"Design '{design}' uses multiple transcript mapping files.")
            transcript_files_by_design[design] = mapping_path

    for design, mapping_path in transcript_files_by_design.items():
        logger.info(f"Loading transcript mapping for design {design} from {mapping_path}...")
        transcripts_map[design] = load_transcripts_mapping(mapping_path)

    # Create Output Directories
    if os.path.isdir(args.output_dir):
        logger.warning(f"Output directory '{args.output_dir}' already exists; removing it before regenerating.")
        shutil.rmtree(args.output_dir)
    rawdata_dir = os.path.join(args.output_dir, "_RAWdata")
    os.makedirs(rawdata_dir, exist_ok=True)
    raw_variants_path = os.path.join(rawdata_dir, f"{run_name}_variants.txt")
    raw_coverage_path = os.path.join(rawdata_dir, f"{run_name}_coverage.txt")

    # Track distinct patients per design
    total_patients_per_design: dict[str, int] = defaultdict(int)
    for s_obj in samples_map.values():
        total_patients_per_design[s_obj.design] += 1

    # Container for variants per design to compute recurrency
    # design -> patient -> list of VariantRecord
    design_patient_variants: dict[str, dict[str, list[VariantRecord]]] = defaultdict(lambda: defaultdict(list))

    # Containers for generated data per sample
    # sample -> { 'design_coverage': [...], 'panels_coverage': {panel: [...]}, 'design_variants': [...], 'panel_variants': {panel: [...]} }
    processed_samples_data: dict[str, Any] = {}

    # Cache loaded BED files: path -> list of BedRegion
    bed_cache: dict[str, list[BedRegion]] = {}

    def get_bed(bpath: str) -> list[BedRegion]:
        if bpath not in bed_cache:
            bed_cache[bpath] = load_bed_regions(bpath)
        return bed_cache[bpath]

    # Cache loaded panel gene lists: path -> set of gene symbols
    genelist_cache: dict[str, set[str]] = {}

    def get_panel_genes(genelist_path: str) -> set[str]:
        if genelist_path not in genelist_cache:
            genelist_cache[genelist_path] = load_gene_list(genelist_path)
        return genelist_cache[genelist_path]

    # Process each patient's coverage and design variants
    sample_items = sorted(samples_map.items())
    total_samples = len(sample_items)

    logger.info(f"Processing {total_samples} samples...")
    for idx, (sample, s_obj) in enumerate(sample_items, start=1):
        logger.info(f"[{idx}/{total_samples}] Processing sample: {sample} (Design: {s_obj.design})")
        patient_out_dir = os.path.join(args.output_dir, sample)
        os.makedirs(patient_out_dir, exist_ok=True)

        design = s_obj.design
        design_bed_path = s_obj.design_bed
        if not design_bed_path or not os.path.isfile(design_bed_path):
            raise FileNotFoundError(
                f"Design BED file '{design_bed_path}' for sample {sample} is missing."
            )
        design_regions = get_bed(design_bed_path)

        # 1. Coverage Calculation for Full Design
        logger.info(f"  -> Calculating coverage stats for design {design}...")
        design_cov_rows = calculate_coverage_statistics(
            s_obj.coverage_bed_path, design_regions, genome_build="hg38", design=design
        )

        # Write design coverage text file
        des_cov_file = os.path.join(patient_out_dir, f"{sample}_{design}_coverage.txt")
        with open(des_cov_file, "w", encoding="utf-8") as f_cov:
            f_cov.write("\t".join(COVERAGE_HEADER) + "\n")
            for r in design_cov_rows:
                f_cov.write("\t".join(r) + "\n")

        # 2. Variants for Full Design
        logger.info(f"  -> Processing variants for design {design}...")
        design_variant_recs = parse_vcf_variants(
            s_obj.vcf_path,
            design_regions,
            patient=sample,
            design=design,
            api_data=api_data,
            transcripts_map=transcripts_map,
            build=args.build,
        )
        design_patient_variants[design][sample] = design_variant_recs

        # 3. Process Subpanels
        panels_cov_rows: dict[str, list[list[str]]] = {}
        panels_variant_recs: dict[str, list[VariantRecord]] = {}

        for panel in s_obj.panel:
            panel_bed_path = s_obj.panel_bed[panel]
            if panel == design and panel_bed_path == design_bed_path:
                panels_cov_rows[panel] = design_cov_rows
                panels_variant_recs[panel] = design_variant_recs
                continue

            panel_regions = get_bed(panel_bed_path)
            panel_genes = get_panel_genes(s_obj.panel_genelist[panel])

            # Filter coverage rows for panel genes
            p_cov_rows: list[list[str]] = []
            for r in design_cov_rows:
                attr = r[4]
                parts = attr.split(";")
                gene = parts[0]
                if "---" in gene:
                    gene = gene.split("---")[0]
                if "_Exon" in gene:
                    gene = gene.split("_Exon")[0]
                if gene in panel_genes:
                    p_cov_rows.append(r)
            panels_cov_rows[panel] = p_cov_rows

            # Write panel coverage text file
            panel_cov_file = os.path.join(patient_out_dir, f"{sample}_{panel}_coverage.txt")
            with open(panel_cov_file, "w", encoding="utf-8") as f_p_cov:
                f_p_cov.write("\t".join(COVERAGE_HEADER) + "\n")
                for r in p_cov_rows:
                    f_p_cov.write("\t".join(r) + "\n")

            # Parse variants for panel
            p_vars = parse_vcf_variants(
                s_obj.vcf_path,
                panel_regions,
                patient=sample,
                design=design,
                api_data=api_data,
                transcripts_map=transcripts_map,
                build=args.build,
            )
            panels_variant_recs[panel] = p_vars

        processed_samples_data[sample] = {
            "design_cov_rows": design_cov_rows,
            "design_variant_recs": design_variant_recs,
            "panels_cov_rows": panels_cov_rows,
            "panels_variant_recs": panels_variant_recs,
        }

    # 5. Compute Recurrency across all samples on run
    logger.info("Computing recurrency metrics across run...")
    recurrency_map = calculate_recurrency(design_patient_variants, total_patients_per_design)

    # Apply recurrency to panel variant records
    for sample, s_data in processed_samples_data.items():
        design = samples_map[sample].design
        total_pts = total_patients_per_design.get(design, len(samples_map))
        for panel, recs in s_data["panels_variant_recs"].items():
            for r in recs:
                r.data["Recurrency"] = recurrency_map.get(design, {}).get(
                    r.g_nomenclature, f"1/{total_pts}"
                )

    # 6. Write Variant Text Files & RAWdata files
    logger.info("Writing finalized variant files and raw data summaries...")
    raw_var_header = ["Panel"] + VARIANT_HEADERS_OUTPUT
    raw_cov_header = ["Patient", "Panel"] + COVERAGE_HEADER

    with open(raw_variants_path, "w", encoding="utf-8") as f_raw_var, open(
        raw_coverage_path, "w", encoding="utf-8"
    ) as f_raw_cov:
        f_raw_var.write("\t".join(raw_var_header) + "\n")
        f_raw_cov.write("\t".join(raw_cov_header) + "\n")

        for sample, s_obj in sample_items:
            patient_out_dir = os.path.join(args.output_dir, sample)
            s_data = processed_samples_data[sample]
            design = s_obj.design

            # Write design variants text file
            des_var_file = os.path.join(patient_out_dir, f"{sample}_{design}_variants.txt")
            with open(des_var_file, "w", encoding="utf-8") as f_des_var:
                if s_data["design_variant_recs"]:
                    f_des_var.write("\t".join(VARIANT_HEADERS_OUTPUT) + "\n")
                    for rec in s_data["design_variant_recs"]:
                        row_vals = [rec.data.get(h, ".") for h in VARIANT_HEADERS_OUTPUT]
                        f_des_var.write("\t".join(row_vals) + "\n")
                else:
                    f_des_var.write("Geen varianten gevonden\n")

            # Write panels variants & raw data
            for panel in s_obj.panel:
                p_vars = s_data["panels_variant_recs"][panel]
                p_covs = s_data["panels_cov_rows"][panel]

                panel_var_file = os.path.join(patient_out_dir, f"{sample}_{panel}_variants.txt")
                with open(panel_var_file, "w", encoding="utf-8") as f_p_var:
                    if p_vars:
                        f_p_var.write("\t".join(VARIANT_HEADERS_OUTPUT) + "\n")
                        for rec in p_vars:
                            row_vals = [rec.data.get(h, ".") for h in VARIANT_HEADERS_OUTPUT]
                            f_p_var.write("\t".join(row_vals) + "\n")
                            # Add to raw data variants
                            f_raw_var.write(panel + "\t" + "\t".join(row_vals) + "\n")
                    else:
                        f_p_var.write("Geen varianten gevonden\n")

                # Add to raw data coverage
                for r in p_covs:
                    f_raw_cov.write(sample + "\t" + panel + "\t" + "\t".join(r) + "\n")

    # 7. Generate Excel Workbooks
    logger.info("Generating Excel workbooks with exact formatting...")
    for sample, s_obj in sample_items:
        patient_out_dir = os.path.join(args.output_dir, sample)
        s_data = processed_samples_data[sample]
        design = s_obj.design

        design_genes = get_panel_genes(s_obj.design_genelist)
        design_gene_count = len(design_genes)

        for panel in s_obj.panel:
            p_covs = s_data["panels_cov_rows"][panel]
            p_vars = s_data["panels_variant_recs"][panel]

            panel_bed_path = s_obj.panel_bed[panel]
            panel_genes = get_panel_genes(s_obj.panel_genelist[panel])
            panel_gene_count = len(panel_genes)

            # Keep the design workbook in the patient folder and also in the
            # results folder, matching smallvariantsToExcel_no_job.pl.
            excel_paths = [
                os.path.join(patient_out_dir, f"_{sample}_{panel}.xlsx")
                if panel == design
                else os.path.join(args.output_dir, f"{sample}_{panel}.xlsx")
            ]
            if panel == design:
                excel_paths.append(os.path.join(args.output_dir, f"{sample}_{panel}.xlsx"))

            for excel_path in excel_paths:
                write_patient_excel(
                    workbook_path=excel_path,
                    patient=sample,
                    panel=panel,
                    design=design,
                    run_name=run_name,
                    design_gene_count=design_gene_count,
                    panel_gene_count=panel_gene_count,
                    coverage_rows=p_covs,
                    variant_records=p_vars,
                    threshold_coverage=args.threshold_coverage,
                    api_data=api_data,
                    build=args.build,
                    runtype=args.runtype,
                )

    # 8. Generate Consolidated RAWdata Excel
    rawdata_excel_path = os.path.join(args.output_dir, f"_RAWdata_{run_name}.xlsx")
    logger.info(f"Writing RAWdata Excel workbook: {rawdata_excel_path}...")
    write_rawdata_excel(
        workbook_path=rawdata_excel_path,
        run_name=run_name,
        variants_raw_path=raw_variants_path,
        coverage_raw_path=raw_coverage_path,
    )

    logger.info("smallvariants_to_excel analysis completed successfully!")


if __name__ == "__main__":
    main()
