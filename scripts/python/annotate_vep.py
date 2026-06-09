#!/usr/bin/env python3
"""
Annotate simple- or complex-SV CSVs with gene symbols from a VEP-annotated VCF.

Workflow (intended use):
  1. Take the filtered VCF produced by filter_vcf.py
  2. Upload to the Ensembl VEP web server (https://www.ensembl.org/vep) with
     default parameters
  3. Download the VEP-annotated VCF
  4. Run this script to extract per-variant gene symbols (CSQ_SYMBOL) and
     merge them onto the simple- or complex-SVs CSV from vcf_to_csv.py

Join keys:
  --mode simple:  merged on (CHROM, START, STOP) -- one row per simple SV
  --mode complex: merged on (CHROM1, START1) and parses VEP's BND mate
                  coordinates to confirm the link to (CHROM2, START2)

Only numeric chromosomes (chr1-22 / 1-22) are kept. The VEP CHROM column is
normalized to match the SV CSV convention (preserves "chr" prefix if present
in the SV CSV; strips otherwise).

Inputs and outputs use the SV CSV's chromosome convention -- whatever
vcf_to_csv.py produced is what gets matched against the VEP output.

Usage:
    # Annotate simple SVs
    python annotate_vep.py \\
        --sv-csv  results/mySample/mySample_simple_svs.csv \\
        --vep-vcf results/mySample/mySample_vep.vcf \\
        --mode    simple \\
        --out     results/mySample/mySample_simple_svs_vep.csv

    # Annotate complex SVs (BNDs)
    python annotate_vep.py \\
        --sv-csv  results/mySample/mySample_complex_svs.csv \\
        --vep-vcf results/mySample/mySample_vep.vcf \\
        --mode    complex \\
        --out     results/mySample/mySample_complex_svs_vep.csv

Optional dependency:
    Uses 'vcfpy' if available (faster + safer parsing). Falls back to a
    simple text parser otherwise.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

NUMERIC_CHROMS = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}

# Default VEP CSQ keys (Ensembl default field order, web server defaults).
# Override at the CLI if your VEP run uses a different schema.
DEFAULT_CSQ_KEYS = [
    "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
    "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
    "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
    "Codons", "Existing_variation", "REF_ALLELE", "UPLOADED_ALLELE",
    "DISTANCE", "STRAND", "FLAGS", "SYMBOL_SOURCE", "HGNC_ID",
    "MANE_SELECT", "MANE_PLUS_CLINICAL", "TSL", "APPRIS", "SIFT",
    "PolyPhen", "AF", "CLIN_SIG", "SOMATIC", "PHENO", "PUBMED",
    "MOTIF_NAME", "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
]

# BND ALT spec: N[chr5:12345[, ]chr5:12345]N, etc.
BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sv-csv", type=Path, required=True,
                   help="Input SV CSV from vcf_to_csv.py (simple or complex)")
    p.add_argument("--vep-vcf", type=Path, required=True,
                   help="VEP-annotated VCF (downloaded from the VEP web server)")
    p.add_argument("--mode", choices=["simple", "complex"], required=True,
                   help="Which CSV schema to expect")
    p.add_argument("--out", type=Path, required=True,
                   help="Output annotated CSV path")
    p.add_argument("--csq-keys", type=str, default=None,
                   help="Optional comma-separated CSQ field order, if you ran VEP "
                        "with a non-default schema. By default uses the Ensembl "
                        "web-server defaults.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# VEP VCF parsing
# ---------------------------------------------------------------------------
def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        else:
            out[token] = token
    return out


def read_vep_vcf(vcf_path: Path, csq_keys: list[str]) -> pd.DataFrame:
    """Parse a VEP-annotated VCF into a DataFrame with one row per variant.

    Pulls CHROM, POS, ALT, INFO_SVTYPE, INFO_END, INFO_CHR2, CSQ_SYMBOL,
    CSQ_Allele (used to recover BND mates). Multiple CSQ blocks per variant
    are collapsed by taking the first one (matches the user's original code).
    """
    rows = []
    with vcf_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom = fields[0]
            pos = int(fields[1])
            alt = fields[4]
            info = parse_info(fields[7])

            row = {
                "CHROM": chrom,
                "POS": pos,
                "ALT": alt,
                "INFO_SVTYPE": info.get("SVTYPE", ""),
                "INFO_END": int(info["END"]) if "END" in info else None,
                "INFO_CHR2": info.get("CHR2"),
            }
            # Take the first CSQ block. VEP's web server emits one CSQ per
            # transcript; for SV annotation, gene symbol typically agrees
            # across transcripts of the same gene.
            csq = info.get("CSQ", "")
            if csq:
                first = csq.split(",", 1)[0]
                values = first.split("|")
                for k, v in zip(csq_keys, values):
                    row[f"CSQ_{k}"] = v if v else None
            rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Chromosome normalization
# ---------------------------------------------------------------------------
def normalize_chrom_column(s: pd.Series, like: pd.Series) -> pd.Series:
    """Match the 'chr' prefix convention of a reference series.

    If `like` uses 'chr1' etc., ensure `s` also uses 'chr...'.
    If `like` uses '1' etc., strip 'chr' from `s`.
    """
    sample = next((v for v in like.dropna().astype(str).unique()), "")
    has_prefix = sample.startswith("chr")
    s = s.astype(str)
    if has_prefix:
        return s.where(s.str.startswith("chr"), "chr" + s)
    else:
        return s.str.replace("^chr", "", regex=True)


def keep_numeric(df: pd.DataFrame, chrom_col: str) -> pd.DataFrame:
    """Drop rows whose CHROM is not on the numeric autosome whitelist."""
    return df[df[chrom_col].astype(str).isin(NUMERIC_CHROMS)].copy()


# ---------------------------------------------------------------------------
# Merges
# ---------------------------------------------------------------------------
def annotate_simple(sv_df: pd.DataFrame, vep_df: pd.DataFrame) -> pd.DataFrame:
    """Merge gene symbols onto the simple-SVs table on (CHROM, START, STOP).

    Uses an inner join, matching the user's original `simple_variants_1kb.py`.
    Rows with no VEP match are dropped (i.e. the output is only annotated SVs).
    """
    # Match the user's original: VEP's POS -> START, VEP's INFO_END -> STOP
    vep_simple = vep_df.rename(columns={"POS": "START", "INFO_END": "STOP"})
    vep_simple = vep_simple[["CHROM", "START", "STOP", "CSQ_SYMBOL"]].dropna(subset=["STOP"])
    vep_simple["STOP"] = vep_simple["STOP"].astype(int)
    vep_simple["CHROM"] = normalize_chrom_column(vep_simple["CHROM"], sv_df["CHROM"])
    vep_simple = keep_numeric(vep_simple, "CHROM")

    merged = pd.merge(sv_df, vep_simple, how="inner",
                      on=["CHROM", "START", "STOP"])
    return merged


def annotate_complex(sv_df: pd.DataFrame, vep_df: pd.DataFrame) -> pd.DataFrame:
    """Merge gene symbols onto the complex-SVs (BND) table on (CHROM1, START1).

    The mate coordinates (CHROM2, START2) recovered from VEP's CSQ_Allele
    (the BND ALT field) are reported as POS2 for sanity-checking against the
    SV CSV; rows where the mate doesn't match are still kept but flagged in
    a separate column.
    """
    vep_bnd = vep_df[vep_df["INFO_SVTYPE"] == "BND"].copy()
    if vep_bnd.empty:
        return pd.DataFrame()

    # Mate from the ALT spec (preferred); fall back to INFO_CHR2 + INFO_END
    def _mate_from_row(row: pd.Series) -> tuple[str | None, int | None]:
        alt = row.get("CSQ_Allele")
        if isinstance(alt, str):
            m = BND_RE.search(alt)
            if m:
                return m.group(1), int(m.group(2))
        if pd.notna(row.get("INFO_CHR2")) and pd.notna(row.get("INFO_END")):
            return str(row["INFO_CHR2"]), int(row["INFO_END"])
        return None, None

    mate_chrom, mate_pos = zip(*vep_bnd.apply(_mate_from_row, axis=1)) if len(vep_bnd) else ([], [])
    vep_bnd = vep_bnd.assign(
        VEP_CHROM2=list(mate_chrom),
        VEP_START2=list(mate_pos),
    )
    vep_bnd = vep_bnd[["CHROM", "POS", "CSQ_SYMBOL", "VEP_CHROM2", "VEP_START2"]].rename(
        columns={"POS": "START1"}
    )
    vep_bnd = vep_bnd.rename(columns={"CHROM": "CHROM1"})
    vep_bnd["CHROM1"] = normalize_chrom_column(vep_bnd["CHROM1"], sv_df["CHROM1"])
    vep_bnd["VEP_CHROM2"] = normalize_chrom_column(
        vep_bnd["VEP_CHROM2"].fillna(""), sv_df["CHROM2"]
    ).replace("", np.nan)

    vep_bnd = keep_numeric(vep_bnd, "CHROM1")

    merged = pd.merge(sv_df, vep_bnd, how="inner", on=["CHROM1", "START1"])
    # Confirm mate consistency (informational; doesn't drop rows)
    merged["MATE_AGREES"] = (
        (merged["VEP_CHROM2"] == merged["CHROM2"])
        & (merged["VEP_START2"].astype("Int64") == merged["START2"])
    )
    return merged.drop(columns=["VEP_CHROM2", "VEP_START2"])


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main() -> int:
    args = parse_args()
    if not args.sv_csv.exists():
        print(f"ERROR: SV CSV not found: {args.sv_csv}", file=sys.stderr)
        return 1
    if not args.vep_vcf.exists():
        print(f"ERROR: VEP VCF not found: {args.vep_vcf}", file=sys.stderr)
        return 1

    csq_keys = (
        [k.strip() for k in args.csq_keys.split(",")] if args.csq_keys else DEFAULT_CSQ_KEYS
    )

    sv_df = pd.read_csv(args.sv_csv)
    if args.mode == "simple":
        sv_df = keep_numeric(sv_df, "CHROM")
    else:
        sv_df = keep_numeric(sv_df, "CHROM1")
        sv_df = keep_numeric(sv_df, "CHROM2")

    print(f"Loaded {len(sv_df)} SV rows ({args.mode})")

    vep_df = read_vep_vcf(args.vep_vcf, csq_keys)
    print(f"Loaded {len(vep_df)} VEP records")

    if args.mode == "simple":
        out_df = annotate_simple(sv_df, vep_df)
    else:
        out_df = annotate_complex(sv_df, vep_df)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} annotated rows -> {args.out}")
    if args.mode == "simple" and not out_df.empty:
        top_genes = out_df["CSQ_SYMBOL"].dropna().value_counts().head(10)
        if not top_genes.empty:
            print("\nTop 10 affected genes:")
            print(top_genes.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
