#!/usr/bin/env python3
"""
Convert a filtered Sniffles2 VCF into two CSV files:

  <sample>_simple_svs.csv   -> DEL, INS, DUP, INV with SIZE > 1 kb
      columns: CHROM, START, STOP, SIZE, SVTYPE, STRAND, SIZE_kb
  <sample>_complex_svs.csv  -> BND (translocations / complex rearrangements)
      columns: CHROM1, START1, CHROM2, START2, SVTYPE

Why a 1 kb filter only on simple SVs:
  Sniffles reports SVLEN / END for DEL/INS/DUP/INV, so a size threshold is
  meaningful. BNDs describe two breakend coordinates and have no intrinsic
  length, so the same filter cannot be applied.

Why numeric chromosomes only:
  Sex chromosomes (chrX, chrY), mitochondrion (chrM), and decoy contigs
  (chr*_random, chrUn_*) behave differently from autosomes. The filter
  step should have already dropped them; this is a second line of defence.

These CSVs feed directly into:
  - plot_histograms.py (chromosome-size-normalized SV burden)
  - annotate_vep.py    (merge with VEP web-server output)
  - Circa / Circos     (see docs/circos_setup.md)

Usage:
    python vcf_to_csv.py input.filtered.vcf --sample mySample --outdir results/mySample/

    # Custom size threshold for simple SVs (default 1000 bp):
    python vcf_to_csv.py input.filtered.vcf --sample mySample --min-size 500
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

# Match Sniffles2 contigs of the form "chr1".."chr22" or "1".."22".
# Sex chromosomes, mitochondrion, alt/decoy contigs all excluded.
NUMERIC_CHROMS = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}

SIMPLE_SVTYPES = {"DEL", "INS", "DUP", "INV"}
COMPLEX_SVTYPES = {"BND", "TRA"}  # TRA included for callers that emit it

# Matches the BND ALT spec: N[chr5:12345[, ]chr5:12345]N, etc.
BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", type=Path, help="Input filtered VCF")
    p.add_argument("--sample", required=True, help="Sample name (used in output filenames)")
    p.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
    p.add_argument(
        "--min-size",
        type=int,
        default=1000,
        help="Minimum SVLEN (bp) to keep for simple SVs. Default 1000 (= 1 kb). Not applied to BNDs.",
    )
    return p.parse_args()


def parse_info(info: str) -> dict[str, str]:
    """Parse a VCF INFO column into a dict. Flags (no '=') map to themselves."""
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


def parse_bnd_alt(alt: str) -> Optional[tuple[str, int]]:
    """Extract (mate_chrom, mate_pos) from a BND ALT string."""
    m = BND_RE.search(alt)
    if not m:
        return None
    return m.group(1), int(m.group(2))


def read_vcf_records(vcf_path: Path) -> pd.DataFrame:
    """Load VCF data lines into a DataFrame. Header lines are skipped here
    (they're parsed separately by other scripts when needed)."""
    rows = []
    with vcf_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            rows.append({
                "CHROM": fields[0],
                "POS": int(fields[1]),
                "ID": fields[2],
                "REF": fields[3],
                "ALT": fields[4],
                "QUAL": fields[5],
                "FILTER": fields[6],
                "INFO": fields[7],
            })
    return pd.DataFrame(rows)


def build_simple(df: pd.DataFrame, min_size_bp: int) -> pd.DataFrame:
    """Build the simple-SV dataframe (DEL/INS/DUP/INV), filtered to size >= min_size_bp.

    Columns:
        CHROM, START, STOP, SIZE, SVTYPE, STRAND, SIZE_kb

    SIZE is always positive (length of the SV in bp). The VCF spec convention
    reports DEL SVLEN as negative; we store the absolute value here since the
    SVTYPE column already encodes the direction (DEL vs INS/DUP/INV). This
    keeps downstream Excel/pandas inspection unambiguous.
    """
    records = []
    for _, row in df.iterrows():
        info = parse_info(row["INFO"])
        svtype = info.get("SVTYPE", "")
        if svtype not in SIMPLE_SVTYPES:
            continue

        start = int(row["POS"])
        end = int(info["END"]) if "END" in info else start + abs(int(info.get("SVLEN", 0)))
        # SIZE is the absolute length. The VCF spec stores DEL SVLEN as
        # negative; we take abs() because SVTYPE already encodes direction.
        svlen = abs(int(info.get("SVLEN", end - start)))
        if svlen < min_size_bp:
            continue

        strand = info.get("STRANDS", info.get("STRAND", "."))

        records.append({
            "CHROM": row["CHROM"],
            "START": start,
            "STOP": end,
            "SIZE": svlen,
            "SVTYPE": svtype,
            "STRAND": strand,
            "SIZE_kb": round(svlen / 1000, 3),
        })
    return pd.DataFrame(records, columns=["CHROM", "START", "STOP", "SIZE", "SVTYPE", "STRAND", "SIZE_kb"])


def build_complex(df: pd.DataFrame) -> pd.DataFrame:
    """Build the complex-SV (BND) dataframe.

    BNDs have no intrinsic length, so no size filter is applied. Each row
    represents one breakend pair: (CHROM1, START1) ⟷ (CHROM2, START2).

    A BND record can be on a non-numeric chromosome at one end but link to a
    numeric chromosome at the other end. We require both ends to be on the
    numeric autosome whitelist to keep the table self-consistent with the
    simple-SV table.
    """
    records = []
    for _, row in df.iterrows():
        info = parse_info(row["INFO"])
        svtype = info.get("SVTYPE", "")
        if svtype not in COMPLEX_SVTYPES:
            continue

        mate = parse_bnd_alt(row["ALT"])
        if mate is None:
            # Some callers put the mate in INFO instead
            if "CHR2" in info and "END" in info:
                mate = (info["CHR2"], int(info["END"]))
            else:
                continue
        mate_chrom, mate_pos = mate

        if mate_chrom not in NUMERIC_CHROMS:
            continue

        records.append({
            "CHROM1": row["CHROM"],
            "START1": int(row["POS"]),
            "CHROM2": mate_chrom,
            "START2": mate_pos,
            "SVTYPE": svtype,
        })
    return pd.DataFrame(records, columns=["CHROM1", "START1", "CHROM2", "START2", "SVTYPE"])


def main() -> int:
    args = parse_args()
    if not args.vcf.exists():
        print(f"ERROR: VCF not found: {args.vcf}", file=sys.stderr)
        return 1

    args.outdir.mkdir(parents=True, exist_ok=True)

    df = read_vcf_records(args.vcf)
    if df.empty:
        print("WARNING: VCF contained no data lines.")
        return 0

    # Inspect SVTYPEs present (transparency for the user)
    svtype_counts = df["INFO"].apply(lambda s: parse_info(s).get("SVTYPE", "?")).value_counts()
    print("Unique SVTYPEs in file:")
    print(svtype_counts.to_string())
    print()

    # Drop non-numeric chromosomes (belt-and-suspenders; filter_vcf.py should
    # have done this already, but VCFs from other pipelines might not have).
    n_before = len(df)
    df = df[df["CHROM"].isin(NUMERIC_CHROMS)].reset_index(drop=True)
    n_dropped = n_before - len(df)
    if n_dropped:
        print(f"Dropped {n_dropped} records on non-numeric chromosomes (chrX/Y/M/decoys).")

    simple_df = build_simple(df, min_size_bp=args.min_size)
    complex_df = build_complex(df)

    simple_path = args.outdir / f"{args.sample}_simple_svs.csv"
    complex_path = args.outdir / f"{args.sample}_complex_svs.csv"

    simple_df.to_csv(simple_path, index=False)
    complex_df.to_csv(complex_path, index=False)

    print(f"Wrote {len(simple_df):>6d} simple SVs (|SIZE| >= {args.min_size} bp)  -> {simple_path}")
    print(f"Wrote {len(complex_df):>6d} complex SVs                              -> {complex_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
