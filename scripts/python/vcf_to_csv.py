#!/usr/bin/env python3
"""
Convert a filtered Sniffles2 VCF into two CSV files:

  <sample>_simple_svs.csv   -> DEL, INS, DUP, INV
      columns: chrom, start, stop, size, svtype, strand
  <sample>_complex_svs.csv  -> BND (translocations / complex rearrangements)
      columns: chrom1, start1, end1, chrom2, start2, end2, svtype

These CSVs are designed to feed directly into Circa / Circos
(see docs/circos_setup.md) and the histogram script.

Usage:
    python vcf_to_csv.py input.filtered.vcf --sample mySample --outdir results/mySample/

Pseudo-code (from the report):
  - Inspect the unique SVTYPEs present in the file
  - Drop non-standard chromosomes
  - Iterate over ALT and INFO columns, parsing SVTYPE / SVLEN / END / STRAND
  - Build the two dataframes and write them out
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
STANDARD_CHROMS |= {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}

SIMPLE_SVTYPES = {"DEL", "INS", "DUP", "INV"}
COMPLEX_SVTYPES = {"BND", "TRA"}  # TRA included for tools that emit it

# Matches the BND ALT spec, e.g. N[chr5:12345[, ]chr5:12345]N, etc.
BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", type=Path, help="Input filtered VCF")
    p.add_argument("--sample", required=True, help="Sample name (used in output filenames)")
    p.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
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
    """Load VCF data lines into a DataFrame."""
    rows = []
    with vcf_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            rows.append(
                {
                    "CHROM": fields[0],
                    "POS": int(fields[1]),
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "QUAL": fields[5],
                    "FILTER": fields[6],
                    "INFO": fields[7],
                }
            )
    return pd.DataFrame(rows)


def build_simple(df: pd.DataFrame) -> pd.DataFrame:
    """Build the simple-SV dataframe (DEL/INS/DUP/INV)."""
    records = []
    for _, row in df.iterrows():
        info = parse_info(row["INFO"])
        svtype = info.get("SVTYPE", "")
        if svtype not in SIMPLE_SVTYPES:
            continue

        start = int(row["POS"])
        # END for symbolic SVs; fall back to POS + |SVLEN|
        end = int(info["END"]) if "END" in info else start + abs(int(info.get("SVLEN", 0)))
        size = abs(int(info.get("SVLEN", end - start)))
        strand = info.get("STRANDS", info.get("STRAND", "."))

        records.append(
            {
                "chrom": row["CHROM"],
                "start": start,
                "stop": end,
                "size": size,
                "svtype": svtype,
                "strand": strand,
            }
        )
    return pd.DataFrame(records, columns=["chrom", "start", "stop", "size", "svtype", "strand"])


def build_complex(df: pd.DataFrame) -> pd.DataFrame:
    """Build the complex-SV (BND) dataframe."""
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

        # BNDs have no intrinsic length; record both breakend coordinates
        records.append(
            {
                "chrom1": row["CHROM"],
                "start1": int(row["POS"]),
                "end1": int(row["POS"]),
                "chrom2": mate_chrom,
                "start2": mate_pos,
                "end2": mate_pos,
                "svtype": svtype,
            }
        )
    return pd.DataFrame(
        records,
        columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "svtype"],
    )


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

    # Inspect SVTYPEs present (per the report's pseudo-code)
    svtype_counts = df["INFO"].apply(lambda s: parse_info(s).get("SVTYPE", "?")).value_counts()
    print("Unique SVTYPEs in file:")
    print(svtype_counts.to_string())
    print()

    # Drop non-standard chromosomes
    n_before = len(df)
    df = df[df["CHROM"].isin(STANDARD_CHROMS)].reset_index(drop=True)
    print(f"Dropped {n_before - len(df)} records on non-standard chromosomes.")

    simple_df = build_simple(df)
    complex_df = build_complex(df)

    simple_path = args.outdir / f"{args.sample}_simple_svs.csv"
    complex_path = args.outdir / f"{args.sample}_complex_svs.csv"

    simple_df.to_csv(simple_path, index=False)
    complex_df.to_csv(complex_path, index=False)

    print(f"Wrote {len(simple_df):>6d} simple SVs   -> {simple_path}")
    print(f"Wrote {len(complex_df):>6d} complex SVs  -> {complex_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
