#!/usr/bin/env python3
"""
Plot per-chromosome SV burden normalized by chromosome size.

For each autosome (chr1..chr22), this script computes the fraction of the
chromosome affected by each simple SV type (DEL, INS, DUP, INV):

    affected_fraction[svtype, chrom] = sum(|SIZE| for that svtype on that chrom)
                                       / chrom_size

It then produces a stacked horizontal bar plot, one bar per chromosome, with
segments coloured by SVTYPE and a trailing "Other" segment representing the
unaffected fraction (1 - sum of SV fractions). This visualizes the relative
SV burden across chromosomes regardless of their size.

Chromosome sizes are parsed automatically from the input VCF's ##contig=<...>
header lines (one source of truth — no hardcoded reference assumptions).

Why this is more useful than raw SV counts:
  - Larger chromosomes naturally have more SVs. Normalizing by chromosome
    size lets you compare burden meaningfully across chromosomes.
  - Stacking by SVTYPE shows the relative contribution of each variant
    class to the total affected fraction.

Input:
  - Simple-SVs CSV from vcf_to_csv.py (columns: CHROM, START, STOP, SIZE,
    SVTYPE, STRAND, SIZE_kb). The 1 kb filter should already be applied.
  - The Sniffles VCF used to produce it (for chromosome-size lookup).

Output:
  - One PNG: <sample>_sv_burden.png

Usage:
    python plot_histograms.py \\
        results/mySample/mySample_simple_svs.csv \\
        --vcf results/mySample/mySample.sniffles.vcf \\
        --sample mySample \\
        --outdir results/mySample/

    # Restrict to specific chromosomes (still normalized correctly):
    python plot_histograms.py ... --chroms chr9,chr17,chr18
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")  # non-interactive backend; safe for HPC and CI

# Colour palette consistent with the Circos plots (docs/circos_setup.md).
SVTYPE_COLOURS = {
    "DEL": "#D7263D",  # red
    "INS": "#F4D35E",  # yellow
    "DUP": "#7B2D26",  # burgundy
    "INV": "#2E8B57",  # green
}
OTHER_COLOUR = "#E0E0E0"  # neutral grey for the unaffected fraction

# Order SV types are stacked left-to-right
SVTYPE_ORDER = ["DEL", "INS", "DUP", "INV"]

# Matches ##contig=<ID=chr1,length=248956422,...>
CONTIG_RE = re.compile(r"##contig=<.*?ID=([^,>]+).*?length=(\d+)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", type=Path, help="Simple-SVs CSV from vcf_to_csv.py (1 kb filtered)")
    p.add_argument("--vcf", type=Path, required=True,
                   help="Source Sniffles VCF (used to read ##contig header lines for chromosome sizes)")
    p.add_argument("--sample", required=True, help="Sample name (used in titles and filenames)")
    p.add_argument("--outdir", type=Path, default=Path("."), help="Output directory")
    p.add_argument("--chroms", type=str, default=None,
                   help="Comma-separated chromosomes to plot (default: all autosomes present)")
    return p.parse_args()


def parse_contig_sizes(vcf_path: Path) -> dict[str, int]:
    """Parse `##contig=<ID=...,length=...>` lines from a VCF header.

    Returns a dict mapping contig name to length (bp). Reading stops at the
    first non-header line for efficiency on large VCFs.
    """
    sizes: dict[str, int] = {}
    with vcf_path.open() as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            m = CONTIG_RE.search(line)
            if m:
                sizes[m.group(1)] = int(m.group(2))
    return sizes


def build_burden_table(
    sv_df: pd.DataFrame,
    chrom_sizes: dict[str, int],
    wanted_chroms: list[str] | None,
) -> pd.DataFrame:
    """Compute per-chromosome SV burden, normalized by chromosome size.

    Returns a DataFrame indexed by chromosome with columns:
        DEL, INS, DUP, INV, Sum_SVs, Other, Chromosome_size_Mbp

    Each SVTYPE column is the affected fraction (0-1). 'Other' is
    1 - Sum_SVs, capped at 0 in case overlapping SVs push the total above 1.
    """
    # SIZE is positive by construction (see vcf_to_csv.py). Aggregate sum
    # of bases affected per (chrom, svtype).
    pivot = (
        sv_df.groupby(["CHROM", "SVTYPE"])["SIZE"]
        .sum()
        .unstack(fill_value=0)
        .reindex(columns=SVTYPE_ORDER, fill_value=0)
    )

    # Determine which chromosomes to keep. Default: any autosome present in
    # the VCF header. The CSV itself might not cover every chromosome (sample
    # may have no SVs on some), so we drive the order from the header.
    autosomes_in_header = [c for c in chrom_sizes
                           if c in {f"chr{i}" for i in range(1, 23)}
                           or c in {str(i) for i in range(1, 23)}]
    # Sort autosomes numerically (chr1, chr2, ..., chr22)
    def _chrom_key(c: str) -> int:
        return int(c.replace("chr", ""))
    autosomes_in_header.sort(key=_chrom_key)

    if wanted_chroms:
        chroms = [c for c in autosomes_in_header if c in wanted_chroms]
    else:
        chroms = autosomes_in_header

    if not chroms:
        return pd.DataFrame()

    # Ensure every wanted chromosome appears, even with zero SVs
    pivot = pivot.reindex(chroms, fill_value=0)

    # Normalize each row by that chromosome's size (bp)
    sizes_bp = pd.Series({c: chrom_sizes[c] for c in chroms}, name="size_bp")
    burden = pivot.div(sizes_bp, axis=0)

    burden["Sum_SVs"] = burden[SVTYPE_ORDER].sum(axis=1)
    burden["Other"] = (1.0 - burden["Sum_SVs"]).clip(lower=0.0)
    burden["Chromosome_size_Mbp"] = (sizes_bp / 1_000_000).round(2)

    return burden


def plot_stacked(burden: pd.DataFrame, sample: str, out_path: Path) -> None:
    """Render the stacked horizontal bar plot."""
    chroms = burden.index.tolist()
    n = len(chroms)
    y = np.arange(n)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.35 * n + 1)))

    left = np.zeros(n)
    for svt in SVTYPE_ORDER:
        widths = burden[svt].values
        ax.barh(y, widths, left=left, color=SVTYPE_COLOURS[svt],
                edgecolor="white", linewidth=0.4, label=svt)
        left += widths

    # Add the unaffected remainder
    other = burden["Other"].values
    ax.barh(y, other, left=left, color=OTHER_COLOUR,
            edgecolor="white", linewidth=0.4, label="Other (unaffected)")

    ax.set_yticks(y)
    ax.set_yticklabels(chroms)
    ax.invert_yaxis()  # chr1 at the top
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Fraction of chromosome affected by SV type")
    ax.set_title(f"{sample} — SV burden per chromosome (size-normalized, |SIZE| > 1 kb)",
                 fontsize=12)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="lower right", frameon=False, fontsize=9, ncol=5,
              bbox_to_anchor=(1.0, -0.18 / max(1, n / 10)))

    # Annotate each bar with its Sum_SVs as a percentage, for quick reading
    for yi, total in zip(y, burden["Sum_SVs"].values):
        if total > 0:
            ax.text(min(total + 0.005, 0.99), yi, f"{total*100:.2f}%",
                    va="center", ha="left", fontsize=8, color="#333333")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    if not args.csv.exists():
        print(f"ERROR: CSV not found: {args.csv}", file=sys.stderr)
        return 1
    if not args.vcf.exists():
        print(f"ERROR: VCF not found (needed for chromosome sizes): {args.vcf}", file=sys.stderr)
        return 1

    sv_df = pd.read_csv(args.csv)
    required = {"CHROM", "SIZE", "SVTYPE"}
    missing = required - set(sv_df.columns)
    if missing:
        print(f"ERROR: CSV is missing required columns: {missing}", file=sys.stderr)
        return 2

    chrom_sizes = parse_contig_sizes(args.vcf)
    if not chrom_sizes:
        print(f"ERROR: no ##contig lines found in VCF header: {args.vcf}", file=sys.stderr)
        return 3
    print(f"Loaded {len(chrom_sizes)} contig sizes from VCF header.")

    wanted = None
    if args.chroms:
        wanted = [c.strip() for c in args.chroms.split(",") if c.strip()]

    burden = build_burden_table(sv_df, chrom_sizes, wanted)
    if burden.empty:
        print("WARNING: no chromosomes to plot after filtering.")
        return 0

    suffix = "_" + "-".join(burden.index.tolist()) if wanted else ""
    out = args.outdir / f"{args.sample}{suffix}_sv_burden.png"
    plot_stacked(burden, args.sample, out)
    print(f"Wrote {out}")
    print()
    print("Summary (fraction of each chromosome affected):")
    print(burden[SVTYPE_ORDER + ["Sum_SVs", "Other"]].round(4).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
