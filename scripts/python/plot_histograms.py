#!/usr/bin/env python3
"""
Generate histograms of structural-variant sizes from a simple-SVs CSV.

For each (sample, cutoff) pair this writes one PNG showing size distribution
faceted by SVTYPE. Cutoffs follow the report: no_filter, >1kb, >10kb.

BNDs are intentionally not handled here — they have no intrinsic length.
Use the complex-SV CSV with Circos instead.

Usage:
    python plot_histograms.py results/mySample/mySample_simple_svs.csv \
        --sample mySample --outdir results/mySample/histplots/

    # restrict to specific chromosomes
    python plot_histograms.py ... --chroms chr9,chr17,chr18
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.use("Agg")  # non-interactive backend, safe for HPC

# Same colour palette referenced in the report's Circos description, so plots
# stay visually consistent across histograms and Circos.
SVTYPE_COLOURS = {
    "DEL": "#D7263D",  # Maraschino red
    "INS": "#F4D35E",  # Banana yellow
    "DUP": "#7B2D26",  # Plum / burgundy
    "INV": "#2E8B57",  # Fern green
}

CUTOFFS = {
    "no_filter": 0,
    "gt1kb": 1_000,
    "gt10kb": 10_000,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", type=Path, help="Common-SVs CSV from vcf_to_csv.py")
    p.add_argument("--sample", required=True, help="Sample name (used in titles/filenames)")
    p.add_argument("--outdir", type=Path, default=Path("histplots"), help="Output directory")
    p.add_argument(
        "--chroms",
        type=str,
        default=None,
        help="Optional comma-separated chromosomes to restrict to, e.g. chr9,chr17,chr18.",
    )
    p.add_argument("--bins", type=int, default=50, help="Number of histogram bins (default 50)")
    p.add_argument("--log-x", action="store_true", help="Use log10 scale on size axis")
    return p.parse_args()


def plot_one(
    df: pd.DataFrame,
    sample: str,
    cutoff_name: str,
    cutoff_value: int,
    outdir: Path,
    bins: int,
    log_x: bool,
) -> Path:
    """Plot a single (sample, cutoff) histogram."""
    sub = df[df["size"] > cutoff_value].copy()

    svtypes_present = [s for s in SVTYPE_COLOURS if (sub["svtype"] == s).any()]
    n = max(len(svtypes_present), 1)

    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, svt in zip(axes, svtypes_present):
        sizes = sub.loc[sub["svtype"] == svt, "size"].to_numpy()
        if log_x and sizes.size > 0:
            sizes = np.log10(sizes + 1)
            ax.set_xlabel("log10(size + 1)")
        else:
            ax.set_xlabel("size (bp)")

        ax.hist(sizes, bins=bins, color=SVTYPE_COLOURS[svt], edgecolor="black", linewidth=0.4)
        ax.set_title(f"{svt}  (n={len(sizes)})")
        ax.set_ylabel("count")
        ax.spines[["top", "right"]].set_visible(False)

    if not svtypes_present:
        axes[0].text(
            0.5, 0.5, "no SVs above cutoff",
            ha="center", va="center", transform=axes[0].transAxes,
        )
        axes[0].set_axis_off()

    cutoff_label = {"no_filter": "no filter", "gt1kb": "size > 1 kb", "gt10kb": "size > 10 kb"}[cutoff_name]
    fig.suptitle(f"{sample} — SV sizes ({cutoff_label})", fontsize=13)
    fig.tight_layout()

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"{sample}_{cutoff_name}.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    args = parse_args()
    if not args.csv.exists():
        print(f"ERROR: CSV not found: {args.csv}", file=sys.stderr)
        return 1

    df = pd.read_csv(args.csv)
    if df.empty:
        print("WARNING: CSV is empty; nothing to plot.")
        return 0

    required = {"chrom", "size", "svtype"}
    missing = required - set(df.columns)
    if missing:
        print(f"ERROR: CSV is missing columns: {missing}", file=sys.stderr)
        return 2

    if args.chroms:
        wanted = {c.strip() for c in args.chroms.split(",") if c.strip()}
        df = df[df["chrom"].isin(wanted)].copy()
        suffix = "_" + "-".join(sorted(wanted))
    else:
        suffix = ""

    print(f"Loaded {len(df)} simple SVs from {args.csv}")
    for cutoff_name, cutoff_value in CUTOFFS.items():
        out = plot_one(
            df,
            sample=args.sample + suffix,
            cutoff_name=cutoff_name,
            cutoff_value=cutoff_value,
            outdir=args.outdir,
            bins=args.bins,
            log_x=args.log_x,
        )
        print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
