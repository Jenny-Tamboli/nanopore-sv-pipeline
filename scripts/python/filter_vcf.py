#!/usr/bin/env python3
"""
Filter a Sniffles2 VCF to keep only:
  - rows with FILTER == PASS
  - rows with PRECISE in the INFO field (i.e. precise breakpoints)

The header is preserved so the original sample ID and file identifier can be
back-traced. Non-standard chromosomes are dropped.

Usage:
    python filter_vcf.py input.vcf -o output.filtered.vcf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable

# Standard human chromosomes (GRCh37/38). Override with --chroms if needed.
STANDARD_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
# Also accept the no-"chr" convention
STANDARD_CHROMS |= {str(i) for i in range(1, 23)} | {"X", "Y", "MT", "M"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", type=Path, help="Input VCF (uncompressed)")
    p.add_argument("-o", "--output", type=Path, required=True, help="Output filtered VCF")
    p.add_argument(
        "--chroms",
        type=str,
        default=None,
        help="Comma-separated list of allowed chromosomes (default: chr1-22, chrX, chrY, chrM).",
    )
    p.add_argument(
        "--keep-imprecise",
        action="store_true",
        help="Keep IMPRECISE calls (by default only PRECISE rows are kept).",
    )
    return p.parse_args()


def info_has_flag(info_field: str, flag: str) -> bool:
    """Return True if a VCF INFO field contains a flag (e.g. PRECISE)."""
    return any(part == flag for part in info_field.split(";"))


def iter_filtered(
    lines: Iterable[str],
    allowed_chroms: set[str],
    keep_imprecise: bool,
) -> Iterable[tuple[str, bool]]:
    """Yield (line, is_kept) tuples. Header lines are always kept."""
    for raw in lines:
        line = raw.rstrip("\n")
        if not line:
            continue
        if line.startswith("#"):
            yield line, True
            continue

        fields = line.split("\t")
        if len(fields) < 8:
            # Malformed; skip
            yield line, False
            continue

        chrom, _pos, _id, _ref, _alt, _qual, filt, info = fields[:8]

        if chrom not in allowed_chroms:
            yield line, False
            continue
        if filt != "PASS":
            yield line, False
            continue
        if not keep_imprecise and not info_has_flag(info, "PRECISE"):
            yield line, False
            continue

        yield line, True


def main() -> int:
    args = parse_args()

    if args.chroms:
        allowed = set(c.strip() for c in args.chroms.split(",") if c.strip())
    else:
        allowed = STANDARD_CHROMS

    if not args.vcf.exists():
        print(f"ERROR: input VCF not found: {args.vcf}", file=sys.stderr)
        return 1

    kept = 0
    dropped = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.vcf.open() as fh_in, args.output.open("w") as fh_out:
        for line, is_kept in iter_filtered(fh_in, allowed, args.keep_imprecise):
            if line.startswith("#"):
                fh_out.write(line + "\n")
                continue
            if is_kept:
                fh_out.write(line + "\n")
                kept += 1
            else:
                dropped += 1

    print(f"Filtered VCF written to {args.output}")
    print(f"  Kept:    {kept}")
    print(f"  Dropped: {dropped}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
