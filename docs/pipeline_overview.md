# Pipeline overview

Detailed notes on each stage. The high-level diagram lives in the top-level README.

## 1. Alignment — NGMLR

NGMLR is a long-read aligner explicitly designed for ONT/PacBio reads and is the
recommended upstream for Sniffles. Key flag:

```bash
ngmlr -x ont -r reference.fa -q reads.fastq -o sample.sam
```

`-x ont` tunes scoring parameters for Nanopore error profiles. For PacBio HiFi
use `-x pacbio` instead.

Output is SAM; we convert to BAM immediately to save disk.

## 2. Sort & index — samtools

Sniffles2 requires a coordinate-sorted, indexed BAM:

```bash
samtools sort -@ 8 -o sample.sorted.bam sample.unsorted.bam
samtools index sample.sorted.bam
```

`flagstat` output is written alongside the BAM as a quick QC sanity check.

## 3. SV calling — Sniffles2

```bash
sniffles --input sample.sorted.bam \
         --reference reference.fa \
         --vcf sample.sniffles.vcf \
         --sample-id sample \
         --minsupport auto \
         --minsvlen 50
```

`--minsvlen 50` is the Sniffles2 default; raise it (e.g. `100`) if you only
care about larger events. `--minsupport auto` lets the caller adapt to coverage.

## 4. Filtering — `filter_vcf.py`

Two criteria, both from the report's spec:

| Criterion | Implementation |
|-----------|----------------|
| `FILTER == PASS` | Column 7 of the VCF |
| `PRECISE` flag in INFO | Token in column 8 (case-sensitive) |

Header lines are passed through unchanged so the original sample-ID and
file-identifier provenance is preserved.

Non-standard contigs (`chr18_gl000207_random` etc.) are dropped by default.
Override with `--chroms`.

## 5. VCF → CSV — `vcf_to_csv.py`

Splits records into two tables based on SVTYPE:

- **Simple SVs** (DEL, INS, DUP, INV) → one row per event with a meaningful
  `size` field.
- **Complex SVs** (BND, optionally TRA) → one row per breakend pair. No `size`
  is recorded because BNDs only describe two coordinates, possibly on different
  chromosomes.

INFO parsing handles both `=value` pairs and bare flags. The BND `ALT` regex
matches all four orientations of the VCF 4.x BND spec (`N[chr:pos[`,
`N]chr:pos]`, `[chr:pos[N`, `]chr:pos]N`).

## 6. Histograms — `plot_histograms.py`

For each cutoff (`no_filter`, `>1kb`, `>10kb`), one figure per sample. SVTYPE
panels share the report's Circos colour palette so histograms and Circos plots
look consistent side-by-side.

`--chroms chr9,chr17,chr18` restricts to a subset; the chosen chromosomes are
appended to the filename so the same outdir can hold both global and
per-chromosome figures.

## 7. Circos / Circa

See `circos_setup.md`. The CSVs written in step 5 are already in the layout
expected by Circa's rectangle, scatter, and connection tracks.

## Why BNDs aren't size-filtered

A BND record describes a single breakend (or a pair via the mate ID). The
event has no inherent length — it is two genomic coordinates joined by
evidence. `>1kb` / `>10kb` cutoffs are therefore only meaningful for the
simple-SV table.

For Circos, use the `no_filter` CSVs so translocations are still represented.
