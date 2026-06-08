# Circos / Circa setup

The CSVs produced by `vcf_to_csv.py` are designed to drop straight into
[Circa](http://omgenomics.com/circa/) (a Circos-plot web UI). The same files
also work as input for the classic `circos` Perl tool with minimal reshaping.

## Track design

| Track | Source CSV | Meaning |
|-------|------------|---------|
| Rectangle | `<sample>_simple_svs.csv` | Span of each SV within a chromosome |
| Scatter   | `<sample>_simple_svs.csv` | Ratio of SVTYPEs per chromosome |
| Connection | `<sample>_complex_svs.csv` | BND / translocation links |

## SVTYPE colour code

To stay consistent with the histograms:

| SVTYPE | Colour | Hex |
|--------|--------|-----|
| DEL | Maraschino red | `#D7263D` |
| INS | Banana yellow  | `#F4D35E` |
| DUP | Plum / burgundy | `#7B2D26` |
| INV | Fern green     | `#2E8B57` |

## Which CSV to use

Use the **`no_filter`** versions (no size cutoff applied). Translocations
(BNDs) have no length and would be excluded by a size threshold, leaving
the connection track empty.

If you want a "high-confidence only" view, filter the simple-SV CSV by
`size > 1000` *before* loading into Circa, but keep the complex-SV CSV
unfiltered.

## Circa loading order

1. Open Circa, set genome to GRCh37 or GRCh38 to match your reference.
2. Add a **rectangle** track → load `<sample>_simple_svs.csv`, map
   `chrom`/`start`/`stop` to the corresponding fields, colour by `svtype`.
3. Add a **scatter** track → same CSV, compute SVTYPE proportions per
   chromosome bin (1 Mb is a reasonable default).
4. Add a **connection** track → load `<sample>_complex_svs.csv`, map
   `chrom1`/`start1` and `chrom2`/`start2` as the two endpoints.

## Sanity check before plotting

The report's cross-verification step is worth doing in code too:

```python
import pandas as pd
df = pd.read_csv("results/<sample>/<sample>_simple_svs.csv")
print(df["svtype"].value_counts())          # expected types only
print(df.groupby("chrom")["svtype"].count()) # spread across chroms
```

If a chromosome shows a wildly different SV count from neighbours of similar
size, that's usually a sign of a missing region in the BAM or a contig name
mismatch.
