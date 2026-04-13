#!/usr/bin/env python3
"""
Convert a Bracken *_bracken.report to MetaPhlAn3 format for use with HUMAnN3
--taxonomic-profile.

Usage:
    python bracken_to_mpa.py <sample_bracken.report> <sample.mpa>

Input: Bracken-style report (same column layout as Kraken2 report, with
       species-level read counts redistributed by Bracken).

Output: MetaPhlAn3 tab-separated format:
    #mpa_v3_CHOCOPhlAn_201901
    k__Bacteria\t96.61000
    k__Bacteria|p__Bacteroidota\t16.93000
    ...
    k__...|s__Bacteroides_nordii\t3.19000

Only the 7 standard NCBI ranks (D/P/C/O/F/G/S) are emitted; intermediate
ranks (K, K1, K2, R, R1, etc.) are skipped. Percentages are taken directly
from column 1 of the Bracken report (% of all reads in the sample, including
unclassified) — no re-normalisation is applied.

Data loss notes:
- Unclassified reads are not included (Kraken2 never assigns them to a taxon).
- Reads assigned above species level in Kraken2 are redistributed to species
  by Bracken before this script runs, so using _bracken.report (not the raw
  Kraken2 report) minimises species-level gaps.
- Reads at sub-species ranks (S1, S2, …) are already rolled into their parent
  S-level count in the Bracken report, so they are not lost here.
- Species whose names do not appear in HUMAnN3's ChocoPhlAn database will
  fall back to the translated (UniRef90 diamond) search only; they still
  contribute to pathway/gene-family tables but without nucleotide-level hits.
  
Credit: This file was partially generated using using Claude 3.5 Sonnet by Anthropic - AH
"""

import sys

# Kraken2 rank code -> MetaPhlAn prefix character
RANK_PREFIXES = {
    'D': 'k',  # Domain / Superkingdom -> k__ (MetaPhlAn convention)
    'P': 'p',
    'C': 'c',
    'O': 'o',
    'F': 'f',
    'G': 'g',
    'S': 's',
}

RANK_ORDER = ['D', 'P', 'C', 'O', 'F', 'G', 'S']


def convert(report_path, output_path):
    lineage = {}   # rank -> "prefix__Name" for the current tree path
    entries = []   # list of (clade_string, percentage)

    with open(report_path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue

            pct  = float(parts[0].strip())
            rank = parts[3].strip()
            name = parts[5].strip()   # strip leading indentation spaces

            if rank not in RANK_PREFIXES:
                continue

            prefix     = RANK_PREFIXES[rank]
            clean_name = name.replace(' ', '_')
            tag        = f'{prefix}__{clean_name}'

            # Set current rank in lineage and clear all deeper ranks.
            # Because Kraken reports are in DFS (parent-before-child) order,
            # this correctly builds the path as we traverse the tree.
            rank_idx = RANK_ORDER.index(rank)
            lineage[rank] = tag
            for deeper in RANK_ORDER[rank_idx + 1:]:
                lineage.pop(deeper, None)

            clade = '|'.join(lineage[r] for r in RANK_ORDER if r in lineage)
            entries.append((clade, pct))

    with open(output_path, 'w') as fh:
        # Header must contain "v3" for HUMAnN3's version check.
        # Format matches MetaPhlAn3 3-column output: clade_name\tNCBI_tax_id\trelative_abundance
        fh.write('#v30_CHOCOPhlAn_201901\n')
        for clade, pct in entries:
            fh.write(f'{clade}\t1\t{pct:.5f}\n')

    return len(entries)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'Usage: {sys.argv[0]} <bracken_report> <output_mpa>')
        sys.exit(1)
    n = convert(sys.argv[1], sys.argv[2])
    print(f'Wrote {n} entries to {sys.argv[2]}')
