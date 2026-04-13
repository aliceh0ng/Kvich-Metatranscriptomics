#!/usr/bin/env python3
"""
Parse gene_id, gene name, and product from a GTF file.
- gene_id and gene name come from 'gene' feature lines
- product comes from 'CDS' feature lines (joined by gene_id)
"""

import re
import sys


def extract_attr(attrs: str, key: str) -> str:
    """Return the value of a GTF attribute, or empty string if absent."""
    m = re.search(rf'{key} "([^"]+)"', attrs)
    return m.group(1) if m else ""


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <gtf_file>", file=sys.stderr)
        sys.exit(1)

    gtf_path = sys.argv[1]

    # Pass 1: collect gene_id -> product from CDS lines
    # A gene may have multiple CDS lines; take the first non-empty product seen.
    cds_product: dict[str, str] = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = fields[8]
            gene_id = extract_attr(attrs, "gene_id")
            product = extract_attr(attrs, "product")
            if gene_id and product and gene_id not in cds_product:
                cds_product[gene_id] = product

    # Pass 2: emit one row per gene line, joining product from CDS
    print("gene_id\tgene\tproduct")
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            gene_id = extract_attr(attrs, "gene_id")
            gene    = extract_attr(attrs, "gene")
            product = cds_product.get(gene_id, "")
            print(f"{gene_id}\t{gene}\t{product}")


if __name__ == "__main__":
    main()

# Credit: This file was partially generated using using Claude 3.5 Sonnet by Anthropic - AH
