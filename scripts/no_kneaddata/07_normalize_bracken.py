"""
Normalize Bracken abundance estimates following Kvich et al.:
  - Scaling factor = min_reads_across_samples / reads_in_sample
  - Scaled counts = bracken new_est_reads * scaling_factor
  - Threshold: remove taxa with log10(scaled_count) < 0.9 across all samples

Outputs:
  counts_unscaled.tsv   -- raw bracken new_est_reads (for within-sample comparisons)
  counts_scaled.tsv     -- depth-normalized counts (for across-sample comparisons)
  counts_filtered.tsv   -- scaled counts with low-abundance taxa removed (threshold = 0.9)
  scaling_factors.tsv   -- per-sample scaling factors

Usage (run from bmeg524/ root):
  python scripts/no_kneaddata/07_normalize_bracken.py [--bracken-dir DIR] [--out-dir DIR]

Defaults:
  --bracken-dir  no_kneaddata/community_composition/bracken
  --out-dir      no_kneaddata/community_composition/normalized
  
Credit: This file was partially generated using using Claude 3.5 Sonnet by Anthropic - AH
"""

import os
import glob
import argparse
import pandas as pd
import numpy as np

THRESHOLD = 0.9   # log10 scaled count cutoff from Kvich et al. Fig S1


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--bracken-dir", default="no_kneaddata/community_composition/bracken")
    p.add_argument("--out-dir",     default="no_kneaddata/community_composition/normalized")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ── 1. Load all Bracken outputs ────────────────────────────────────────────
    bracken_files = sorted(glob.glob(os.path.join(args.bracken_dir, "*.bracken")))
    if not bracken_files:
        raise FileNotFoundError(f"No .bracken files found in {args.bracken_dir}")

    frames = {}
    for f in bracken_files:
        sample = os.path.basename(f).replace(".bracken", "")
        df = pd.read_csv(f, sep="\t")
        frames[sample] = df.set_index("name")["new_est_reads"]

    counts = pd.DataFrame(frames).fillna(0)   # taxa x samples
    samples = counts.columns.tolist()
    print(f"Loaded {len(samples)} samples, {len(counts)} taxa")

    # ── 2. Calculate read counts per sample ───────────────────────────────────
    # Uses R1 only (paired-end, so total fragments = R1 read count)
    # input/ is symlinked at bmeg524 root → bfrag4/input (same files as bfrag3/input)
    read_counts = {}
    for sample in samples:
        r1 = os.path.join("input", f"{sample}_unmapped_R1.fq")
        result = os.popen(f"wc -l < {r1}").read().strip()
        read_counts[sample] = int(result) // 4

    print("\nRead counts per sample (sorted):")
    for s, n in sorted(read_counts.items(), key=lambda x: x[1]):
        print(f"  {s}: {n:,}")

    # ── 3. Scaling factors ─────────────────────────────────────────────────────
    min_reads = min(read_counts.values())
    min_sample = [s for s, n in read_counts.items() if n == min_reads][0]
    print(f"\nMinimum reads: {min_reads:,} ({min_sample})")

    scaling_factors = {s: min_reads / read_counts[s] for s in samples}

    # ── 4. Apply scaling ───────────────────────────────────────────────────────
    scaled = counts.copy()
    for sample in samples:
        scaled[sample] = counts[sample] * scaling_factors[sample]

    # ── 5. Apply threshold (log10 scaled count >= 0.9 in at least one sample) ──
    log10_scaled = np.log10(scaled + 1)   # +1 to avoid log10(0)
    keep = (log10_scaled >= THRESHOLD).any(axis=1)
    filtered = scaled[keep]

    print(f"\nTaxa before threshold: {len(counts)}")
    print(f"Taxa after threshold:  {len(filtered)}")

    # ── 6. Save outputs ────────────────────────────────────────────────────────
    counts.to_csv(os.path.join(args.out_dir, "counts_unscaled.tsv"), sep="\t")
    scaled.to_csv(os.path.join(args.out_dir, "counts_scaled.tsv"), sep="\t")
    filtered.to_csv(os.path.join(args.out_dir, "counts_filtered.tsv"), sep="\t")

    sf_df = pd.DataFrame({
        "sample":         list(scaling_factors.keys()),
        "total_reads":    [read_counts[s] for s in scaling_factors],
        "scaling_factor": list(scaling_factors.values()),
    })
    sf_df.to_csv(os.path.join(args.out_dir, "scaling_factors.tsv"), sep="\t", index=False)

    print(f"\nOutput files:")
    print(f"  {args.out_dir}/counts_unscaled.tsv  -- raw bracken counts")
    print(f"  {args.out_dir}/counts_scaled.tsv    -- depth-normalized counts")
    print(f"  {args.out_dir}/counts_filtered.tsv  -- scaled, threshold applied (log10 >= {THRESHOLD})")
    print(f"  {args.out_dir}/scaling_factors.tsv  -- per-sample scaling factors")


if __name__ == "__main__":
    main()
