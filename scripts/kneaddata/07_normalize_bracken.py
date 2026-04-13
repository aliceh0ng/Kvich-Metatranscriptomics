"""
Normalize Bracken abundance estimates following Kvich et al.:
  - Scaling factor = min_reads_across_samples / reads_in_sample
  - Scaled counts = bracken new_est_reads * scaling_factor
  - Threshold: remove taxa with log10(scaled_count) < 0.9 across all samples

Outputs:
  output/bracken/counts_unscaled.tsv   -- raw bracken new_est_reads (for within-sample comparisons)
  output/bracken/counts_scaled.tsv     -- depth-normalized counts (for across-sample comparisons)
  output/bracken/counts_filtered.tsv   -- scaled counts with low-abundance taxa removed (threshold = 0.9)

Usage:
  python scripts/normalize_bracken.py
  
Credit: This file was partially generated using using Claude 3.5 Sonnet by Anthropic - AH
"""

import os
import glob
import pandas as pd
import numpy as np

PROJ        = "/scratch/st-ctropini-1/ahong/bfrag4"
BRACKEN_DIR = os.path.join(PROJ, "output/bracken")
KD_DIR      = os.path.join(PROJ, "kneaddata_output")
SAMPLES_FILE = os.path.join(PROJ, "samples_filtered.txt")
THRESHOLD   = 0.9   # log10 scaled count cutoff from Kvich et al. Fig S1

# ── 1. Load filtered sample list ───────────────────────────────────────────────

with open(SAMPLES_FILE) as f:
    samples_filtered = [line.strip() for line in f if line.strip()]

# ── 2. Load all Bracken outputs ────────────────────────────────────────────────

frames = {}
for sample in samples_filtered:
    f = os.path.join(BRACKEN_DIR, f"{sample}.bracken")
    if not os.path.exists(f):
        print(f"WARNING: missing {f}, skipping")
        continue
    df = pd.read_csv(f, sep="\t")
    frames[sample] = df.set_index("name")["new_est_reads"]

counts = pd.DataFrame(frames).fillna(0)   # taxa x samples
samples = counts.columns.tolist()

# ── 3. Calculate read counts per sample (from KneadData paired R1 output) ──────
# Uses paired_1 (R1) only — total fragments = R1 read count

read_counts = {}
for sample in samples:
    r1 = os.path.join(KD_DIR, sample, f"{sample}_unmapped_R1_kneaddata_paired_1.fastq")
    result = os.popen(f"wc -l < {r1}").read().strip()
    read_counts[sample] = int(result) // 4

print("Read counts per sample (post-KneadData):")
for s, n in sorted(read_counts.items()):
    print(f"  {s}: {n:,}")

# ── 4. Scaling factors ─────────────────────────────────────────────────────────

min_reads = min(read_counts.values())
min_sample = [s for s, n in read_counts.items() if n == min_reads][0]
print(f"\nMinimum reads: {min_reads:,} ({min_sample})")

scaling_factors = {s: min_reads / read_counts[s] for s in samples}

# ── 5. Apply scaling ───────────────────────────────────────────────────────────

scaled = counts.copy()
for sample in samples:
    scaled[sample] = counts[sample] * scaling_factors[sample]

# ── 6. Apply threshold (log10 scaled count >= 0.9 in at least one sample) ──────

log10_scaled = np.log10(scaled + 1)   # +1 to avoid log10(0)
keep = (log10_scaled >= THRESHOLD).any(axis=1)
filtered = scaled[keep]

print(f"\nTaxa before threshold: {len(counts)}")
print(f"Taxa after threshold:  {len(filtered)}")

# ── 7. Save outputs ────────────────────────────────────────────────────────────

counts.to_csv(os.path.join(BRACKEN_DIR, "counts_unscaled.tsv"), sep="\t")
scaled.to_csv(os.path.join(BRACKEN_DIR, "counts_scaled.tsv"), sep="\t")
filtered.to_csv(os.path.join(BRACKEN_DIR, "counts_filtered.tsv"), sep="\t")

sf_df = pd.DataFrame({
    "sample": list(scaling_factors.keys()),
    "total_reads": [read_counts[s] for s in scaling_factors],
    "scaling_factor": list(scaling_factors.values()),
})
sf_df.to_csv(os.path.join(BRACKEN_DIR, "scaling_factors.tsv"), sep="\t", index=False)

print("\nOutput files:")
print(f"  {BRACKEN_DIR}/counts_unscaled.tsv  -- raw bracken counts")
print(f"  {BRACKEN_DIR}/counts_scaled.tsv    -- depth-normalized counts")
print(f"  {BRACKEN_DIR}/counts_filtered.tsv  -- scaled, threshold applied (log10 >= {THRESHOLD})")
print(f"  {BRACKEN_DIR}/scaling_factors.tsv  -- per-sample scaling factors")
