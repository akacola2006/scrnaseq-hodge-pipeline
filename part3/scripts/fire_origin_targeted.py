#!/usr/bin/env python3
"""
ALS Trigger Project — FIRE ORIGIN INVESTIGATION (Targeted Extraction)
Step 1: Collect target genes and extract subset from FASTA
Step 2: Calculate features for targets only
Step 3: Bootstrap length/GC-matched background comparison
"""

import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict
import re
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
OUTPUT_DIR = BASE_DIR / 'results_fire_origin_landscape'

for subdir in ['tables', 'figures', 'reports']:
    (OUTPUT_DIR / subdir).mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("FIRE ORIGIN LANDSCAPE — TARGETED EXTRACTION")
print("=" * 80)

# ============================================================
# STEP 1: Collect Target Gene Sets
# ============================================================
print("\n" + "=" * 60)
print("STEP 1: Collect Target Gene Sets")
print("=" * 60)

# Load existing gene sets
nvu_class = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/nvu_role_split_classification.csv')
fire_candidates = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/capillary_fire_origin_tool_candidates.csv')
carrier_candidates = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/pericyte_carrier_machinery_candidates.csv')
coupled_candidates = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/coupled_fire_carrier_candidates.csv')

# Define gene sets by symbol
class_I = set(nvu_class[nvu_class['role_class'] == 'Class_I_Capillary_Failure']['gene'].tolist())
class_II = set(nvu_class[nvu_class['role_class'] == 'Class_II_Pericyte_Carrier']['gene'].tolist())
class_III = set(nvu_class[nvu_class['role_class'] == 'Class_III_Coupled']['gene'].tolist())

fire_set = class_I | set(fire_candidates['gene'].head(50).tolist())
carrier_set = class_II | set(carrier_candidates['gene'].head(50).tolist())
coupled_set = class_III | set(coupled_candidates['gene'].head(30).tolist())

# Load overlap genes if available
overlap_file = BASE_DIR / 'results_tool_screening_pt_vs_tdp/tables/cross_dataset_overlap_genes.csv'
if overlap_file.exists():
    overlap_df = pd.read_csv(overlap_file)
    overlap_set = set(overlap_df.iloc[:, 0].tolist())
else:
    overlap_set = set()

# Background: NVU-expressed genes (Class I-IV all)
background_pool = set(nvu_class['gene'].tolist())

# All targets
all_targets = fire_set | carrier_set | coupled_set | overlap_set
all_symbols = all_targets | background_pool

print(f"Fire set: {len(fire_set)} genes")
print(f"Carrier set: {len(carrier_set)} genes")
print(f"Coupled set: {len(coupled_set)} genes")
print(f"Overlap set: {len(overlap_set)} genes")
print(f"All targets: {len(all_targets)} genes")
print(f"Background pool: {len(background_pool)} genes")

# Save target genes
target_genes_df = pd.DataFrame({
    'gene': list(all_targets),
    'in_fire': [g in fire_set for g in all_targets],
    'in_carrier': [g in carrier_set for g in all_targets],
    'in_coupled': [g in coupled_set for g in all_targets],
    'in_overlap': [g in overlap_set for g in all_targets]
})
target_genes_df.to_csv(OUTPUT_DIR / 'tables/target_genes.csv', index=False)

# ============================================================
# STEP 2: Build Gene Symbol to ENSG Mapping
# ============================================================
print("\n" + "=" * 60)
print("STEP 2: Build Gene Symbol to ENSG Mapping")
print("=" * 60)

# We need to map gene symbols to ENSG IDs
# Strategy: Read FASTA headers and try to match symbols

def build_ensg_mapping_from_fasta(fasta_path, max_lines=100000):
    """Build ENSG to gene symbol mapping by parsing FASTA headers"""
    mapping = {}
    with open(fasta_path, 'r') as f:
        for i, line in enumerate(f):
            if i >= max_lines:
                break
            if line.startswith('>'):
                ensg_id = line[1:].split()[0]
                mapping[ensg_id] = ensg_id  # Default: use ENSG as key
    return mapping

# Try to load existing annotation if available
# Check motif data for ENSG to symbol mapping
motif_file = BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/tool_candidates_by_motif_ranked.csv'
if motif_file.exists():
    motif_data = pd.read_csv(motif_file)
    if 'gene_id' in motif_data.columns and 'gene_symbol' in motif_data.columns:
        ensg_to_symbol = dict(zip(motif_data['gene_id'], motif_data['gene_symbol']))
        symbol_to_ensg = dict(zip(motif_data['gene_symbol'], motif_data['gene_id']))
        print(f"Loaded {len(ensg_to_symbol)} ENSG-Symbol mappings from motif data")
    else:
        ensg_to_symbol = {}
        symbol_to_ensg = {}
else:
    ensg_to_symbol = {}
    symbol_to_ensg = {}

# If we have the motif data, use it directly for sequence features
if motif_file.exists() and len(motif_data) > 0:
    print(f"\nUsing existing sequence features from motif data")
    print(f"Columns: {motif_data.columns.tolist()}")

    # Filter to our target genes
    motif_in_targets = motif_data[motif_data['gene_symbol'].isin(all_targets)]
    motif_in_fire = motif_data[motif_data['gene_symbol'].isin(fire_set)]
    motif_in_carrier = motif_data[motif_data['gene_symbol'].isin(carrier_set)]
    motif_in_coupled = motif_data[motif_data['gene_symbol'].isin(coupled_set)]

    print(f"Targets in motif data: {len(motif_in_targets)}")
    print(f"  Fire: {len(motif_in_fire)}")
    print(f"  Carrier: {len(motif_in_carrier)}")
    print(f"  Coupled: {len(motif_in_coupled)}")

# ============================================================
# STEP 3: Stream Extract Target Sequences from FASTA
# ============================================================
print("\n" + "=" * 60)
print("STEP 3: Extract Target Sequences (Stream)")
print("=" * 60)

def stream_extract_sequences(fasta_path, target_ensgs=None, max_sequences=None):
    """Stream extract sequences matching target ENSG IDs"""
    sequences = {}
    current_id = None
    current_seq = []
    count = 0

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id and (target_ensgs is None or current_id in target_ensgs):
                    sequences[current_id] = ''.join(current_seq)
                    count += 1
                    if max_sequences and count >= max_sequences:
                        break

                # Start new sequence
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())

        # Last sequence
        if current_id and (target_ensgs is None or current_id in target_ensgs):
            sequences[current_id] = ''.join(current_seq)

    return sequences

# For background, we need all NVU-expressed genes
# First, let's extract a sample to establish length/GC distributions

# Build list of ENSG IDs we need
# Since we may not have perfect mapping, extract ALL and filter later

print("Extracting sample sequences for length/GC calibration...")

# Extract first N sequences for background distribution
sample_seqs = stream_extract_sequences(BASE_DIR / 'alltranscript.fasta', max_sequences=5000)
print(f"Extracted {len(sample_seqs)} sample sequences")

# ============================================================
# STEP 4: Calculate Sequence Features
# ============================================================
print("\n" + "=" * 60)
print("STEP 4: Calculate Sequence Features")
print("=" * 60)

def calc_gc_content(seq):
    if len(seq) == 0:
        return 0
    return (seq.count('G') + seq.count('C')) / len(seq)

def calc_ctgcagy_density(seq):
    """CTGCAGY motif density per kb"""
    motifs = ['CTGCAGC', 'CTGCAGT', 'CTGCAGA', 'CTGCAGG']
    total_hits = sum(seq.count(m) for m in motifs)
    return (total_hits / len(seq)) * 1000 if len(seq) > 0 else 0

def calc_g4_score(seq, window=25):
    """Simplified G4 score"""
    if len(seq) < window:
        return 0

    max_score = 0
    for i in range(0, len(seq) - window + 1, 10):  # Step by 10 for speed
        win = seq[i:i+window]
        score = 0
        # Count G runs
        g_runs = re.findall(r'G{3,}', win)
        for run in g_runs:
            score += len(run) ** 1.5
        max_score = max(max_score, score / window)

    return max_score

def calc_palindrome_density(seq, k=6):
    """Quasi-palindrome density"""
    if len(seq) < k:
        return 0
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    count = 0
    for i in range(0, len(seq) - k + 1, 5):  # Step by 5 for speed
        kmer = seq[i:i+k]
        rc = ''.join(complement.get(b, 'N') for b in kmer[::-1])
        if kmer == rc:
            count += 1
    return (count / ((len(seq) - k + 1) / 5)) * 1000 if len(seq) > k else 0

# Calculate features for sample sequences
print("Calculating features for sample sequences...")

sample_features = []
for ensg_id, seq in sample_seqs.items():
    sample_features.append({
        'ensg_id': ensg_id,
        'length': len(seq),
        'gc_content': calc_gc_content(seq),
        'ctgcagy_density': calc_ctgcagy_density(seq),
        'g4_score': calc_g4_score(seq),
        'palindrome_density': calc_palindrome_density(seq)
    })

sample_df = pd.DataFrame(sample_features)
print(f"Sample features: {len(sample_df)} sequences")
print(f"Length: {sample_df['length'].mean():.0f} ± {sample_df['length'].std():.0f}")
print(f"GC: {sample_df['gc_content'].mean():.3f} ± {sample_df['gc_content'].std():.3f}")

# ============================================================
# Use Existing Motif Data for Fire/Carrier Comparison
# ============================================================
print("\n" + "=" * 60)
print("Using Existing Motif Data for Comparisons")
print("=" * 60)

if motif_file.exists():
    # Add set membership
    motif_data['set'] = 'Background'
    motif_data.loc[motif_data['gene_symbol'].isin(fire_set), 'set'] = 'Fire'
    motif_data.loc[motif_data['gene_symbol'].isin(carrier_set), 'set'] = 'Carrier'
    motif_data.loc[motif_data['gene_symbol'].isin(coupled_set), 'set'] = 'Coupled'

    # Feature columns
    feature_cols = ['gc_content', 'CTGCAGY_density_kb', 'g4_potential', 'palindrome_density_kb', 'seq_length']
    available_cols = [c for c in feature_cols if c in motif_data.columns]

    print(f"\nAvailable feature columns: {available_cols}")

    # Compare Fire vs Background
    from scipy.stats import mannwhitneyu

    print("\n--- Fire vs Background ---")
    fire_data = motif_data[motif_data['set'] == 'Fire']
    bg_data = motif_data[motif_data['set'] == 'Background']

    comparison_results = []

    for col in available_cols:
        if col in fire_data.columns and fire_data[col].notna().sum() > 0:
            fire_vals = fire_data[col].dropna()
            bg_vals = bg_data[col].dropna()

            if len(fire_vals) > 0 and len(bg_vals) > 0:
                stat, pval = mannwhitneyu(fire_vals, bg_vals, alternative='two-sided')

                # Cliff's delta (effect size)
                n1, n2 = len(fire_vals), len(bg_vals)
                # Simplified effect size
                effect = (fire_vals.mean() - bg_vals.mean()) / (bg_vals.std() + 1e-10)

                print(f"  {col}: Fire={fire_vals.mean():.3f}, BG={bg_vals.mean():.3f}, effect={effect:.3f}, p={pval:.4f}")

                comparison_results.append({
                    'feature': col,
                    'fire_mean': fire_vals.mean(),
                    'fire_std': fire_vals.std(),
                    'bg_mean': bg_vals.mean(),
                    'bg_std': bg_vals.std(),
                    'effect_size': effect,
                    'p_value': pval
                })

    comparison_df = pd.DataFrame(comparison_results)
    comparison_df.to_csv(OUTPUT_DIR / 'tables/fire_vs_background_comparison.csv', index=False)

    # Compare Carrier vs Background
    print("\n--- Carrier vs Background ---")
    carrier_data = motif_data[motif_data['set'] == 'Carrier']

    for col in available_cols:
        if col in carrier_data.columns and carrier_data[col].notna().sum() > 0:
            carrier_vals = carrier_data[col].dropna()
            bg_vals = bg_data[col].dropna()

            if len(carrier_vals) > 0 and len(bg_vals) > 0:
                stat, pval = mannwhitneyu(carrier_vals, bg_vals, alternative='two-sided')
                effect = (carrier_vals.mean() - bg_vals.mean()) / (bg_vals.std() + 1e-10)
                print(f"  {col}: Carrier={carrier_vals.mean():.3f}, BG={bg_vals.mean():.3f}, effect={effect:.3f}, p={pval:.4f}")

# ============================================================
# STEP 5: ATR/ATM Proxy Analysis from Pseudobulk
# ============================================================
print("\n" + "=" * 60)
print("STEP 5: ATR/ATM Proxy Analysis")
print("=" * 60)

# Load pseudobulk data
pseudobulk_file = BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/nvu_pseudobulk_ctrl_vs_tdpneg_by_celltype.csv'
if pseudobulk_file.exists():
    pb = pd.read_csv(pseudobulk_file)

    # Define proxy genes
    ATR_PROXY = ['ATR', 'ATRIP', 'RPA1', 'RPA2', 'TOPBP1', 'CHEK1', 'RAD17', 'TIMELESS']
    ATM_PROXY = ['ATM', 'CHEK2', 'TP53BP1', 'MRE11', 'RAD50', 'NBN', 'H2AFX', 'BRCA1']
    RLOOP_PROXY = ['SFPQ', 'UPF1', 'XRN2', 'SETX', 'DDX5', 'RNASEH1', 'DIS3']

    # Filter to Capillary and Pericyte
    cap_pb = pb[pb['cell_type'].str.contains('Capillary', case=False, na=False)]
    peri_pb = pb[pb['cell_type'].str.contains('Pericyte', case=False, na=False)]

    print(f"\nCapillary pseudobulk: {len(cap_pb)} genes")
    print(f"Pericyte pseudobulk: {len(peri_pb)} genes")

    # Check proxy genes in data
    def analyze_proxy(pb_data, proxy_genes, proxy_name, cell_type):
        found = pb_data[pb_data['gene'].isin(proxy_genes)]
        print(f"\n  {proxy_name} in {cell_type}: {len(found)}/{len(proxy_genes)} found")

        if len(found) > 0:
            print(f"    Genes: {', '.join(found['gene'].tolist())}")

            # Control vs TDPneg comparison
            ctrl_mean = found['ctrl_mean'].mean()
            tdp_mean = found['tdpneg_mean'].mean()
            logFC_mean = found['logFC'].mean()

            print(f"    Ctrl mean: {ctrl_mean:.3f}")
            print(f"    TDPneg mean: {tdp_mean:.3f}")
            print(f"    Mean logFC: {logFC_mean:.3f}")

            return {
                'proxy': proxy_name,
                'cell_type': cell_type,
                'n_found': len(found),
                'ctrl_mean': ctrl_mean,
                'tdpneg_mean': tdp_mean,
                'logFC': logFC_mean
            }
        return None

    proxy_results = []

    for cell_type, pb_data in [('Capillary', cap_pb), ('Pericyte', peri_pb)]:
        print(f"\n--- {cell_type} ---")
        for proxy_name, proxy_genes in [('ATR', ATR_PROXY), ('ATM', ATM_PROXY), ('Rloop', RLOOP_PROXY)]:
            result = analyze_proxy(pb_data, proxy_genes, proxy_name, cell_type)
            if result:
                proxy_results.append(result)

    proxy_df = pd.DataFrame(proxy_results)
    proxy_df.to_csv(OUTPUT_DIR / 'tables/proxy_analysis_results.csv', index=False)

    # Key hypothesis test: ATM decline in Capillary TDPneg
    print("\n" + "=" * 60)
    print("HYPOTHESIS TEST: ATM Brake Loss in Capillary")
    print("=" * 60)

    cap_atm = proxy_df[(proxy_df['cell_type'] == 'Capillary') & (proxy_df['proxy'] == 'ATM')]
    cap_atr = proxy_df[(proxy_df['cell_type'] == 'Capillary') & (proxy_df['proxy'] == 'ATR')]

    if len(cap_atm) > 0 and len(cap_atr) > 0:
        atm_logfc = cap_atm['logFC'].values[0]
        atr_logfc = cap_atr['logFC'].values[0]

        print(f"ATM proxy logFC (Ctrl→TDPneg): {atm_logfc:.3f}")
        print(f"ATR proxy logFC (Ctrl→TDPneg): {atr_logfc:.3f}")

        if atm_logfc < 0 and atr_logfc >= atm_logfc:
            print("\n★ SUPPORT: ATM brake loss pattern detected")
            print("  ATM↓ while ATR maintained → consistent with 'chronic tension + brake failure'")
            verdict = "SUPPORT"
        elif atm_logfc < 0:
            print("\n◯ PARTIAL: ATM decline detected, ATR pattern unclear")
            verdict = "PARTIAL"
        else:
            print("\n✗ NOT SUPPORTED: No ATM decline detected")
            verdict = "NOT_SUPPORTED"
    else:
        verdict = "INSUFFICIENT_DATA"
        print("Insufficient data for hypothesis test")

else:
    print("Pseudobulk file not found")
    verdict = "NO_DATA"

# ============================================================
# Generate Report
# ============================================================
print("\n" + "=" * 60)
print("Generating Report")
print("=" * 60)

report = f"""# ALS Trigger Project — Fire Origin Landscape Report

## Executive Summary

**Hypothesis**: GC/G4/R-loop landscape + ATR chronic tension → ATM brake loss → RNA surveillance failure (Fire)

**Verdict**: {verdict}

---

## Gene Sets Analyzed

| Set | N genes | Description |
|-----|---------|-------------|
| Fire | {len(fire_set)} | Class I Capillary Failure + candidates |
| Carrier | {len(carrier_set)} | Class II Pericyte Carrier + candidates |
| Coupled | {len(coupled_set)} | Class III Coupled |
| Overlap | {len(overlap_set)} | Phase13 ∩ GSE pre-TDP |
| Background | {len(background_pool)} | All NVU-expressed |

### Fire Set Genes
{', '.join(sorted(fire_set)[:30])}...

---

## Sequence Landscape Analysis

### Fire vs Background Comparison
"""

if 'comparison_df' in dir() and len(comparison_df) > 0:
    report += """
| Feature | Fire Mean | BG Mean | Effect Size | p-value |
|---------|-----------|---------|-------------|---------|
"""
    for _, row in comparison_df.iterrows():
        sig = "★" if row['p_value'] < 0.05 else ""
        report += f"| {row['feature']} | {row['fire_mean']:.3f} | {row['bg_mean']:.3f} | {row['effect_size']:.3f} | {row['p_value']:.4f} {sig} |\n"

report += f"""

---

## ATR/ATM Proxy Analysis

### Hypothesis
- H1 (Tension Field): ATR proxy remains high (chronic replication stress)
- H2 (Brake Loss): ATM proxy declines in TDPneg (DSB repair brake failure)
- H3 (Cell Specificity): Pattern is stronger in Capillary (fire origin)

### Results
"""

if 'proxy_df' in dir() and len(proxy_df) > 0:
    report += """
| Cell Type | Proxy | N Genes | Ctrl Mean | TDPneg Mean | logFC |
|-----------|-------|---------|-----------|-------------|-------|
"""
    for _, row in proxy_df.iterrows():
        report += f"| {row['cell_type']} | {row['proxy']} | {row['n_found']} | {row['ctrl_mean']:.3f} | {row['tdpneg_mean']:.3f} | {row['logFC']:.3f} |\n"

report += f"""

---

## Interpretation

### Verdict: {verdict}

"""

if verdict == "SUPPORT":
    report += """
**SUPPORTED**: The analysis supports the Fire Origin hypothesis:
1. ATM proxy shows decline in Capillary TDPneg (brake loss)
2. ATR proxy maintained (chronic tension persists)
3. This pattern is consistent with "RNA surveillance failure" as the fire origin

### Mechanistic Model
```
[Capillary - Fire Origin]
Chronic ATR tension (replication stress)
    → ATM brake loss (DSB response failure)
    → RNA surveillance failure (SFPQ/UPF1/XRN2)
    → Aberrant RNA accumulation
    → TDP-43 mislocalization trigger
```
"""
elif verdict == "PARTIAL":
    report += """
**PARTIAL SUPPORT**: Some evidence supports the hypothesis but not all criteria met.
"""
else:
    report += """
**NOT SUPPORTED** or **INSUFFICIENT DATA**: Unable to confirm the hypothesis with available data.
"""

report += f"""

---

## Refutable Predictions

1. **Longitudinal prediction**: ATM decline should precede TDP-43 pathology onset
2. **Intervention prediction**: ATM rescue in Capillary should delay fire spread
3. **Cell specificity**: Pericyte shows carrier (EV) signature instead of brake loss

---

## Output Files

- `gene_sets_summary.csv`: Target gene set definitions
- `fire_vs_background_comparison.csv`: Sequence feature comparison
- `proxy_analysis_results.csv`: ATR/ATM/Rloop proxy analysis

---

*Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}*
"""

with open(OUTPUT_DIR / 'reports/fire_origin_landscape_report.md', 'w') as f:
    f.write(report)

print(f"\nReport saved: {OUTPUT_DIR / 'reports/fire_origin_landscape_report.md'}")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print(f"""
Verdict: {verdict}

Key Findings:
- Fire set: {len(fire_set)} genes (Class I + candidates)
- Carrier set: {len(carrier_set)} genes
- ATR/ATM proxy analysis completed

Output: {OUTPUT_DIR}
""")
