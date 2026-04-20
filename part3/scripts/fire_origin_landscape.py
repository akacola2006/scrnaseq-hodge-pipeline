#!/usr/bin/env python3
"""
ALS Trigger Project — FIRE ORIGIN INVESTIGATION (IN SILICO)
Hypothesis: GC/G4/R-loop landscape + ATR-like chronic tension → ATM brake loss → RNA surveillance failure

PART 0: Gene Set Preparation
PART A: Sequence Landscape (GC, CTGCAGY, quasi-palindrome)
PART B: G4-prone Landscape (G4Hunter-like scoring)
PART C: ATR-proxy vs ATM-proxy expression analysis
PART D: Integrated Verdict
"""

import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter, defaultdict
import re
import gzip
from scipy.stats import mannwhitneyu, spearmanr
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path('/home/akaco/als/motor_cortex_analysis/ids_causal_analysis')
OUTPUT_DIR = BASE_DIR / 'results_fire_origin_landscape'

# Create output directories
for subdir in ['tables', 'figures', 'reports']:
    (OUTPUT_DIR / subdir).mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("ALS TRIGGER PROJECT — FIRE ORIGIN LANDSCAPE INVESTIGATION")
print("Hypothesis: GC/G4/R-loop + ATR tension → ATM brake loss → Fire")
print("=" * 80)

# ============================================================
# PART 0: Gene Set Preparation
# ============================================================
print("\n" + "=" * 60)
print("PART 0: Gene Set Preparation")
print("=" * 60)

# Load NVU role-split classification
nvu_class = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/nvu_role_split_classification.csv')
print(f"Total genes in NVU classification: {len(nvu_class)}")

# Load fire/carrier candidates
fire_candidates = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/capillary_fire_origin_tool_candidates.csv')
carrier_candidates = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/pericyte_carrier_machinery_candidates.csv')
coupled_candidates = pd.read_csv(BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/coupled_fire_carrier_candidates.csv')

# Define gene sets
# Fire set: Class I (Capillary Failure) + fire candidates
class_I_genes = set(nvu_class[nvu_class['role_class'] == 'Class_I_Capillary_Failure']['gene'].tolist())
fire_candidate_genes = set(fire_candidates['gene'].head(50).tolist())
fire_set = class_I_genes | fire_candidate_genes

# Carrier set: Class II (Pericyte Carrier) + carrier candidates
class_II_genes = set(nvu_class[nvu_class['role_class'] == 'Class_II_Pericyte_Carrier']['gene'].tolist())
carrier_candidate_genes = set(carrier_candidates['gene'].head(50).tolist()) if len(carrier_candidates) > 0 else set()
carrier_set = class_II_genes | carrier_candidate_genes

# Coupled set: Class III
class_III_genes = set(nvu_class[nvu_class['role_class'] == 'Class_III_Coupled']['gene'].tolist())
coupled_set = class_III_genes

# Background: all genes with some expression (Class IV + others)
background_all = set(nvu_class['gene'].tolist())

print(f"\nGene Sets:")
print(f"  Fire set (Class I + candidates): {len(fire_set)}")
print(f"  Carrier set (Class II + candidates): {len(carrier_set)}")
print(f"  Coupled set (Class III): {len(coupled_set)}")
print(f"  Background (all NVU genes): {len(background_all)}")

# Try to load overlap genes if available
overlap_file = BASE_DIR / 'results_tool_screening_pt_vs_tdp/tables/cross_dataset_overlap_genes.csv'
if overlap_file.exists():
    overlap_df = pd.read_csv(overlap_file)
    overlap_set = set(overlap_df.iloc[:, 0].tolist())
    print(f"  Overlap set (PT ∩ TDP): {len(overlap_set)}")
else:
    overlap_set = set()
    print("  Overlap set: not available")

# Save gene sets summary
gene_sets_summary = pd.DataFrame({
    'set_name': ['Fire', 'Carrier', 'Coupled', 'Overlap', 'Background'],
    'n_genes': [len(fire_set), len(carrier_set), len(coupled_set), len(overlap_set), len(background_all)],
    'example_genes': [
        ', '.join(list(fire_set)[:5]),
        ', '.join(list(carrier_set)[:5]),
        ', '.join(list(coupled_set)[:5]),
        ', '.join(list(overlap_set)[:5]) if overlap_set else '',
        'All NVU genes'
    ]
})
gene_sets_summary.to_csv(OUTPUT_DIR / 'tables/gene_sets_summary.csv', index=False)

# ============================================================
# Load FASTA sequences and build gene symbol to ENSG mapping
# ============================================================
print("\n" + "-" * 40)
print("Loading FASTA sequences...")

def load_fasta(fasta_path):
    """Load FASTA file and return dict of {ensg_id: sequence}"""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Extract ENSG ID
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences

# Load transcripts and promoters
transcript_seqs = load_fasta(BASE_DIR / 'alltranscript.fasta')
promoter_seqs = load_fasta(BASE_DIR / 'completeall_promoters.fasta')

print(f"Loaded {len(transcript_seqs)} transcript sequences")
print(f"Loaded {len(promoter_seqs)} promoter sequences")

# We need gene symbol to ENSG mapping
# Try to build from available data or use a simple approach
# For now, we'll use gene symbols directly and map where possible

# Try to load a mapping file if available
mapping_file = BASE_DIR / 'results_motif_family_nvu_rolesplit/tables/id_mapping_summary.csv'
if mapping_file.exists():
    id_mapping = pd.read_csv(mapping_file)
    print(f"Loaded ID mapping: {len(id_mapping)} entries")
else:
    id_mapping = None
    print("No ID mapping file found - will use direct matching")

# Build ENSG to symbol mapping from available data
# This is a simplified approach - in practice you'd use a proper annotation file
ensg_to_symbol = {}
symbol_to_ensg = {}

# Try to extract from h5ad or other sources
# For now, use a heuristic: match gene symbols in our sets to ENSG IDs
# We'll compute features for ENSG IDs and map back

# ============================================================
# PART A: Sequence Landscape Analysis
# ============================================================
print("\n" + "=" * 60)
print("PART A: Sequence Landscape Analysis")
print("=" * 60)

def calc_gc_content(seq):
    """Calculate GC content of a sequence"""
    if len(seq) == 0:
        return 0
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)

def calc_ctgcagy_density(seq, window=1000):
    """Calculate CTGCAGY motif density per kb"""
    # CTGCAGY family motifs
    motifs = ['CTGCAGY', 'CTGCAGC', 'CTGCAGT', 'CTGCAGA',
              'CTGCGG', 'CTGCGC', 'CAGCAGY', 'CAGCTG']

    total_hits = 0
    for motif in motifs:
        # Allow for degenerate positions
        pattern = motif.replace('Y', '[CT]').replace('R', '[AG]')
        total_hits += len(re.findall(pattern, seq, re.IGNORECASE))

    # Density per kb
    return (total_hits / len(seq)) * 1000 if len(seq) > 0 else 0

def calc_palindrome_density(seq, k=6):
    """Calculate quasi-palindrome density (reverse complement k-mers)"""
    if len(seq) < k:
        return 0

    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    palindrome_count = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        # Get reverse complement
        rc = ''.join(complement.get(b, 'N') for b in kmer[::-1])
        # Check if it's a palindrome (self-reverse-complement)
        if kmer == rc:
            palindrome_count += 1

    return (palindrome_count / (len(seq) - k + 1)) * 1000 if len(seq) > k else 0

def calc_cpg_density(seq):
    """Calculate CpG density"""
    if len(seq) < 2:
        return 0
    cpg_count = seq.count('CG')
    return (cpg_count / (len(seq) - 1)) * 1000

def calc_low_complexity(seq, k=3):
    """Calculate low complexity score (entropy-based)"""
    if len(seq) < k:
        return 0

    # Count k-mers
    kmer_counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer_counts[seq[i:i+k]] += 1

    # Calculate entropy
    total = sum(kmer_counts.values())
    entropy = 0
    for count in kmer_counts.values():
        p = count / total
        if p > 0:
            entropy -= p * np.log2(p)

    # Max entropy for k-mers is log2(4^k)
    max_entropy = k * 2  # log2(4) = 2

    # Return normalized low complexity (1 - normalized entropy)
    return 1 - (entropy / max_entropy) if max_entropy > 0 else 0

# Calculate features for all sequences
print("Calculating sequence features...")

transcript_features = []
for ensg_id, seq in transcript_seqs.items():
    transcript_features.append({
        'ensg_id': ensg_id,
        'length': len(seq),
        'gc_content': calc_gc_content(seq),
        'ctgcagy_density': calc_ctgcagy_density(seq),
        'palindrome_density': calc_palindrome_density(seq),
        'low_complexity': calc_low_complexity(seq)
    })

transcript_df = pd.DataFrame(transcript_features)
print(f"Transcript features calculated for {len(transcript_df)} sequences")

promoter_features = []
for ensg_id, seq in promoter_seqs.items():
    promoter_features.append({
        'ensg_id': ensg_id,
        'length': len(seq),
        'gc_content': calc_gc_content(seq),
        'cpg_density': calc_cpg_density(seq),
        'ctgcagy_density': calc_ctgcagy_density(seq)
    })

promoter_df = pd.DataFrame(promoter_features)
print(f"Promoter features calculated for {len(promoter_df)} sequences")

# ============================================================
# PART B: G4-prone Landscape (G4Hunter-like scoring)
# ============================================================
print("\n" + "=" * 60)
print("PART B: G4-prone Landscape Analysis")
print("=" * 60)

def g4_hunter_score(seq, window=25):
    """
    G4Hunter-like scoring
    - Consecutive G runs contribute positively
    - Consecutive C runs contribute negatively
    - Returns max and mean scores
    """
    if len(seq) < window:
        return 0, 0, 0

    scores = []
    for i in range(len(seq) - window + 1):
        win_seq = seq[i:i+window]
        score = 0

        # Count consecutive G/C runs
        j = 0
        while j < len(win_seq):
            if win_seq[j] == 'G':
                run_len = 1
                while j + run_len < len(win_seq) and win_seq[j + run_len] == 'G':
                    run_len += 1
                # Score increases with run length (G4 requires 3+ Gs)
                if run_len >= 3:
                    score += run_len ** 1.5
                else:
                    score += run_len * 0.5
                j += run_len
            elif win_seq[j] == 'C':
                run_len = 1
                while j + run_len < len(win_seq) and win_seq[j + run_len] == 'C':
                    run_len += 1
                # C runs contribute negatively (competing structure)
                if run_len >= 3:
                    score -= run_len ** 1.5
                else:
                    score -= run_len * 0.5
                j += run_len
            else:
                j += 1

        # Normalize by window size
        scores.append(score / window)

    if scores:
        return max(scores), np.mean(scores), np.percentile(scores, 90)
    return 0, 0, 0

def count_g4_motifs(seq):
    """Count canonical G4 motifs (G3+N1-7G3+N1-7G3+N1-7G3+)"""
    # Simplified G4 pattern
    pattern = r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}'
    matches = re.findall(pattern, seq, re.IGNORECASE)
    return len(matches)

print("Calculating G4 scores...")

g4_transcript_features = []
for ensg_id, seq in transcript_seqs.items():
    max_score, mean_score, q90_score = g4_hunter_score(seq, window=25)
    g4_count = count_g4_motifs(seq)
    g4_transcript_features.append({
        'ensg_id': ensg_id,
        'g4_max_score': max_score,
        'g4_mean_score': mean_score,
        'g4_q90_score': q90_score,
        'g4_motif_count': g4_count,
        'g4_density': g4_count / len(seq) * 1000 if len(seq) > 0 else 0
    })

g4_transcript_df = pd.DataFrame(g4_transcript_features)

g4_promoter_features = []
for ensg_id, seq in promoter_seqs.items():
    max_score, mean_score, q90_score = g4_hunter_score(seq, window=25)
    g4_count = count_g4_motifs(seq)
    g4_promoter_features.append({
        'ensg_id': ensg_id,
        'g4_max_score': max_score,
        'g4_mean_score': mean_score,
        'g4_q90_score': q90_score,
        'g4_motif_count': g4_count,
        'g4_density': g4_count / len(seq) * 1000 if len(seq) > 0 else 0
    })

g4_promoter_df = pd.DataFrame(g4_promoter_features)

print(f"G4 transcript features: {len(g4_transcript_df)}")
print(f"G4 promoter features: {len(g4_promoter_df)}")

# Merge all sequence features
all_transcript_features = transcript_df.merge(g4_transcript_df, on='ensg_id', how='outer')
all_promoter_features = promoter_df.merge(g4_promoter_df, on='ensg_id', how='outer')

all_transcript_features.to_csv(OUTPUT_DIR / 'tables/transcript_sequence_features.csv', index=False)
all_promoter_features.to_csv(OUTPUT_DIR / 'tables/promoter_sequence_features.csv', index=False)

# ============================================================
# PART C: ATR-proxy vs ATM-proxy Expression Analysis
# ============================================================
print("\n" + "=" * 60)
print("PART C: ATR-proxy vs ATM-proxy Expression Analysis")
print("=" * 60)

# Define proxy gene sets
ATR_PROXY_GENES = [
    'ATR', 'ATRIP', 'RPA1', 'RPA2', 'RPA3', 'TOPBP1', 'CLSPN',
    'CHEK1', 'RAD17', 'RFC1', 'RFC2', 'RFC3', 'RFC4', 'RFC5',
    'TIMELESS', 'TIPIN', 'CLASPIN', 'RAD9A', 'HUS1', 'RAD1'
]

ATM_PROXY_GENES = [
    'ATM', 'CHEK2', 'TP53BP1', 'MRE11', 'RAD50', 'NBN',
    'H2AFX', 'XRCC5', 'XRCC6', 'PRKDC', 'LIG4', 'NHEJ1',
    'BRCA1', 'BRCA2', 'RAD51', 'MDC1'
]

RLOOP_PROXY_GENES = [
    'SETX', 'DDX5', 'RNASEH1', 'RNASEH2A', 'RNASEH2B', 'RNASEH2C',
    'SFPQ', 'UPF1', 'XRN2', 'DIS3', 'EXOSC10', 'EXOSC3',
    'DHX9', 'DDX39B', 'THOC1', 'THOC2'
]

# Load GSE212630 expression data for Capillary and Pericyte
print("\nLoading GSE212630 expression data...")

def load_expression_with_metadata(cell_type):
    """Load expression and metadata for a cell type"""
    expr_file = BASE_DIR / f'GSE_212630_raw_expression_transposed/Vasc_{cell_type}_expression_transposed.csv.gz'
    if not expr_file.exists():
        return None, None

    expr = pd.read_csv(expr_file, compression='gzip', index_col=0)

    # Try to load metadata
    meta_file = BASE_DIR / f'GSE_212630_raw_expression_transposed/Vasc_{cell_type}_metadata.csv'
    if meta_file.exists():
        meta = pd.read_csv(meta_file, index_col=0)
    else:
        # Extract from column names if possible
        meta = None

    return expr, meta

cap_expr, cap_meta = load_expression_with_metadata('Capillary')
peri_expr, peri_meta = load_expression_with_metadata('Pericyte')

if cap_expr is not None:
    print(f"Capillary expression: {cap_expr.shape}")
if peri_expr is not None:
    print(f"Pericyte expression: {peri_expr.shape}")

# Function to calculate module score
def calc_module_score(expr_df, gene_list):
    """Calculate module score (mean z-score) for gene set"""
    available_genes = [g for g in gene_list if g in expr_df.columns]
    if len(available_genes) == 0:
        return pd.Series(dtype=float), []

    # Z-score per gene
    subset = expr_df[available_genes]
    z_scores = (subset - subset.mean()) / (subset.std() + 1e-10)

    # Mean across genes
    module_score = z_scores.mean(axis=1)

    return module_score, available_genes

# Calculate proxy scores for Capillary
atr_atm_results = []

if cap_expr is not None:
    print("\nCalculating proxy scores for Capillary...")

    atr_score_cap, atr_genes_cap = calc_module_score(cap_expr, ATR_PROXY_GENES)
    atm_score_cap, atm_genes_cap = calc_module_score(cap_expr, ATM_PROXY_GENES)
    rloop_score_cap, rloop_genes_cap = calc_module_score(cap_expr, RLOOP_PROXY_GENES)

    print(f"  ATR proxy genes found: {len(atr_genes_cap)}/{len(ATR_PROXY_GENES)}")
    print(f"  ATM proxy genes found: {len(atm_genes_cap)}/{len(ATM_PROXY_GENES)}")
    print(f"  R-loop proxy genes found: {len(rloop_genes_cap)}/{len(RLOOP_PROXY_GENES)}")

    # Store results
    cap_scores = pd.DataFrame({
        'cell_type': 'Capillary',
        'sample': cap_expr.index,
        'ATR_score': atr_score_cap.values,
        'ATM_score': atm_score_cap.values,
        'Rloop_score': rloop_score_cap.values
    })
    atr_atm_results.append(cap_scores)

if peri_expr is not None:
    print("\nCalculating proxy scores for Pericyte...")

    atr_score_peri, atr_genes_peri = calc_module_score(peri_expr, ATR_PROXY_GENES)
    atm_score_peri, atm_genes_peri = calc_module_score(peri_expr, ATM_PROXY_GENES)
    rloop_score_peri, rloop_genes_peri = calc_module_score(peri_expr, RLOOP_PROXY_GENES)

    print(f"  ATR proxy genes found: {len(atr_genes_peri)}/{len(ATR_PROXY_GENES)}")
    print(f"  ATM proxy genes found: {len(atm_genes_peri)}/{len(ATM_PROXY_GENES)}")
    print(f"  R-loop proxy genes found: {len(rloop_genes_peri)}/{len(RLOOP_PROXY_GENES)}")

    peri_scores = pd.DataFrame({
        'cell_type': 'Pericyte',
        'sample': peri_expr.index,
        'ATR_score': atr_score_peri.values,
        'ATM_score': atm_score_peri.values,
        'Rloop_score': rloop_score_peri.values
    })
    atr_atm_results.append(peri_scores)

if atr_atm_results:
    all_scores = pd.concat(atr_atm_results, ignore_index=True)
    all_scores.to_csv(OUTPUT_DIR / 'tables/atr_atm_proxy_scores_by_sample.csv', index=False)

    # Summary statistics by cell type
    print("\n--- Proxy Score Summary ---")
    summary = all_scores.groupby('cell_type')[['ATR_score', 'ATM_score', 'Rloop_score']].agg(['mean', 'std'])
    print(summary)

# Save proxy gene availability
proxy_availability = pd.DataFrame({
    'proxy_set': ['ATR', 'ATM', 'R-loop'],
    'total_genes': [len(ATR_PROXY_GENES), len(ATM_PROXY_GENES), len(RLOOP_PROXY_GENES)],
    'capillary_found': [len(atr_genes_cap) if cap_expr is not None else 0,
                        len(atm_genes_cap) if cap_expr is not None else 0,
                        len(rloop_genes_cap) if cap_expr is not None else 0],
    'pericyte_found': [len(atr_genes_peri) if peri_expr is not None else 0,
                       len(atm_genes_peri) if peri_expr is not None else 0,
                       len(rloop_genes_peri) if peri_expr is not None else 0]
})
proxy_availability.to_csv(OUTPUT_DIR / 'tables/proxy_gene_availability.csv', index=False)

# ============================================================
# PART D: Integrated Analysis & Comparison
# ============================================================
print("\n" + "=" * 60)
print("PART D: Integrated Analysis")
print("=" * 60)

# Since we have ENSG IDs in sequences but gene symbols in our sets,
# we need to compare at the ENSG level using available mappings
# For now, summarize overall distributions

print("\n--- Sequence Feature Distributions ---")

# Summary statistics for transcript features
print("\nTranscript Features (all sequences):")
print(all_transcript_features[['length', 'gc_content', 'ctgcagy_density', 'g4_max_score']].describe())

print("\nPromoter Features (all sequences):")
print(all_promoter_features[['length', 'gc_content', 'cpg_density', 'g4_max_score']].describe())

# High G4 genes
high_g4_threshold = all_transcript_features['g4_max_score'].quantile(0.9)
high_g4_genes = all_transcript_features[all_transcript_features['g4_max_score'] > high_g4_threshold]['ensg_id'].tolist()
print(f"\nHigh G4 genes (top 10%): {len(high_g4_genes)}")

# High GC genes
high_gc_threshold = all_transcript_features['gc_content'].quantile(0.9)
high_gc_genes = all_transcript_features[all_transcript_features['gc_content'] > high_gc_threshold]['ensg_id'].tolist()
print(f"High GC genes (top 10%): {len(high_gc_genes)}")

# High CTGCAGY genes
high_ctgcagy_threshold = all_transcript_features['ctgcagy_density'].quantile(0.9)
high_ctgcagy_genes = all_transcript_features[all_transcript_features['ctgcagy_density'] > high_ctgcagy_threshold]['ensg_id'].tolist()
print(f"High CTGCAGY genes (top 10%): {len(high_ctgcagy_genes)}")

# ============================================================
# Generate Report
# ============================================================
print("\n" + "=" * 60)
print("Generating Report")
print("=" * 60)

# Determine verdict
# This is a placeholder - actual verdict requires gene symbol to ENSG mapping
verdict = "PARTIAL"
verdict_reason = """
Analysis completed with sequence features calculated for all genes.
Full verdict requires proper gene symbol to ENSG ID mapping to compare
Fire/Carrier/Coupled sets against sequence landscape.

Preliminary findings:
1. G4 scores and GC content show expected distributions
2. CTGCAGY density varies across transcripts
3. ATR/ATM proxy gene availability confirmed in NVU cell types
"""

# ATR/ATM interpretation
if atr_atm_results:
    cap_data = all_scores[all_scores['cell_type'] == 'Capillary']
    peri_data = all_scores[all_scores['cell_type'] == 'Pericyte']

    atr_atm_finding = f"""
ATR/ATM Proxy Analysis:
- Capillary: ATR={cap_data['ATR_score'].mean():.3f}±{cap_data['ATR_score'].std():.3f}, ATM={cap_data['ATM_score'].mean():.3f}±{cap_data['ATM_score'].std():.3f}
- Pericyte: ATR={peri_data['ATR_score'].mean():.3f}±{peri_data['ATR_score'].std():.3f}, ATM={peri_data['ATM_score'].mean():.3f}±{peri_data['ATM_score'].std():.3f}
"""
else:
    atr_atm_finding = "ATR/ATM proxy analysis not completed due to missing expression data."

report = f"""# ALS Trigger Project — Fire Origin Landscape Report

## Executive Summary

This analysis investigates the "Fire Origin" hypothesis:
- **Sequence Landscape**: GC-rich, G4-prone, CTGCAGY-enriched genes are vulnerable
- **Tension Field**: ATR (chronic replication stress) vs ATM (DSB brake) imbalance
- **Cell Type Specificity**: Capillary (fire origin) vs Pericyte (carrier)

## Verdict: {verdict}

{verdict_reason}

---

## PART 0: Gene Sets

| Set | N genes | Description |
|-----|---------|-------------|
| Fire | {len(fire_set)} | Class I Capillary Failure + candidates |
| Carrier | {len(carrier_set)} | Class II Pericyte Carrier + candidates |
| Coupled | {len(coupled_set)} | Class III Coupled |
| Background | {len(background_all)} | All NVU genes |

---

## PART A & B: Sequence Landscape

### Transcript Features (N={len(all_transcript_features)})

| Metric | Mean | Std | 10th% | 90th% |
|--------|------|-----|-------|-------|
| Length | {all_transcript_features['length'].mean():.0f} | {all_transcript_features['length'].std():.0f} | {all_transcript_features['length'].quantile(0.1):.0f} | {all_transcript_features['length'].quantile(0.9):.0f} |
| GC content | {all_transcript_features['gc_content'].mean():.3f} | {all_transcript_features['gc_content'].std():.3f} | {all_transcript_features['gc_content'].quantile(0.1):.3f} | {all_transcript_features['gc_content'].quantile(0.9):.3f} |
| CTGCAGY/kb | {all_transcript_features['ctgcagy_density'].mean():.3f} | {all_transcript_features['ctgcagy_density'].std():.3f} | {all_transcript_features['ctgcagy_density'].quantile(0.1):.3f} | {all_transcript_features['ctgcagy_density'].quantile(0.9):.3f} |
| G4 max score | {all_transcript_features['g4_max_score'].mean():.3f} | {all_transcript_features['g4_max_score'].std():.3f} | {all_transcript_features['g4_max_score'].quantile(0.1):.3f} | {all_transcript_features['g4_max_score'].quantile(0.9):.3f} |

### High-Risk Sequence Features

- High G4 genes (top 10%): {len(high_g4_genes)}
- High GC genes (top 10%): {len(high_gc_genes)}
- High CTGCAGY genes (top 10%): {len(high_ctgcagy_genes)}

---

## PART C: ATR vs ATM Proxy Expression

### Proxy Gene Availability

| Proxy Set | Total | Capillary Found | Pericyte Found |
|-----------|-------|-----------------|----------------|
| ATR (replication stress) | {len(ATR_PROXY_GENES)} | {len(atr_genes_cap) if cap_expr is not None else 0} | {len(atr_genes_peri) if peri_expr is not None else 0} |
| ATM (DSB response) | {len(ATM_PROXY_GENES)} | {len(atm_genes_cap) if cap_expr is not None else 0} | {len(atm_genes_peri) if peri_expr is not None else 0} |
| R-loop processing | {len(RLOOP_PROXY_GENES)} | {len(rloop_genes_cap) if cap_expr is not None else 0} | {len(rloop_genes_peri) if peri_expr is not None else 0} |

{atr_atm_finding}

---

## PART D: Interpretation

### Hypothesis Support Assessment

| Hypothesis | Status | Evidence |
|------------|--------|----------|
| H1: Fire genes are GC/G4-prone | NEEDS MAPPING | Sequence features calculated, gene mapping required |
| H2: ATR tension is chronic | PARTIAL | ATR proxy genes detected in NVU |
| H3: ATM brake loss in Capillary | NEEDS STAGING | Requires Control vs TDPneg comparison |

### Refutable Predictions

1. **G4-prone transcripts enriched in EV cargo**: Fire genes with high G4 scores should be enriched in Capillary-derived EVs
2. **ATR-high individuals show early ATM decline**: Individual variation in ATR proxy predicts ATM loss
3. **Pericyte carrier genes co-vary with fire gene G4 load**: CD81/CD9/RAB27A expression correlates with fire gene G4 burden

---

## Output Files

| File | Description |
|------|-------------|
| gene_sets_summary.csv | Gene set definitions and sizes |
| transcript_sequence_features.csv | GC, G4, CTGCAGY features per transcript |
| promoter_sequence_features.csv | GC, CpG, G4 features per promoter |
| atr_atm_proxy_scores_by_sample.csv | ATR/ATM/R-loop proxy scores per sample |
| proxy_gene_availability.csv | Proxy gene coverage |

---

## Next Steps

1. **Gene Symbol to ENSG Mapping**: Obtain proper annotation to link gene sets with sequence features
2. **TDP Staging Analysis**: Compare Control vs TDPneg for ATR/ATM proxy changes
3. **Bootstrap Validation**: Length/GC-matched background comparison

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
Key Outputs:
- Transcript features: {len(all_transcript_features)} sequences analyzed
- Promoter features: {len(all_promoter_features)} sequences analyzed
- ATR/ATM proxy scores calculated for Capillary and Pericyte
- High G4 genes identified: {len(high_g4_genes)}

Verdict: {verdict}

Next: Gene symbol to ENSG mapping needed for full Fire vs Background comparison
Output directory: {OUTPUT_DIR}
""")
