"""
Step 9: Translation gene insertion experiment for glioma
========================================================
Inserts translation-related genes (RPL, RPS, EIF, EEF) one at a time
into the 2000-gene set, runs Hodge decomposition, and records where
each inserted gene lands in the phi ranking.

This directly addresses the criticism:
  "Ribosomal genes aren't in the gene set, so of course translation
   doesn't come up as upstream."

If inserted ribosomal genes land at the phi-top (upstream), the sALS
translation-upstream finding could be an artifact. If they do NOT land
upstream, it confirms disease-specificity.

Usage:
    python step9_gene_insertion.py
"""
import gzip
import json
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.covariance import LedoitWolf

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

BASE = Path(r"D:\Projects\IDS paper writing\glioma_hodge")
DATA_DIR = BASE / "data"
OUT_DIR = BASE / "results" / "gene_insertion"


# ── SPD / Hodge utilities (from step2) ────────────────────────

def spd_log(A):
    A_sym = (A + A.T) / 2.0
    eigvals, eigvecs = np.linalg.eigh(A_sym)
    eigvals = np.maximum(eigvals, 1e-15)
    log_A = eigvecs @ np.diag(np.log(eigvals)) @ eigvecs.T
    return (log_A + log_A.T) / 2.0


def compute_correlation_lw(X):
    lw = LedoitWolf()
    lw.fit(X)
    cov = lw.covariance_
    std = np.sqrt(np.diag(cov))
    std = np.maximum(std, 1e-15)
    corr = cov / np.outer(std, std)
    np.fill_diagonal(corr, 1.0)
    eigvals, eigvecs = np.linalg.eigh(corr)
    eigvals = np.maximum(eigvals, 1e-10)
    corr = (eigvecs * eigvals) @ eigvecs.T
    return (corr + corr.T) / 2.0


def hodge_phi_kn(flow, n_nodes):
    edges_i, edges_j = np.triu_indices(n_nodes, k=1)
    edges = np.stack([edges_i, edges_j], axis=1).astype(np.int32)
    div = np.zeros(n_nodes, dtype=np.float64)
    np.add.at(div, edges[:, 0], -flow)
    np.add.at(div, edges[:, 1], +flow)
    phi = (div - div.mean()) / n_nodes
    gradient = phi[edges[:, 1]] - phi[edges[:, 0]]
    flow_sq = np.dot(flow, flow)
    grad_sq = np.dot(gradient, gradient)
    gf = float(grad_sq / flow_sq) if flow_sq > 1e-30 else 0.0
    return phi, gf


def build_flow_edge_weight(delta):
    n_genes = delta.shape[0]
    gene_d = np.linalg.norm(delta, axis=1)
    n_edges = n_genes * (n_genes - 1) // 2
    flow = np.empty(n_edges, dtype=np.float64)
    idx = 0
    for i in range(n_genes):
        js = np.arange(i + 1, n_genes)
        delta_ij = np.abs(delta[i, js])
        diff = gene_d[i] - gene_d[js]
        flow[idx:idx + len(js)] = delta_ij * np.sign(diff)
        idx += len(js)
    flow_std = float(np.std(flow))
    if flow_std > 1e-15:
        flow = flow / flow_std
    return flow


# ── Data loading ──────────────────────────────────────────────

def load_full_expression():
    """Load full TOIL expression matrix and map ENSG IDs to gene symbols."""
    logger.info("Loading probemap for ENSG → gene symbol mapping...")
    probe = pd.read_csv(DATA_DIR / "raw" / "gencode.v23.annotation.gene.probemap", sep='\t')
    ensg_to_gene = dict(zip(probe['id'], probe['gene']))

    logger.info("Loading full expression matrix (this may take a while)...")
    # Read the gzipped TSV
    df = pd.read_csv(
        DATA_DIR / "raw" / "TcgaTargetGtex_gene_expected_count.gz",
        sep='\t', index_col=0, compression='gzip'
    )
    logger.info("  Raw shape: %s", df.shape)

    # Map ENSG IDs to gene symbols
    df.index = df.index.map(lambda x: ensg_to_gene.get(x, x))

    # Remove duplicates (keep first)
    df = df[~df.index.duplicated(keep='first')]
    logger.info("  After dedup: %d genes × %d samples", df.shape[0], df.shape[1])

    return df


def main():
    t0 = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("Step 9: Translation gene insertion experiment")
    logger.info("=" * 60)

    # Load existing 2000-gene expression data and metadata
    logger.info("Loading 2000-gene expression data...")
    expr_2k = pd.read_parquet(DATA_DIR / "expression_gbmlgg_gtex.parquet")
    meta = pd.read_csv(DATA_DIR / "metadata_gbmlgg_gtex.csv", index_col=0)
    gene_names_2k = list(expr_2k.index)
    sample_ids = list(expr_2k.columns)
    logger.info("  2000 genes, %d samples", len(sample_ids))

    # Align samples with metadata
    common_meta = [s for s in sample_ids if s in meta.index]
    mask_known = meta.loc[common_meta, "window_label"] != "LGG IDH-mut (grade unknown)"
    sample_ids_clean = [s for s, m in zip(common_meta, mask_known) if m]
    logger.info("  After excluding grade-unknown: %d samples", len(sample_ids_clean))

    # Load full expression data
    full_expr = load_full_expression()

    # Find translation genes to insert
    import re
    def is_pseudogene(name):
        return bool(re.search(r'P\d+$', name))

    all_full_genes = set(full_expr.index)
    existing_genes = set(gene_names_2k)

    # Core ribosomal proteins + translation factors
    target_prefixes = {
        'RPL': 'ribosomal_large',
        'RPS': 'ribosomal_small',
        'EIF': 'initiation_factor',
        'EEF': 'elongation_factor',
    }

    candidates = []
    for prefix, category in target_prefixes.items():
        genes = [g for g in all_full_genes
                 if g.startswith(prefix)
                 and not is_pseudogene(g)
                 and g not in existing_genes
                 and not g.endswith(('-AS1', '-AS2'))  # antisense
                 and 'MIR' not in g
                 and '-' not in g  # fusion genes
                 ]
        for g in sorted(genes):
            candidates.append((g, category))

    logger.info("  Candidate translation genes for insertion: %d", len(candidates))
    for cat in set(c[1] for c in candidates):
        genes_in_cat = [c[0] for c in candidates if c[1] == cat]
        logger.info("    %s: %d — %s", cat, len(genes_in_cat), genes_in_cat[:10])

    # Check which candidates have expression data for our samples
    valid_candidates = []
    for gene, cat in candidates:
        if gene in full_expr.index:
            # Check overlap with our sample IDs
            available = [s for s in sample_ids_clean if s in full_expr.columns]
            if len(available) > 600:  # need most samples
                valid_candidates.append((gene, cat))

    logger.info("  Valid candidates (with expression data): %d", len(valid_candidates))

    # Find common samples between full_expr and our sample set
    common_samples = [s for s in sample_ids_clean if s in full_expr.columns]
    logger.info("  Common samples: %d / %d", len(common_samples), len(sample_ids_clean))

    if len(common_samples) < 600:
        logger.error("Too few common samples. Aborting.")
        return

    # Prepare base expression matrix (2000 genes × common samples)
    X_base = expr_2k[common_samples].values.astype(np.float64)  # (2000, n_samples)
    meta_common = meta.loc[common_samples]

    # Window definitions (clinical grade)
    windows = sorted(meta_common['window'].unique())
    logger.info("  Windows: %s", windows)

    # Run insertion experiment for each candidate gene
    results = []
    n_genes_base = len(gene_names_2k)

    for gi, (gene, category) in enumerate(valid_candidates):
        if gi % 10 == 0:
            logger.info("  Processing %d/%d: %s (%s)", gi+1, len(valid_candidates), gene, category)

        # Get expression vector for this gene
        gene_expr = full_expr.loc[gene, common_samples].values.astype(np.float64)

        # Stack: 2001 genes (original 2000 + inserted gene at the end)
        X_aug = np.vstack([X_base, gene_expr.reshape(1, -1)])  # (2001, n_samples)
        n_genes_aug = n_genes_base + 1
        inserted_idx = n_genes_aug - 1

        # Run Hodge on each transition with augmented gene set
        phi_percentiles = []
        phi_values = []

        for t in range(len(windows) - 1):
            w_from, w_to = windows[t], windows[t + 1]

            mask_from = (meta_common['window'] == w_from).values
            mask_to = (meta_common['window'] == w_to).values

            X_from = X_aug[:, mask_from].T  # (n_samples_from, 2001)
            X_to = X_aug[:, mask_to].T

            if X_from.shape[0] < 30 or X_to.shape[0] < 30:
                continue

            # Correlation matrices
            corr_from = compute_correlation_lw(X_from)
            corr_to = compute_correlation_lw(X_to)

            # SPD log + delta
            log_from = spd_log(corr_from)
            log_to = spd_log(corr_to)
            delta = log_to - log_from

            # Flow + Hodge
            flow = build_flow_edge_weight(delta)
            phi, gf = hodge_phi_kn(flow, n_genes_aug)

            # Percentile of inserted gene
            ranks = stats.rankdata(phi)
            pctile = ranks[inserted_idx] / n_genes_aug * 100
            phi_percentiles.append(pctile)
            phi_values.append(float(phi[inserted_idx]))

        if phi_percentiles:
            mean_pctile = np.mean(phi_percentiles)
            mean_phi = np.mean(phi_values)
            results.append({
                'gene': gene,
                'category': category,
                'mean_percentile': float(mean_pctile),
                'mean_phi': float(mean_phi),
                'per_transition_percentile': [float(p) for p in phi_percentiles],
                'per_transition_phi': [float(p) for p in phi_values],
                'n_transitions': len(phi_percentiles),
            })

    # Sort by mean percentile
    results.sort(key=lambda x: x['mean_percentile'], reverse=True)

    # Summary
    logger.info("")
    logger.info("=" * 60)
    logger.info("RESULTS: Translation gene insertion")
    logger.info("=" * 60)

    for r in results:
        logger.info("  %s (%s): mean pctile=%.1f, mean phi=%+.6f",
                     r['gene'], r['category'], r['mean_percentile'], r['mean_phi'])

    # Category summary
    logger.info("")
    logger.info("Category summary:")
    for cat in sorted(set(r['category'] for r in results)):
        cat_results = [r for r in results if r['category'] == cat]
        pctiles = [r['mean_percentile'] for r in cat_results]
        logger.info("  %s (n=%d): mean=%.1f, median=%.1f, range=[%.1f, %.1f]",
                     cat, len(pctiles), np.mean(pctiles), np.median(pctiles),
                     min(pctiles), max(pctiles))

    # Overall: are translation genes upstream?
    all_pctiles = [r['mean_percentile'] for r in results]
    logger.info("")
    logger.info("Overall (n=%d): mean=%.1f, median=%.1f",
                 len(all_pctiles), np.mean(all_pctiles), np.median(all_pctiles))
    stat_w, p_w = stats.wilcoxon([p - 50 for p in all_pctiles])
    logger.info("Wilcoxon signed-rank (vs 50th percentile): p=%.4e", p_w)

    # Compare with sALS (where ribosomal genes are at ~13-26th percentile = upstream)
    n_upstream = sum(1 for p in all_pctiles if p > 75)
    n_downstream = sum(1 for p in all_pctiles if p < 25)
    logger.info("  Upstream (>75th pctile): %d/%d", n_upstream, len(all_pctiles))
    logger.info("  Downstream (<25th pctile): %d/%d", n_downstream, len(all_pctiles))

    # Save results
    output = {
        'gene_results': results,
        'summary': {
            'n_genes_tested': len(results),
            'mean_percentile': float(np.mean(all_pctiles)),
            'median_percentile': float(np.median(all_pctiles)),
            'wilcoxon_p': float(p_w),
            'n_upstream_75': n_upstream,
            'n_downstream_25': n_downstream,
        }
    }

    with open(OUT_DIR / "insertion_results.json", "w") as f:
        json.dump(output, f, indent=2)

    # Also save as CSV for easy viewing
    df_results = pd.DataFrame(results)
    df_results.to_csv(OUT_DIR / "insertion_results.csv", index=False)

    elapsed = time.time() - t0
    logger.info("")
    logger.info("Step 9 COMPLETE (%.1f s)", elapsed)
    logger.info("Results saved to %s", OUT_DIR)


if __name__ == "__main__":
    main()
