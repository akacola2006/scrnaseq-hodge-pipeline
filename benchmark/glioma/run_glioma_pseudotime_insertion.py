"""
All-gene insertion φ for glioma using pseudotime transitions (T0→T1, T1→T2, T2→T3, T3→T4).
GPU accelerated. Base=2000 genes, all extra genes inserted one-by-one.

For each inserted gene:
  - Compute LW log_corr for 5 windows (shared across transitions)
  - Derive 4 Δlog_corr matrices
  - Compute phi_rank for inserted gene in each transition
  - Record per-transition phi_cond_pct and phi_static_pct (static = window 0)
  - Average phi across transitions for final ranking
"""
import sys, os
os.environ["PYTHONIOENCODING"] = "utf-8"
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

import numpy as np
import pandas as pd
import torch
import json
from scipy.stats import spearmanr
from pathlib import Path
import time

BASE_DIR = Path("D:/Projects/glioma_hodge")
OUTPUT = BASE_DIR / "condition_phi" / "pseudotime_insertion"
OUTPUT.mkdir(parents=True, exist_ok=True)

N_BASE = 2000
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")
if device.type == "cuda":
    print(f"GPU: {torch.cuda.get_device_name(0)}")


def gpu_log_corr(X_torch):
    """Ledoit-Wolf shrinkage + log(corr) on GPU."""
    n, p = X_torch.shape
    Xc = X_torch - X_torch.mean(dim=0)
    S = (Xc.T @ Xc) / (n - 1)
    trS2 = (S * S).sum().item()
    trS = S.trace().item()
    denom = (n + 2) * (trS2 - trS ** 2 / p)
    if abs(denom) < 1e-15:
        alpha = 0.0
    else:
        alpha = min(1.0, max(0.0,
            ((n - 2) / n * trS2 + trS ** 2) / denom))
    S = (1 - alpha) * S + alpha * torch.eye(p, device=device, dtype=torch.float32) * (trS / p)
    std = torch.sqrt(torch.diag(S)).clamp(min=1e-15)
    corr = S / torch.outer(std, std)
    corr.fill_diagonal_(1.0)
    eigvals, eigvecs = torch.linalg.eigh(corr)
    eigvals = eigvals.clamp(min=1e-10)
    log_c = (eigvecs * torch.log(eigvals)) @ eigvecs.T
    return ((log_c + log_c.T) / 2).cpu().numpy().astype(np.float64)


def fast_phi_rank(M, inserted_idx):
    """Rank of inserted gene in phi (rank 1 = highest phi)."""
    d = np.linalg.norm(M, axis=1)
    abs_M = np.abs(M)
    sign_d = np.sign(d[np.newaxis, :] - d[:, np.newaxis])
    np.fill_diagonal(sign_d, 0)
    div = (abs_M * sign_d).sum(axis=1)
    div_ins = div[inserted_idx]
    return int((div > div_ins).sum()) + 1


def compute_phi_full(M, n):
    d = np.linalg.norm(M, axis=1)
    abs_M = np.abs(M)
    sd = np.sign(d[np.newaxis, :] - d[:, np.newaxis])
    np.fill_diagonal(sd, 0)
    div = (abs_M * sd).sum(axis=1)
    return (div - div.mean()) / n, d


# ── Load data ──
print("Loading glioma expression data...")
expr_base = pd.read_parquet(BASE_DIR / "data/expression_gbmlgg_gtex.parquet")  # 2000 base genes
expr_full = pd.read_parquet(BASE_DIR / "data/expression_full_gbmlgg_gtex.parquet")  # all genes
meta = pd.read_csv(BASE_DIR / "data/metadata_gbmlgg_gtex.csv")

with open(BASE_DIR / "data/gene_names.json") as f:
    gn = json.load(f)
base_genes = gn["hvg_genes"]

# Combine: base genes from filtered file, extra genes from full file
# Use filtered expr for base (guaranteed to have all 2000)
all_extra_genes = [g for g in expr_full.index if g not in set(base_genes)]
all_genes = base_genes + all_extra_genes
n_all = len(all_genes)
print(f"Total genes: {n_all} (Base: {len(base_genes)}, Extra: {len(all_extra_genes)})")

# Build combined expression matrix (genes × samples)
# Ensure same sample columns
common_samples = [s for s in expr_base.columns if s in expr_full.columns]
expr = pd.concat([expr_base[common_samples], expr_full.loc[all_extra_genes, common_samples]])
print(f"Combined expression: {expr.shape[0]} genes x {expr.shape[1]} samples")

# Window sample IDs
window_ids = {}
for w in range(5):
    ids = [s for s in meta[meta["window"] == w]["sample_id"] if s in expr.columns]
    window_ids[w] = ids
    print(f"  Window {w}: {len(ids)} samples")

# Gene indices
base_set = set(base_genes)
extra_genes = all_extra_genes
print(f"Extra genes for insertion: {len(extra_genes)}")

gene_to_col = {g: i for i, g in enumerate(all_genes)}
base_col_idx = np.array([gene_to_col[g] for g in base_genes])

# ── Upload to GPU ──
print("Uploading to GPU...")
# All expression: genes × samples → we need samples × genes per window
# Pre-extract per-window matrices
X_all = expr.values.astype(np.float32)  # genes × all_samples
sample_to_col = {s: i for i, s in enumerate(expr.columns)}

window_gpu_base = {}
window_gpu_all = {}
for w in range(5):
    col_indices = [sample_to_col[s] for s in window_ids[w]]
    X_w = X_all[:, col_indices].T  # samples × genes
    window_gpu_base[w] = torch.from_numpy(X_w[:, base_col_idx].copy()).to(device)
    window_gpu_all[w] = torch.from_numpy(X_w.copy()).to(device)
    print(f"  W{w}: {window_gpu_base[w].shape} (base), {window_gpu_all[w].shape} (all)")

if device.type == "cuda":
    print(f"GPU memory: {torch.cuda.memory_allocated()/1e9:.2f} GB")

# ── Base φ per transition ──
print("\nComputing base log-corr for all 5 windows...")
log_corr_base = {}
for w in range(5):
    log_corr_base[w] = gpu_log_corr(window_gpu_base[w])
    print(f"  W{w}: done")

transitions = [(0, 1), (1, 2), (2, 3), (3, 4)]
base_phi = {}
base_phi_static, _ = compute_phi_full(log_corr_base[0], N_BASE)

for w_from, w_to in transitions:
    delta = log_corr_base[w_to] - log_corr_base[w_from]
    phi, _ = compute_phi_full(delta, N_BASE)
    base_phi[(w_from, w_to)] = phi
    label = f"T{w_from}_{w_to}"
    rho, _ = spearmanr(phi, base_phi_static)
    print(f"  {label}: R2(phi~static) = {rho**2:.3f}")

# Base gene percentiles (per transition + static)
all_results = {}  # gene → dict of values
for i, g in enumerate(base_genes):
    d = {"gene": g, "in_base": True}
    for w_from, w_to in transitions:
        phi = base_phi[(w_from, w_to)]
        rank = int((phi > phi[i]).sum()) + 1
        d[f"phi_T{w_from}{w_to}_pct"] = rank / N_BASE * 100
    rank_s = int((base_phi_static > base_phi_static[i]).sum()) + 1
    d["phi_static_pct"] = rank_s / N_BASE * 100
    all_results[g] = d

# ── Resume from checkpoint ──
cp_path = OUTPUT / "checkpoint.npz"
if cp_path.exists():
    cp = np.load(cp_path, allow_pickle=True)
    saved = cp["results"].item()
    all_results.update(saved)
    print(f"Resumed from checkpoint: {len(all_results)} genes done")

extra_genes = [g for g in extra_genes if g not in all_results]
print(f"\nInserting {len(extra_genes)} remaining genes...")

# ── Main insertion loop ──
t0 = time.time()
n_done = 0
n_failed = 0
n_ins = N_BASE + 1

for gi, gene in enumerate(extra_genes):
    col_idx = gene_to_col[gene]

    try:
        # Compute LW log_corr for 5 windows (augmented with inserted gene)
        log_corr_ins = {}
        for w in range(5):
            X_ins = torch.cat([window_gpu_base[w], window_gpu_all[w][:, col_idx:col_idx+1]], dim=1)
            log_corr_ins[w] = gpu_log_corr(X_ins)
    except Exception:
        n_failed += 1
        continue

    d = {"gene": gene, "in_base": False}

    # phi per transition
    for w_from, w_to in transitions:
        delta_ins = log_corr_ins[w_to] - log_corr_ins[w_from]
        rank_c = fast_phi_rank(delta_ins, N_BASE)
        d[f"phi_T{w_from}{w_to}_pct"] = rank_c / n_ins * 100

    # phi_static (window 0 only)
    rank_s = fast_phi_rank(log_corr_ins[0], N_BASE)
    d["phi_static_pct"] = rank_s / n_ins * 100

    all_results[gene] = d
    n_done += 1

    total = n_done + n_failed
    if total % 100 == 0 or n_done == 1:
        elapsed = time.time() - t0
        rate = total / elapsed if elapsed > 0 else 0
        eta = (len(extra_genes) - total) / rate / 60 if rate > 0 else 0
        print(f"  {total}/{len(extra_genes)} ({n_done}ok {n_failed}fail), {rate:.1f}/s, ETA {eta:.0f}min")

    if n_done % 1000 == 0:
        np.savez(cp_path, results=all_results)

elapsed_total = time.time() - t0
print(f"\nDone: {n_done}/{len(extra_genes)} in {elapsed_total/60:.1f}min ({n_done/elapsed_total:.1f}/s)")

# ── Build result table ──
rows = []
for g in all_genes:
    if g not in all_results:
        continue
    d = all_results[g]
    row = {
        "gene": d["gene"],
        "in_base": d["in_base"],
        "phi_static_pct": d["phi_static_pct"],
    }
    # Per-transition phi
    phi_pcts = []
    for w_from, w_to in transitions:
        key = f"phi_T{w_from}{w_to}_pct"
        row[key] = d[key]
        phi_pcts.append(d[key])
    # Mean phi across transitions
    row["phi_cond_mean_pct"] = np.mean(phi_pcts)
    rows.append(row)

result = pd.DataFrame(rows)
print(f"Result: {len(result)} genes")

# Inverse CDF for mean phi_cond and phi_static
sorted_phi_static = np.sort(base_phi_static)[::-1]
# For mean phi_cond: use average of base phi across transitions
base_phi_mean = np.mean([base_phi[t] for t in transitions], axis=0)
sorted_phi_cond_mean = np.sort(base_phi_mean)[::-1]

def inv_cdf(pct, arr):
    n = len(arr)
    pos = np.clip(pct / 100.0 * n, 0, n - 1)
    lo = int(np.floor(pos)); hi = min(lo + 1, n - 1); frac = pos - lo
    return arr[lo] * (1 - frac) + arr[hi] * frac

result["raw_phi_cond"] = result["phi_cond_mean_pct"].apply(lambda p: inv_cdf(p, sorted_phi_cond_mean))
result["raw_phi_static"] = result["phi_static_pct"].apply(lambda p: inv_cdf(p, sorted_phi_static))

# Poly3 residual
valid = result.dropna(subset=["raw_phi_cond", "raw_phi_static"])
coeffs = np.polyfit(valid["raw_phi_static"], valid["raw_phi_cond"], 3)
result["resid"] = result["raw_phi_cond"] - np.polyval(coeffs, result["raw_phi_static"])
result["z_resid"] = (result["resid"] - result["resid"].mean()) / result["resid"].std()

result.to_csv(OUTPUT / "glioma_pt_allgene_insertion.csv", index=False)
print(f"Saved: {OUTPUT / 'glioma_pt_allgene_insertion.csv'}")

# ── Summary ──
rho, _ = spearmanr(result["raw_phi_cond"].dropna(), result["raw_phi_static"].dropna())
n_rew = (result["z_resid"] > 1).sum()
n_col = (result["z_resid"] < -1).sum()
print(f"\nR2(cond_mean~static): {rho**2:.3f}")
print(f"REWIRING (z>+1): {n_rew}, COLLAPSE (z<-1): {n_col}")

# Per-transition stats
for w_from, w_to in transitions:
    key = f"phi_T{w_from}{w_to}_pct"
    print(f"  {key}: mean={result[key].mean():.1f}%, median={result[key].median():.1f}%")

# Key genes
print("\nKey glioma genes:")
for g in ["IDH1","IDH2","TP53","EGFR","PTEN","CDKN2A","TERT","ATRX","CIC",
          "TARDBP","NEMF","RPL6","RPL3","RPS7","EIF4E","MOG","MBP","NEFL",
          "VEGFA","TOP2A","MKI67","CXCR4","PDGFRA","MET","CDK4"]:
    m = result[result["gene"] == g]
    if len(m) == 0:
        continue
    r = m.iloc[0]
    b = "Y" if r["in_base"] else "N"
    print(f"  {g:>10s}: cond_mean={r.phi_cond_mean_pct:5.1f}% stat={r.phi_static_pct:5.1f}% z={r.z_resid:+.2f} [{b}]")

# Translation class test
print("\nTranslation class tests:")
from scipy.stats import ttest_1samp
for cls, prefixes in [("RPL/RPS", ["RPL", "RPS"]), ("EIF/EEF", ["EIF", "EEF"]),
                       ("Translation full", ["RPL", "RPS", "EIF", "EEF"])]:
    m = result[result["gene"].apply(lambda g: any(str(g).startswith(p) for p in prefixes))]
    if len(m) < 5:
        continue
    t, p = ttest_1samp(m["z_resid"].dropna(), 0)
    direction = "REWIRING" if m["z_resid"].mean() > 0 else "COLLAPSE"
    print(f"  {cls:20s}: n={len(m):4d}, mean_z={m['z_resid'].mean():+.3f}, p={p:.3e} [{direction}]")

# Enrichment
print("\nRunning enrichment...")
import gseapy
LIBS = ["GO_Biological_Process_2025", "KEGG_2026", "Reactome_Pathways_2024"]
bg = [g for g in result["gene"].dropna().tolist() if isinstance(g, str) and len(g) > 1]

for label, mask in [("REWIRING", result["z_resid"] > 1), ("COLLAPSE", result["z_resid"] < -1)]:
    genes = [g for g in result[mask]["gene"].dropna().tolist() if isinstance(g, str) and len(g) > 1]
    if len(genes) < 10:
        continue
    try:
        enr = gseapy.enrich(gene_list=genes, gene_sets=LIBS, background=bg, outdir=None, no_plot=True, cutoff=0.05)
        res = enr.results
        if len(res) > 0:
            res = res[res["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value").drop_duplicates("Term", keep="first")
            print(f"\n{label} ({len(genes)} genes): {len(res)} sig terms")
            for _, r in res.head(15).iterrows():
                gs = r["Gene_set"].replace("GO_Biological_Process_2025","GO:BP").replace("KEGG_2026","KEGG").replace("Reactome_Pathways_2024","React")
                print(f"  {r['Term'][:60]:60s} p={r['Adjusted P-value']:.2e} [{gs}]")
        else:
            print(f"\n{label}: no significant terms")
    except Exception as e:
        print(f"\n{label}: error {e}")

# Save enrichment
res_all = []
for label, mask in [("REWIRING", result["z_resid"] > 1), ("COLLAPSE", result["z_resid"] < -1)]:
    genes = [g for g in result[mask]["gene"].dropna().tolist() if isinstance(g, str) and len(g) > 1]
    if len(genes) < 5:
        continue
    try:
        enr = gseapy.enrich(gene_list=genes, gene_sets=LIBS, background=bg, outdir=None, no_plot=True, cutoff=0.1)
        r = enr.results.copy()
        r["direction"] = label
        res_all.append(r)
    except Exception:
        pass
if res_all:
    pd.concat(res_all).to_csv(OUTPUT / "glioma_pt_enrichment.csv", index=False)

print("\nDone.")
