"""
IDH-stratified condition φ: IDH-mut (W2, G3) vs IDH-wt (W3).
All-gene insertion with GPU acceleration.
Isolates the direct co-expression effect of IDH mutation.
"""
import sys, os
os.environ["PYTHONIOENCODING"] = "utf-8"
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

import numpy as np
import pandas as pd
import torch
import json
from scipy.stats import spearmanr, ttest_1samp
from pathlib import Path
import time

BASE_DIR = Path("D:/Projects/glioma_hodge")
OUTPUT = BASE_DIR / "condition_phi" / "idh_stratified"
OUTPUT.mkdir(parents=True, exist_ok=True)

N_BASE = 2000
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {device}")
if device.type == "cuda":
    print(f"GPU: {torch.cuda.get_device_name(0)}")


def gpu_log_corr(X_torch):
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
print("Loading data...")
expr_base = pd.read_parquet(BASE_DIR / "data/expression_gbmlgg_gtex.parquet")
expr_full = pd.read_parquet(BASE_DIR / "data/expression_full_gbmlgg_gtex.parquet")
meta = pd.read_csv(BASE_DIR / "data/metadata_gbmlgg_gtex.csv")

with open(BASE_DIR / "data/gene_names.json") as f:
    base_genes = json.load(f)["hvg_genes"]

base_set = set(base_genes)
extra_genes_list = [g for g in expr_full.index if g not in base_set]
all_genes = base_genes + extra_genes_list
common_samples = [s for s in expr_base.columns if s in expr_full.columns]
expr = pd.concat([expr_base[common_samples], expr_full.loc[extra_genes_list, common_samples]])
print(f"Combined: {expr.shape[0]} genes x {expr.shape[1]} samples")

# IDH-mut (W2: LGG G3 IDH-mut) vs IDH-wt (W3: LGG IDH-wt)
mut_ids = [s for s in meta[(meta["window"] == 2) & (meta["idh_status"] == "Mutant")]["sample_id"]
           if s in expr.columns]
wt_ids = [s for s in meta[(meta["window"] == 3) & (meta["idh_status"] == "WT")]["sample_id"]
          if s in expr.columns]
print(f"IDH-mut (W2 G3): {len(mut_ids)} samples")
print(f"IDH-wt  (W3):    {len(wt_ids)} samples")

# Gene indices
gene_to_col = {g: i for i, g in enumerate(all_genes)}
base_col_idx = np.array([gene_to_col[g] for g in base_genes])
sample_to_col = {s: i for i, s in enumerate(expr.columns)}

# Upload to GPU
print("Uploading to GPU...")
X_all = expr.values.astype(np.float32)

mut_col_indices = [sample_to_col[s] for s in mut_ids]
wt_col_indices = [sample_to_col[s] for s in wt_ids]

X_mut = torch.from_numpy(X_all[:, mut_col_indices].T.copy()).to(device)  # samples x genes
X_wt = torch.from_numpy(X_all[:, wt_col_indices].T.copy()).to(device)

X_mut_base = X_mut[:, base_col_idx]
X_wt_base = X_wt[:, base_col_idx]

print(f"  IDH-mut base: {X_mut_base.shape}, IDH-wt base: {X_wt_base.shape}")
if device.type == "cuda":
    print(f"  GPU memory: {torch.cuda.memory_allocated()/1e9:.2f} GB")

# ── Base phi ──
print("Computing base log-corr...")
log_mut = gpu_log_corr(X_mut_base)
log_wt = gpu_log_corr(X_wt_base)

# Condition: IDH-mut - IDH-wt (effect of having the mutation)
delta_base = log_mut - log_wt
# Static: IDH-wt only (the "normal" reference)
phi_base_cond, _ = compute_phi_full(delta_base, N_BASE)
phi_base_static, _ = compute_phi_full(log_wt, N_BASE)

rho, _ = spearmanr(phi_base_cond, phi_base_static)
print(f"Base R2(cond~static): {rho**2:.3f}")

# Base gene percentiles
all_phi_cond_pct = {}
all_phi_static_pct = {}
for i, g in enumerate(base_genes):
    rank_c = int((phi_base_cond > phi_base_cond[i]).sum()) + 1
    rank_s = int((phi_base_static > phi_base_static[i]).sum()) + 1
    all_phi_cond_pct[g] = rank_c / N_BASE * 100
    all_phi_static_pct[g] = rank_s / N_BASE * 100

# ── Resume ──
cp_path = OUTPUT / "checkpoint.npz"
if cp_path.exists():
    cp = np.load(cp_path, allow_pickle=True)
    all_phi_cond_pct.update(cp["cond_pct"].item())
    all_phi_static_pct.update(cp["static_pct"].item())
    print(f"Resumed: {len(all_phi_cond_pct)} genes done")

extra_genes = [g for g in extra_genes_list if g not in all_phi_cond_pct]
print(f"\nInserting {len(extra_genes)} remaining genes...")

# ── Insertion loop ──
t0 = time.time()
n_done = 0
n_failed = 0
n_ins = N_BASE + 1

for gi, gene in enumerate(extra_genes):
    col_idx = gene_to_col[gene]

    X_mut_ins = torch.cat([X_mut_base, X_mut[:, col_idx:col_idx+1]], dim=1)
    X_wt_ins = torch.cat([X_wt_base, X_wt[:, col_idx:col_idx+1]], dim=1)

    try:
        log_mut_ins = gpu_log_corr(X_mut_ins)
        log_wt_ins = gpu_log_corr(X_wt_ins)
    except Exception:
        n_failed += 1
        continue

    delta_ins = log_mut_ins - log_wt_ins
    rank_c = fast_phi_rank(delta_ins, N_BASE)
    all_phi_cond_pct[gene] = rank_c / n_ins * 100

    rank_s = fast_phi_rank(log_wt_ins, N_BASE)
    all_phi_static_pct[gene] = rank_s / n_ins * 100

    n_done += 1
    total = n_done + n_failed
    if total % 100 == 0 or n_done == 1:
        elapsed = time.time() - t0
        rate = total / elapsed if elapsed > 0 else 0
        eta = (len(extra_genes) - total) / rate / 60 if rate > 0 else 0
        print(f"  {total}/{len(extra_genes)} ({n_done}ok {n_failed}fail), {rate:.1f}/s, ETA {eta:.0f}min")

    if n_done % 1000 == 0:
        np.savez(cp_path, cond_pct=all_phi_cond_pct, static_pct=all_phi_static_pct)

elapsed_total = time.time() - t0
print(f"\nDone: {n_done}/{len(extra_genes)} in {elapsed_total/60:.1f}min ({n_done/elapsed_total:.1f}/s)")

# ── Build result ──
rows = []
for g in all_genes:
    if g not in all_phi_cond_pct:
        continue
    rows.append({
        "gene": g,
        "phi_cond_pct": all_phi_cond_pct[g],
        "phi_static_pct": all_phi_static_pct.get(g, np.nan),
        "in_base": g in base_set,
    })
result = pd.DataFrame(rows)

sorted_phi_cond = np.sort(phi_base_cond)[::-1]
sorted_phi_static = np.sort(phi_base_static)[::-1]

def inv_cdf(pct, arr):
    n = len(arr)
    pos = np.clip(pct / 100.0 * n, 0, n - 1)
    lo = int(np.floor(pos)); hi = min(lo + 1, n - 1); frac = pos - lo
    return arr[lo] * (1 - frac) + arr[hi] * frac

result["raw_phi_cond"] = result["phi_cond_pct"].apply(lambda p: inv_cdf(p, sorted_phi_cond))
result["raw_phi_static"] = result["phi_static_pct"].apply(lambda p: inv_cdf(p, sorted_phi_static))

valid = result.dropna(subset=["raw_phi_cond", "raw_phi_static"])
coeffs = np.polyfit(valid["raw_phi_static"], valid["raw_phi_cond"], 3)
result["resid"] = result["raw_phi_cond"] - np.polyval(coeffs, result["raw_phi_static"])
result["z_resid"] = (result["resid"] - result["resid"].mean()) / result["resid"].std()

result.to_csv(OUTPUT / "glioma_idh_stratified_allgene.csv", index=False)
print(f"Saved: {len(result)} genes")

# ── Summary ──
rho_final, _ = spearmanr(result["raw_phi_cond"].dropna(), result["raw_phi_static"].dropna())
n_rew = (result["z_resid"] > 1).sum()
n_col = (result["z_resid"] < -1).sum()
print(f"\nR2(cond~static): {rho_final**2:.3f}")
print(f"REWIRING (z>+1): {n_rew}, COLLAPSE (z<-1): {n_col}")

# Key genes
print("\nKey genes (IDH-mut effect):")
for g in ["IDH1","IDH2","TP53","EGFR","PTEN","CDKN2A","ATRX","CIC","FUBP1",
          "NF1","PIK3CA","RB1","MDM2","CDK4","PDGFRA","MET","TERT",
          "TARDBP","NEMF","FUS","SOD1",
          "RPL6","RPL3","RPS7","EIF4E",
          "MOG","MBP","NEFL","VEGFA","TOP2A","MKI67",
          "TET2","DNMT3A","KDM5A","SETD2","EZH2",
          "MGMT","HIF1A","MYC","NOTCH1"]:
    m = result[result["gene"] == g]
    if len(m) == 0:
        continue
    r = m.iloc[0]
    b = "Y" if r["in_base"] else "N"
    flag = ""
    if r.z_resid > 1: flag = " ** REWIRING"
    elif r.z_resid < -1: flag = " ** COLLAPSE"
    print(f"  {g:>10s}: cond={r.phi_cond_pct:5.1f}% stat={r.phi_static_pct:5.1f}% z={r.z_resid:+.2f} [{b}]{flag}")

# Translation
print("\nTranslation class tests:")
for cls, prefixes in [("RPL/RPS", ["RPL", "RPS"]), ("EIF/EEF", ["EIF", "EEF"]),
                       ("Translation full", ["RPL", "RPS", "EIF", "EEF"])]:
    m = result[result["gene"].apply(lambda g: any(str(g).startswith(p) for p in prefixes))]
    if len(m) < 5:
        continue
    t, p = ttest_1samp(m["z_resid"].dropna(), 0)
    d = "REWIRING" if m["z_resid"].mean() > 0 else "COLLAPSE"
    print(f"  {cls:20s}: n={len(m):4d}, mean_z={m['z_resid'].mean():+.3f}, p={p:.3e} [{d}]")

# Enrichment
print("\nEnrichment...")
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

print("\nDone.")
