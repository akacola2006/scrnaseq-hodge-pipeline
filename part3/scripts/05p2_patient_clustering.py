#!/usr/bin/env python3
"""
Phase 5′-2: Patient Clustering (Subtype Extraction)
====================================================

Identifies ALS subtypes based on patient-level progression patterns.

Uses patient × cell type PT_IDS patterns to cluster patients into subtypes
(e.g., Oligo-driven, Upper-layer-driven, Mixed).

Output:
- patient_clusters.csv: Patient cluster assignments
- patient_clusters_PCA.png: 2D visualization
- PATIENT_STRATIFIED_OVERVIEW.md: Cluster characteristics
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

# Paths
base_dir = Path('/home/akaco/als/motor_cortex_analysis')
ids_dir = base_dir / 'ids_causal_analysis'
results_dir = ids_dir / 'results' / 'patient_stratified'

print("="*80)
print("Phase 5′-2: Patient Clustering (Subtype Extraction)")
print("="*80)

# ============================================================================
# 1. Load Patient × Cell Type Summary
# ============================================================================
print("\n[1] Loading patient × cell type summary...")

summary_df = pd.read_csv(results_dir / 'patient_celltype_summary.csv')
progression_matrix = pd.read_csv(results_dir / 'patient_progression_matrix.csv', index_col=0)

print(f"  Summary: {len(summary_df)} patient × cell type combinations")
print(f"  Progression matrix: {progression_matrix.shape}")

# Focus on ALS patients only for subtype analysis
als_patients = progression_matrix[progression_matrix['condition'] == 'ALS'].index.tolist()
print(f"\n  ALS patients for clustering: {len(als_patients)}")

# ============================================================================
# 2. Create Patient Feature Matrix
# ============================================================================
print("\n[2] Creating patient feature matrix...")

# Select key cell types for clustering
key_cell_types = [
    'Ex.L5.VAT1L.EYA4',
    'Ex.L5.VAT1L.THSD4',
    'Glia.Oligo',
    'Glia.Astro.GFAP-neg',
    'Glia.Micro',
    'Ex.L2.L3.CUX2.RASGRF2',
    'Ex.L3.L5.CUX2.RORB',
    'Ex.L4.L6.RORB.LRRK1',
    'Ex.L5.L6.THEMIS.TMEM233',
    'In.SOM.SST.GALNT14'
]

# Extract feature matrix for ALS patients
feature_cols = [col for col in key_cell_types if col in progression_matrix.columns]
print(f"  Using {len(feature_cols)} cell types as features")

X = progression_matrix.loc[als_patients, feature_cols].copy()

# Handle missing values (some patients might not have all cell types)
X_filled = X.fillna(X.mean())  # Fill missing with column mean
print(f"  Feature matrix: {X_filled.shape}")
print(f"  Missing values filled: {X.isna().sum().sum()}")

# ============================================================================
# 3. Standardize and PCA
# ============================================================================
print("\n[3] Standardization and PCA...")

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_filled)

# PCA
pca = PCA(n_components=min(10, X_scaled.shape[1]))
X_pca = pca.fit_transform(X_scaled)

print(f"  PCA components: {X_pca.shape[1]}")
print(f"  Variance explained (first 3 PCs): {pca.explained_variance_ratio_[:3].sum():.2%}")

# ============================================================================
# 4. K-Means Clustering
# ============================================================================
print("\n[4] K-Means clustering...")

# Try different k values
k_values = [3, 4, 5, 6]
inertias = []
silhouette_scores = []

for k in k_values:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=50)
    labels = kmeans.fit_predict(X_pca[:, :5])  # Use first 5 PCs
    inertias.append(kmeans.inertia_)

    # Calculate silhouette score
    from sklearn.metrics import silhouette_score
    score = silhouette_score(X_pca[:, :5], labels)
    silhouette_scores.append(score)
    print(f"  k={k}: inertia={kmeans.inertia_:.2f}, silhouette={score:.3f}")

# Select k based on silhouette score
best_k = k_values[np.argmax(silhouette_scores)]
print(f"\n  Selected k={best_k} (highest silhouette score)")

# Final clustering with best k
kmeans_final = KMeans(n_clusters=best_k, random_state=42, n_init=50)
cluster_labels = kmeans_final.fit_predict(X_pca[:, :5])

# ============================================================================
# 5. Create Cluster Assignment DataFrame
# ============================================================================
print("\n[5] Creating cluster assignments...")

cluster_df = pd.DataFrame({
    'patient_id': als_patients,
    'cluster_id': cluster_labels,
    'condition': 'ALS'
})

# Add global PT and stress from progression matrix
cluster_df = cluster_df.merge(
    progression_matrix[['PT_IDS_global_median', 'stress_global_mean', 'total_cells']],
    left_on='patient_id',
    right_index=True
)

# Add PCA coordinates
cluster_df['PC1'] = X_pca[:, 0]
cluster_df['PC2'] = X_pca[:, 1]

# Save cluster assignments
output_file = results_dir / 'patient_clusters.csv'
cluster_df.to_csv(output_file, index=False)
print(f"\n  Saved: {output_file}")

# Print cluster summary
print(f"\n  Cluster sizes:")
for cluster_id in range(best_k):
    n = (cluster_labels == cluster_id).sum()
    print(f"    Cluster {cluster_id}: {n} patients")

# ============================================================================
# 6. Visualize Clusters
# ============================================================================
print("\n[6] Creating visualizations...")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Plot 1: PCA scatter
ax = axes[0]
for cluster_id in range(best_k):
    mask = cluster_labels == cluster_id
    ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
               label=f'Cluster {cluster_id}',
               s=100, alpha=0.7)

# Add patient labels
for i, patient in enumerate(als_patients):
    ax.annotate(patient, (X_pca[i, 0], X_pca[i, 1]),
                fontsize=8, alpha=0.6)

ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} var)')
ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} var)')
ax.set_title('Patient Clustering (PCA)')
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 2: Global PT by cluster
ax = axes[1]
cluster_pt = [cluster_df[cluster_df['cluster_id'] == i]['PT_IDS_global_median'].values
              for i in range(best_k)]
bp = ax.boxplot(cluster_pt, labels=[f'C{i}' for i in range(best_k)])
ax.set_xlabel('Cluster')
ax.set_ylabel('Global PT_IDS Median')
ax.set_title('Disease Progression by Cluster')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
output_fig = results_dir / 'patient_clusters_PCA.png'
plt.savefig(output_fig, dpi=300, bbox_inches='tight')
plt.close()
print(f"  Saved: {output_fig}")

# ============================================================================
# 7. Analyze Cluster Characteristics
# ============================================================================
print("\n[7] Analyzing cluster characteristics...")

# For each cluster, calculate mean PT_IDS for each cell type
cluster_characteristics = []

for cluster_id in range(best_k):
    cluster_patients = cluster_df[cluster_df['cluster_id'] == cluster_id]['patient_id'].tolist()

    char = {'cluster_id': cluster_id, 'n_patients': len(cluster_patients)}

    # Mean PT for each key cell type
    for celltype in feature_cols:
        pt_values = progression_matrix.loc[cluster_patients, celltype].dropna()
        if len(pt_values) > 0:
            char[f'{celltype}_PT_mean'] = pt_values.mean()
        else:
            char[f'{celltype}_PT_mean'] = np.nan

    # Global metrics
    char['global_PT_mean'] = cluster_df[cluster_df['cluster_id'] == cluster_id]['PT_IDS_global_median'].mean()
    char['global_stress_mean'] = cluster_df[cluster_df['cluster_id'] == cluster_id]['stress_global_mean'].mean()

    cluster_characteristics.append(char)

char_df = pd.DataFrame(cluster_characteristics)

# Identify cluster signatures (which cell types are earliest/latest)
print(f"\n  Cluster signatures:")
for cluster_id in range(best_k):
    row = char_df[char_df['cluster_id'] == cluster_id].iloc[0]

    # Get PT values for cell types
    pt_cols = [col for col in char_df.columns if col.endswith('_PT_mean')]
    pt_values = row[pt_cols].dropna()

    if len(pt_values) > 0:
        # Sort by PT (earliest = lowest PT)
        pt_sorted = pt_values.sort_values()

        earliest_3 = pt_sorted.head(3)
        latest_3 = pt_sorted.tail(3)

        print(f"\n    Cluster {cluster_id} ({int(row['n_patients'])} patients):")
        print(f"      Global PT: {row['global_PT_mean']:.3f}")
        print(f"      Earliest cell types:")
        for celltype_full, pt in earliest_3.items():
            celltype = celltype_full.replace('_PT_mean', '').split('.')[-1]
            print(f"        {celltype}: PT={pt:.3f}")
        print(f"      Latest cell types:")
        for celltype_full, pt in latest_3.items():
            celltype = celltype_full.replace('_PT_mean', '').split('.')[-1]
            print(f"        {celltype}: PT={pt:.3f}")

# ============================================================================
# 8. Create Overview Markdown
# ============================================================================
print("\n[8] Creating overview report...")

with open(results_dir / 'PATIENT_STRATIFIED_OVERVIEW.md', 'w') as f:
    f.write("# Patient Stratification Overview\n")
    f.write("## Phase 5′-2: ALS Subtype Identification\n\n")
    f.write(f"**Date**: 2025-11-23\n")
    f.write(f"**Method**: K-Means clustering (k={best_k}) on patient × cell type PT patterns\n\n")
    f.write("---\n\n")

    f.write("## Summary\n\n")
    f.write(f"- **Total ALS patients**: {len(als_patients)}\n")
    f.write(f"- **Identified clusters**: {best_k}\n")
    f.write(f"- **Features used**: {len(feature_cols)} key cell types\n")
    f.write(f"- **PCA variance explained (PC1+PC2)**: {pca.explained_variance_ratio_[:2].sum():.1%}\n\n")

    f.write("---\n\n")
    f.write("## Cluster Characteristics\n\n")

    for cluster_id in range(best_k):
        cluster_patients = cluster_df[cluster_df['cluster_id'] == cluster_id]['patient_id'].tolist()
        row = char_df[char_df['cluster_id'] == cluster_id].iloc[0]

        f.write(f"### Cluster {cluster_id}\n\n")
        f.write(f"- **Patients**: {len(cluster_patients)}\n")
        f.write(f"  - IDs: {', '.join(cluster_patients)}\n")
        f.write(f"- **Global PT mean**: {row['global_PT_mean']:.3f}\n")
        f.write(f"- **Global stress mean**: {row['global_stress_mean']:.3f}\n\n")

        # Get PT values
        pt_cols = [col for col in char_df.columns if col.endswith('_PT_mean')]
        pt_values = row[pt_cols].dropna()

        if len(pt_values) > 0:
            pt_sorted = pt_values.sort_values()

            f.write("**Earliest cell types** (lowest PT):\n")
            for celltype_full, pt in pt_sorted.head(3).items():
                celltype = celltype_full.replace('_PT_mean', '')
                f.write(f"- {celltype}: PT={pt:.3f}\n")

            f.write("\n**Latest cell types** (highest PT):\n")
            for celltype_full, pt in pt_sorted.tail(3).items():
                celltype = celltype_full.replace('_PT_mean', '')
                f.write(f"- {celltype}: PT={pt:.3f}\n")

        f.write("\n")

    f.write("---\n\n")
    f.write("## Interpretation\n\n")
    f.write("*To be completed after representative patient analysis (Phase 5′-3)*\n\n")

    f.write("---\n\n")
    f.write(f"**Analysis completed**: Phase 5′-2\n")
    f.write(f"**Next step**: Phase 5′-3 (Representative patient state DAG analysis)\n")

print(f"  Saved: {results_dir / 'PATIENT_STRATIFIED_OVERVIEW.md'}")

print("\n" + "="*80)
print(f"Phase 5′-2 completed!")
print(f"Identified {best_k} patient clusters")
print(f"Output files in: {results_dir}")
print("="*80)
