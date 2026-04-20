# ATM Fire Origin Evidence Summary (Scripts + Outputs)

This note summarizes what the identified scripts do, where they write outputs,
and the logic they use to support the "ATM brake loss" fire-origin hypothesis.
Content is descriptive only; no new analysis is run here.

Note: The scripts are written with a Linux-style BASE_DIR:
`/home/akaco/als/motor_cortex_analysis/ids_causal_analysis`.
In this workspace, analogous outputs are under `docs/ids_causal_analysis/...`.

## Scripts and Responsibilities

### 1) `docs/ids_causal_analysis/scripts/motif_family_nvu_rolesplit.py`
Purpose:
- Build k-mer motif families (CTGCAGY family).
- Rank "tool/cargo" candidates by motif features.
- Create donor-level pseudobulk (Control vs TDPneg) for NVU cell types.
- Classify genes into fire-origin (Capillary failure) vs carrier (Pericyte) roles.

Key inputs:
- `results_motif_screening/tables/transcript_kmer_enrichment_k7.csv`
- `results_motif_screening/tables/transcript_ctgcagy_summary.csv`
- `results_motif_screening/tables/transcript_structure_proxy_scores.csv`
- `GSE_212630_raw_expression_transposed/*_expression_transposed.csv.gz`
- `GSE_212630_raw_expression_transposed/*_metadata.csv`
- `ENSG_allgene_allcharacter.csv`

Key outputs (under `results_motif_family_nvu_rolesplit/`):
- `tables/motif_family_definitions.csv`
- `tables/tool_candidates_by_motif_ranked.csv`
- `tables/nvu_pseudobulk_ctrl_vs_tdpneg_by_celltype.csv`
- `tables/nvu_role_split_classification.csv`
- `tables/capillary_fire_origin_tool_candidates.csv`
- `tables/pericyte_carrier_machinery_candidates.csv`
- `tables/coupled_fire_carrier_candidates.csv`
- `reports/motif_family_nvu_rolesplit_report.md`
- `figures/motif_family_network.png`
- `figures/rolesplit_scatter_cap_vs_peri.png`

Logic used as evidence:
- Donor-level pseudobulk avoids pseudo-replication.
- "Fire-origin" = genes down in Capillary (Control -> TDPneg), but not down in Pericyte.
- "Carrier" = genes up in Pericyte (Control -> TDPneg), not up in Capillary.
- DDR/RNA surveillance genes (ATM, XRCC5/XRCC6, UPF1, XRN2) are expected to
  fall into the Capillary failure class if they are fire-origin candidates.

### 2) `docs/ids_causal_analysis/scripts/fire_origin_landscape.py`
Purpose:
- In-silico "landscape" scan of GC/G4/CTGCAGY features.
- ATR vs ATM proxy expression analysis in Capillary vs Pericyte.
- Generate a combined report with a provisional verdict.

Key inputs:
- `alltranscript.fasta`
- `completeall_promoters.fasta`
- `results_motif_family_nvu_rolesplit/tables/nvu_role_split_classification.csv`
- `results_motif_family_nvu_rolesplit/tables/capillary_fire_origin_tool_candidates.csv`
- `results_motif_family_nvu_rolesplit/tables/pericyte_carrier_machinery_candidates.csv`
- `results_motif_family_nvu_rolesplit/tables/coupled_fire_carrier_candidates.csv`
- `GSE_212630_raw_expression_transposed/Vasc_Capillary_expression_transposed.csv.gz`
- `GSE_212630_raw_expression_transposed/Vasc_Pericyte_expression_transposed.csv.gz`

Key outputs (under `results_fire_origin_landscape/`):
- `tables/gene_sets_summary.csv`
- `tables/transcript_sequence_features.csv`
- `tables/promoter_sequence_features.csv`
- `tables/atr_atm_proxy_scores_by_sample.csv`
- `tables/proxy_gene_availability.csv`
- `reports/fire_origin_landscape_report.md`

Logic used as evidence:
- ATR proxy scores represent chronic replication stress ("tension").
- ATM proxy scores represent DSB response ("brake").
- The hypothesis expects ATR to be maintained while ATM declines
  in Capillary (Control -> TDPneg). This is a pattern check,
  not a causal proof.
- The script flags that full "Fire vs Background" sequence enrichment
  needs robust gene symbol -> ENSG mapping.

### 3) `docs/ids_causal_analysis/scripts/fire_origin_targeted.py`
Purpose:
- A targeted version focused on candidate sets + pseudobulk stats.
- Performs direct Control vs TDPneg comparison for ATR/ATM/R-loop proxies.
- Generates a verdict (SUPPORT / PARTIAL / NOT_SUPPORTED / NO_DATA).

Key inputs:
- `results_motif_family_nvu_rolesplit/tables/*.csv` (role split + candidates)
- `results_tool_screening_pt_vs_tdp/tables/cross_dataset_overlap_genes.csv` (optional)
- `results_motif_family_nvu_rolesplit/tables/nvu_pseudobulk_ctrl_vs_tdpneg_by_celltype.csv`
- `alltranscript.fasta` (for optional sequence feature calibration)

Key outputs (under `results_fire_origin_landscape/`):
- `tables/target_genes.csv`
- `tables/fire_vs_background_comparison.csv`
- `tables/proxy_analysis_results.csv`
- `reports/fire_origin_landscape_report.md` (overwrites/updates)

Logic used as evidence:
- Calculates proxy mean logFC (Control -> TDPneg) for Capillary and Pericyte.
- If Capillary ATM proxy logFC < 0 while ATR is not lower than ATM,
  verdict becomes SUPPORT. This is used as "ATM brake loss" evidence.

### 4) `docs/ids_causal_analysis/scripts/nvu_culprit_decisive.py`
Purpose:
- Rank NVU cell types as "culprit" (fire origin) using early failures.
- Separate cause-near FAILURE indicators from damage markers.
- Check bootstrap stability and cause -> mediator coupling.

Key inputs:
- `results_nvu_pbmc_alignment/tables/gse212630_nvu_sender_axes_by_donor.csv`

Key outputs (under `results_nvu_culprit_decisive/`):
- `tables/nvu_axes_by_donor_celltype_stage.csv`
- `tables/nvu_axis_stage_changes.csv`
- `tables/nvu_axis_bootstrap_stability.csv`
- `tables/nvu_culprit_celltype_ranking.csv`
- `tables/nvu_mechanism_chain_ranking.csv`
- `figures/nvu_culprit_decisive.png`
- `reports/nvu_culprit_decisive_report.md`

Logic used as evidence:
- Cause-near axes (RNA_Cleanup, DDR) are inverted to FAILURE.
- Earlyness: Control -> TDPneg delta for FAILURE axes.
- Monotonicity: trend across Control/TDPneg/TDPmed/TDPhigh.
- Stability: bootstrap CI excludes zero.
- Mediator coupling: FAILURE axes correlate with EV export.
- The top-ranked cell type is interpreted as "fire origin" in NVU.

## Evidence Chain (What These Scripts Establish)

1) NVU role split + donor-level pseudobulk:
   - Capillary shows early failure signatures (RNA cleanup / DDR down).
   - Pericyte shows carrier machinery upregulation (EV export).

2) Targeted ATM proxy test:
   - In Capillary, Control -> TDPneg ATM proxy decrease is treated as
     "ATM brake loss".
   - ATR proxy is expected to be maintained, supporting the "tension + brake loss"
     pattern rather than global collapse.

3) Culprit ranking:
   - Early failure + monotonic trend + bootstrap stability and mediator coupling
     rank Capillary as the most likely fire-origin cell type.

## Caveats and Limits

- These are observational analyses. "SUPPORT" is pattern-level evidence,
  not definitive causality.
- Staging is cross-sectional (Control/TDPneg/TDPmed/TDPhigh), not longitudinal.
- Several steps depend on gene mapping or proxy gene lists.
- Donor counts can be small in specific cell types.

## Recommended Cross-Checks

1) Verify ATM proxy logFC direction directly in:
   - `docs/ids_causal_analysis/results_fire_origin_landscape/tables/proxy_analysis_results.csv`
2) Confirm Capillary failure class in:
   - `docs/ids_causal_analysis/results_motif_family_nvu_rolesplit/tables/nvu_role_split_classification.csv`
3) Confirm culprit ranking and chain coupling in:
   - `docs/ids_causal_analysis/results_nvu_culprit_decisive/tables/nvu_culprit_celltype_ranking.csv`
   - `docs/ids_causal_analysis/results_nvu_culprit_decisive/tables/nvu_mechanism_chain_ranking.csv`
