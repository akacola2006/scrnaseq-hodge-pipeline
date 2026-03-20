"""
scRNAseq Hodge Decomposition Pipeline — Streamlit App
======================================================
Interactive web interface for the Hodge decomposition pipeline.

Usage:
    streamlit run app.py
"""
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

# ── Page config ──────────────────────────────────────────────
st.set_page_config(
    page_title="scRNAseq Hodge Pipeline",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

PIPELINE_ROOT = Path(__file__).resolve().parent
DATA_DIR = PIPELINE_ROOT / "data"
H5AD_DIR = DATA_DIR / "h5ad"
METADATA_DIR = DATA_DIR / "metadata"
RESULTS_DIR = PIPELINE_ROOT / "results"
CONFIG_PATH = PIPELINE_ROOT / "project_config.yaml"


# ── Sidebar navigation ──────────────────────────────────────
st.sidebar.title("scRNAseq Hodge Pipeline")
page = st.sidebar.radio(
    "Navigation",
    ["Home", "Data Setup", "Configuration", "Run Pipeline", "Results"],
    index=0,
)

st.sidebar.markdown("---")
st.sidebar.caption("Discrete Hodge decomposition for upstream cell type & gene identification")


# =====================================================================
# HOME
# =====================================================================
if page == "Home":
    st.title("scRNAseq Hodge Decomposition Pipeline")
    st.markdown("""
    Single-cell RNA-seq解析パイプライン。離散Hodge分解を用いて、疾患や条件に関連する
    **上流セルタイプ**と**上流遺伝子**を同定します。

    ### Pipeline Steps
    1. **Data Loading & QC** — h5adファイルの読み込み、セルタイプフィルタリング
    2. **Normalization & Residuals** — log1p(CPM) + 回帰残差 (GPU加速)
    3. **PCA** — ランダムSVD (GPU加速)
    4. **SPD Covariance** — Ledoit-Wolf shrinkage + Log-Euclidean演算
    5. **Pseudotime** — 拡散マップによる疑似時間構築
    6. **Lane A** — セルタイプレベルHodge分解 → 上流セルタイプ同定
    7. **Bootstrap** — ドナーレベルブートストラップ検証
    8. **Lane B** — 遺伝子レベルHodge分解 → 上流遺伝子同定
    9. **Dual-Mode** — K_N vs sparse k-NN 比較
    10. **Enrichment** — 機能モジュールエンリッチメント

    ### Getting Started
    **Data Setup** → **Configuration** → **Run Pipeline** → **Results** の順に進んでください。
    """)

    col1, col2, col3 = st.columns(3)
    with col1:
        h5ad_files = list(H5AD_DIR.glob("*.h5ad")) if H5AD_DIR.exists() else []
        st.metric("h5ad Files", len(h5ad_files))
    with col2:
        sample_info_exists = (METADATA_DIR / "sample_info.csv").exists()
        st.metric("sample_info.csv", "OK" if sample_info_exists else "Missing")
    with col3:
        run_dirs = sorted(RESULTS_DIR.glob("run_*")) if RESULTS_DIR.exists() else []
        st.metric("Completed Runs", len(run_dirs))


# =====================================================================
# DATA SETUP
# =====================================================================
elif page == "Data Setup":
    st.title("Data Setup")

    st.header("1. h5ad Files")
    st.markdown(f"Place h5ad files in: `{H5AD_DIR}`")

    if H5AD_DIR.exists():
        h5ad_files = sorted(H5AD_DIR.glob("*.h5ad"))
        if h5ad_files:
            st.success(f"{len(h5ad_files)} h5ad files found")
            df_files = pd.DataFrame({
                "File": [f.name for f in h5ad_files],
                "Size (MB)": [f"{f.stat().st_size / 1e6:.1f}" for f in h5ad_files],
            })
            st.dataframe(df_files, use_container_width=True)

            # Preview first file
            if st.checkbox("Preview first h5ad file"):
                try:
                    import anndata as ad
                    adata = ad.read_h5ad(h5ad_files[0], backed="r")
                    st.write(f"**{h5ad_files[0].name}**: {adata.n_obs} cells x {adata.n_vars} genes")
                    st.write("**obs columns:**", list(adata.obs.columns))
                    if len(adata.obs) > 0:
                        st.dataframe(adata.obs.head(5))
                except Exception as e:
                    st.error(f"Error reading h5ad: {e}")
        else:
            st.warning("No h5ad files found. Place your files in the h5ad/ folder.")
    else:
        st.error(f"Directory not found: {H5AD_DIR}")
        if st.button("Create data directories"):
            H5AD_DIR.mkdir(parents=True, exist_ok=True)
            METADATA_DIR.mkdir(parents=True, exist_ok=True)
            st.success("Directories created!")
            st.rerun()

    st.header("2. Sample Metadata")
    sample_info_path = METADATA_DIR / "sample_info.csv"

    if sample_info_path.exists():
        df_info = pd.read_csv(sample_info_path)
        st.success(f"sample_info.csv loaded: {len(df_info)} samples")
        st.dataframe(df_info, use_container_width=True)

        if "condition" in df_info.columns:
            st.write("**Condition counts:**")
            st.bar_chart(df_info["condition"].value_counts())
    else:
        st.warning("sample_info.csv not found.")
        st.markdown("""
        Create a CSV with columns: `donor_id`, `file`, `condition`

        Example:
        | donor_id | file | condition |
        |----------|------|-----------|
        | D001 | sample1.h5ad | Disease |
        | D002 | sample2.h5ad | Control |
        """)

        uploaded = st.file_uploader("Upload sample_info.csv", type="csv")
        if uploaded is not None:
            METADATA_DIR.mkdir(parents=True, exist_ok=True)
            df_upload = pd.read_csv(uploaded)
            df_upload.to_csv(sample_info_path, index=False)
            st.success("sample_info.csv saved!")
            st.dataframe(df_upload)
            st.rerun()

    st.header("3. Reference Data")
    col1, col2 = st.columns(2)
    with col1:
        annot_path = METADATA_DIR / "gene_annotation.csv"
        if annot_path.exists():
            n_genes = sum(1 for _ in open(annot_path)) - 1
            st.success(f"gene_annotation.csv: {n_genes:,} genes")
        else:
            st.info("gene_annotation.csv: not found (optional)")
    with col2:
        modules_path = METADATA_DIR / "functional_modules.json"
        if modules_path.exists():
            with open(modules_path) as f:
                mods = json.load(f)
            n_mods = len(mods.get("modules", mods))
            st.success(f"functional_modules.json: {n_mods} modules")
        else:
            st.info("functional_modules.json: not found (optional)")


# =====================================================================
# CONFIGURATION
# =====================================================================
elif page == "Configuration":
    st.title("Pipeline Configuration")

    if not CONFIG_PATH.exists():
        st.error("project_config.yaml not found!")
        st.stop()

    import yaml

    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    # ── Project ──
    st.header("Project")
    cfg["project_name"] = st.text_input("Project Name", cfg.get("project_name", "my_project"))

    # ── Cell Types ──
    st.header("Cell Types")
    ct_cfg = cfg.get("cell_types", {})

    ct_column = st.text_input("Cell type column name (in h5ad .obs)",
                              ct_cfg.get("column_name", "CellType"))
    ct_cfg["column_name"] = ct_column

    # Auto-detect cell types from h5ad
    ct_list = ct_cfg.get("include", [])
    if st.button("Auto-detect from h5ad"):
        h5ad_files = sorted(H5AD_DIR.glob("*.h5ad"))
        if h5ad_files:
            try:
                import anndata as ad
                adata = ad.read_h5ad(h5ad_files[0], backed="r")
                if ct_column in adata.obs.columns:
                    detected = sorted(adata.obs[ct_column].unique().tolist())
                    ct_list = detected
                    st.success(f"Detected {len(detected)} cell types: {detected}")
                else:
                    st.error(f"Column '{ct_column}' not found. Available: {list(adata.obs.columns)}")
            except Exception as e:
                st.error(f"Error: {e}")

    ct_text = st.text_area("Cell types (one per line)",
                           "\n".join(ct_list) if ct_list else "CellType_A\nCellType_B")
    ct_cfg["include"] = [x.strip() for x in ct_text.strip().split("\n") if x.strip()]

    pt_text = st.text_area("Pseudotime cell types (subset, one per line)",
                           "\n".join(ct_cfg.get("pseudotime_set", ct_cfg["include"][:4])))
    ct_cfg["pseudotime_set"] = [x.strip() for x in pt_text.strip().split("\n") if x.strip()]
    cfg["cell_types"] = ct_cfg

    # ── Conditions ──
    st.header("Conditions")
    cond_cfg = cfg.get("conditions", {})
    cond_cfg["control_label"] = st.text_input("Control label", cond_cfg.get("control_label", "Control"))
    disease_text = st.text_input("Disease labels (comma-separated)",
                                 ", ".join(cond_cfg.get("disease_labels", ["Disease"])))
    cond_cfg["disease_labels"] = [x.strip() for x in disease_text.split(",") if x.strip()]
    cfg["conditions"] = cond_cfg

    # ── Hardware ──
    st.header("Hardware")
    hw_cfg = cfg.get("hardware", {})
    hw_cfg["use_gpu"] = st.checkbox("Use GPU (if available)", hw_cfg.get("use_gpu", True))
    hw_cfg["n_cpu_cores"] = st.number_input("CPU cores", 1, 64, hw_cfg.get("n_cpu_cores", 8))
    cfg["hardware"] = hw_cfg

    # ── Key Parameters ──
    st.header("Key Parameters")
    col1, col2 = st.columns(2)
    with col1:
        pca_cfg = cfg.get("pca", {})
        pca_cfg["k"] = st.number_input("PCA components (k)", 10, 500, pca_cfg.get("k", 100))
        cfg["pca"] = pca_cfg

        la_cfg = cfg.get("lane_a", {})
        la_cfg["n_windows"] = st.number_input("Pseudotime windows", 3, 20, la_cfg.get("n_windows", 8))
        cfg["lane_a"] = la_cfg
    with col2:
        boot_cfg = cfg.get("bootstrap", {})
        boot_cfg["n_iterations"] = st.number_input("Bootstrap iterations", 10, 500,
                                                    boot_cfg.get("n_iterations", 100))
        cfg["bootstrap"] = boot_cfg

        prep_cfg = cfg.get("preprocessing", {})
        prep_cfg["min_cells_per_donor"] = st.number_input("Min cells per donor", 10, 500,
                                                           prep_cfg.get("min_cells_per_donor", 100))
        cfg["preprocessing"] = prep_cfg

    # ── Save ──
    st.markdown("---")
    if st.button("Save Configuration", type="primary"):
        with open(CONFIG_PATH, "w", encoding="utf-8") as f:
            yaml.dump(cfg, f, default_flow_style=False, allow_unicode=True, sort_keys=False)
        st.success("Configuration saved!")


# =====================================================================
# RUN PIPELINE
# =====================================================================
elif page == "Run Pipeline":
    st.title("Run Pipeline")

    # Validation
    h5ad_files = list(H5AD_DIR.glob("*.h5ad")) if H5AD_DIR.exists() else []
    sample_info_exists = (METADATA_DIR / "sample_info.csv").exists()

    if not h5ad_files:
        st.error("No h5ad files found. Go to Data Setup first.")
        st.stop()
    if not sample_info_exists:
        st.error("sample_info.csv not found. Go to Data Setup first.")
        st.stop()

    st.success(f"Ready: {len(h5ad_files)} h5ad files, sample_info.csv found")

    # Step selection
    all_steps = [
        "validate", "load_data", "residuals", "pca", "spd",
        "pseudotime", "lane_a", "bootstrap", "lane_b",
        "gene_hodge", "dual_mode", "enrichment",
    ]

    run_mode = st.radio("Run mode", ["Full Pipeline", "From specific step"])

    start_step = None
    if run_mode == "From specific step":
        start_step = st.selectbox("Start from step", all_steps)

    # GPU status
    try:
        import torch
        if torch.cuda.is_available():
            gpu_name = torch.cuda.get_device_name(0)
            gpu_mem = torch.cuda.get_device_properties(0).total_mem / 1e9
            st.info(f"GPU: {gpu_name} ({gpu_mem:.1f} GB)")
        else:
            st.info("GPU: not available (will use CPU)")
    except ImportError:
        st.info("PyTorch not installed (will use CPU)")

    # Run button
    st.markdown("---")
    if st.button("Run Pipeline", type="primary"):
        cmd = [sys.executable, str(PIPELINE_ROOT / "run_pipeline.py")]
        if start_step:
            cmd += ["--step", start_step]

        st.markdown("### Pipeline Output")
        output_area = st.empty()
        progress_bar = st.progress(0)
        status_text = st.empty()

        log_lines = []

        with subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            cwd=str(PIPELINE_ROOT),
        ) as proc:
            for line in proc.stdout:
                log_lines.append(line.rstrip())

                # Update progress based on step detection
                for i, step_name in enumerate(all_steps):
                    if f"[STEP" in line or step_name.upper() in line.upper():
                        progress_bar.progress(min((i + 1) / len(all_steps), 1.0))
                        status_text.text(f"Running: {step_name}")

                # Show last 30 lines
                display_text = "\n".join(log_lines[-30:])
                output_area.code(display_text, language="text")

            proc.wait()

        progress_bar.progress(1.0)

        if proc.returncode == 0:
            st.success("Pipeline completed successfully!")
            status_text.text("Done!")
        else:
            st.error(f"Pipeline failed (exit code {proc.returncode})")
            status_text.text("Failed")

        # Show full log
        with st.expander("Full log"):
            st.code("\n".join(log_lines), language="text")


# =====================================================================
# RESULTS
# =====================================================================
elif page == "Results":
    st.title("Results")

    if not RESULTS_DIR.exists():
        st.info("No results yet. Run the pipeline first.")
        st.stop()

    run_dirs = sorted(RESULTS_DIR.glob("run_*"), reverse=True)
    if not run_dirs:
        st.info("No completed runs found.")
        st.stop()

    selected_run = st.selectbox(
        "Select run",
        run_dirs,
        format_func=lambda x: x.name,
    )

    run_dir = Path(selected_run)

    # ── Summary ──
    summary_path = run_dir / "summary.json"
    if summary_path.exists():
        with open(summary_path) as f:
            summary = json.load(f)

        st.header("Run Summary")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Cell Types", summary.get("n_celltypes", "?"))
        with col2:
            st.metric("Samples", summary.get("n_samples", "?"))
        with col3:
            st.metric("PCA k", summary.get("pca_k", "?"))
        with col4:
            gf = summary.get("gradient_fraction")
            st.metric("Gradient Fraction", f"{gf:.3f}" if gf else "?")

        if "upstream_ranking" in summary:
            st.subheader("Upstream Cell Type Ranking")
            for i, ct in enumerate(summary["upstream_ranking"], 1):
                st.write(f"**{i}.** {ct}")

    # ── Lane A ──
    lane_a_path = run_dir / "laneA" / "lane_a_summary.json"
    if lane_a_path.exists():
        with open(lane_a_path) as f:
            la = json.load(f)

        st.header("Lane A: Cell-type Upstream Analysis")

        if "phi_scores" in la:
            phi_df = pd.DataFrame([
                {"Cell Type": ct, "phi": val}
                for ct, val in sorted(la["phi_scores"].items(), key=lambda x: -x[1])
            ])
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("phi Scores")
                st.dataframe(phi_df, use_container_width=True)
            with col2:
                st.subheader("phi Ranking")
                st.bar_chart(phi_df.set_index("Cell Type")["phi"])

        # Hodge decomposition pie chart
        gf = la.get("gradient_fraction", 0)
        cf = la.get("curl_fraction", 0)
        hf = la.get("harmonic_fraction", 0)
        if gf + cf + hf > 0:
            st.subheader("Hodge Decomposition")
            import plotly.graph_objects as go
            fig = go.Figure(data=[go.Pie(
                labels=["Gradient", "Curl", "Harmonic"],
                values=[gf, cf, hf],
                marker_colors=["#2ecc71", "#e74c3c", "#95a5a6"],
                textinfo="label+percent",
                hole=0.3,
            )])
            fig.update_layout(height=350, margin=dict(t=30, b=30))
            st.plotly_chart(fig, use_container_width=True)

    # ── Bootstrap ──
    boot_path = run_dir / "bootstrap" / "bootstrap_summary.json"
    if boot_path.exists():
        with open(boot_path) as f:
            boot = json.load(f)

        st.header("Bootstrap Validation")

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Successful iterations", boot.get("n_success", "?"))
        with col2:
            gf_ci = boot.get("gradient_fraction_ci", [None, None])
            if gf_ci[0] is not None:
                st.metric("GF 95% CI", f"[{gf_ci[0]:.3f}, {gf_ci[1]:.3f}]")

        if "top3_reproducibility" in boot:
            repro_df = pd.DataFrame([
                {"Cell Type": ct, "Reproducibility": val}
                for ct, val in sorted(boot["top3_reproducibility"].items(), key=lambda x: -x[1])
            ])
            st.subheader("Top-3 Reproducibility")
            st.bar_chart(repro_df.set_index("Cell Type")["Reproducibility"])

        if "stability_classification" in boot:
            st.subheader("Stability Classification")
            stab_df = pd.DataFrame([
                {"Cell Type": ct, "Status": status}
                for ct, status in boot["stability_classification"].items()
            ])
            st.dataframe(stab_df, use_container_width=True)

    # ── Gene Hodge ──
    gene_scores_path = run_dir / "laneB" / "gene_hodge" / "gene_phi_scores.csv"
    if gene_scores_path.exists():
        st.header("Gene-Level Hodge Results")
        gene_df = pd.read_csv(gene_scores_path)

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total Genes", len(gene_df))
        with col2:
            n_high = (gene_df["classification"] == "High").sum() if "classification" in gene_df.columns else "?"
            st.metric("High-tier", n_high)
        with col3:
            n_med = (gene_df["classification"] == "Medium").sum() if "classification" in gene_df.columns else "?"
            st.metric("Medium-tier", n_med)

        # Top genes table
        st.subheader("Top Upstream Genes")
        n_show = st.slider("Show top N genes", 10, 100, 30)
        st.dataframe(gene_df.head(n_show), use_container_width=True)

        # phi distribution
        st.subheader("phi Score Distribution")
        import plotly.express as px
        if "classification" in gene_df.columns:
            fig = px.histogram(gene_df, x="phi", color="classification",
                              color_discrete_map={"High": "#e74c3c", "Medium": "#f39c12", "Low": "#95a5a6"},
                              nbins=50, barmode="overlay", opacity=0.7)
        else:
            fig = px.histogram(gene_df, x="phi", nbins=50)
        fig.update_layout(height=350)
        st.plotly_chart(fig, use_container_width=True)

    # ── Gene Hodge Summary ──
    gh_summary_path = run_dir / "laneB" / "gene_hodge" / "gene_hodge_summary.json"
    if gh_summary_path.exists():
        with open(gh_summary_path) as f:
            gh = json.load(f)
        st.subheader("Gene Hodge Summary")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Gradient Fraction", f"{gh.get('gradient_fraction', 0):.4f}")
        with col2:
            st.metric("p-value", f"{gh.get('p_value', 1):.4f}")
        with col3:
            st.metric("Flow Mode", gh.get("flow_mode", "?"))

    # ── Dual Mode ──
    dual_path = run_dir / "laneB" / "gene_hodge" / "dual_mode" / "dual_mode_summary.json"
    if dual_path.exists():
        with open(dual_path) as f:
            dual = json.load(f)
        st.header("Dual-Mode: K_N vs Sparse k-NN")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("K_N GF", f"{dual.get('kn_gf', 0):.4f}")
        with col2:
            st.metric("Optimal k", dual.get("optimal_k", "?"))
        with col3:
            st.metric("phi Spearman", f"{dual.get('optimal_phi_rho', 0):.4f}")

        if "sweep" in dual:
            sweep_df = pd.DataFrame(dual["sweep"])
            st.subheader("k-NN Sweep")
            st.dataframe(sweep_df, use_container_width=True)

            fig = px.line(sweep_df, x="k", y="gf", markers=True,
                         labels={"k": "k (nearest neighbors)", "gf": "Gradient Fraction"})
            fig.add_hline(y=dual.get("kn_gf", 0), line_dash="dash",
                         annotation_text="K_N baseline", line_color="red")
            fig.update_layout(height=350)
            st.plotly_chart(fig, use_container_width=True)

    # ── Enrichment ──
    enrichment_dir = run_dir / "laneB" / "gene_hodge" / "enrichment"
    if enrichment_dir.exists():
        for enrich_file in sorted(enrichment_dir.glob("enrichment_*.csv")):
            st.header(f"Enrichment: {enrich_file.stem.replace('enrichment_', '')}")
            enrich_df = pd.read_csv(enrich_file)
            sig_df = enrich_df[enrich_df["significant"] == True] if "significant" in enrich_df.columns else enrich_df
            if len(sig_df) > 0:
                st.success(f"{len(sig_df)} significant modules")
                st.dataframe(sig_df[["module", "n_overlap", "odds_ratio", "fdr"]].head(20),
                            use_container_width=True)

                fig = px.bar(sig_df.head(15), x="module", y="odds_ratio",
                            color="fdr", color_continuous_scale="Reds_r",
                            labels={"odds_ratio": "Odds Ratio", "fdr": "FDR"})
                fig.update_layout(height=400, xaxis_tickangle=-45)
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No significant modules found")

    # ── Random Baseline ──
    rb_path = run_dir / "laneB" / "random_baseline" / "random_baseline_summary.json"
    if rb_path.exists():
        with open(rb_path) as f:
            rb = json.load(f)
        st.header("Random Matrix Baseline")
        comp = rb.get("comparison", {})
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Observed GF", f"{comp.get('observed_gf', 0):.4f}")
        with col2:
            st.metric("Null Mean GF", f"{comp.get('null_mean', 0):.4f}")
        with col3:
            signal = comp.get("signal_above_null", False)
            st.metric("Signal", "YES" if signal else "NO")
        st.info(comp.get("interpretation", ""))

    # ── Directional ──
    dir_path = run_dir / "laneB" / "directional" / "directional_summary.json"
    if dir_path.exists():
        with open(dir_path) as f:
            dir_data = json.load(f)
        st.header("Directional Decomposition (Δ⁺ / Δ⁻)")
        energy = dir_data.get("energy_fraction", {})
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Gain energy (Δ⁺)", f"{energy.get('gain', 0):.1%}")
        with col2:
            st.metric("Loss energy (Δ⁻)", f"{energy.get('loss', 0):.1%}")

        cross = dir_data.get("cross_correlations", {})
        if cross:
            st.write("**Cross-correlations (Spearman):**")
            for k, v in cross.items():
                st.write(f"  {k}: ρ = {v['rho']:.3f}")

    # ── Multi-Transition ──
    mt_path = run_dir / "laneB" / "multi_transition" / "multi_transition_summary.json"
    if mt_path.exists():
        with open(mt_path) as f:
            mt = json.load(f)
        st.header("Multi-Transition Integration")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Valid Transitions", mt.get("n_transitions", "?"))
        with col2:
            st.metric("Mean Concordance", f"{mt.get('mean_concordance', 0):.3f}")
        with col3:
            st.metric("Best GF", f"{mt.get('best_gf', 0):.4f}")

    # ── Two-Axis ──
    ta_path = run_dir / "laneB" / "two_axis" / "two_axis_trajectory.csv"
    if ta_path.exists():
        st.header("2-Axis Model (TRS x MSS)")
        ta_df = pd.read_csv(ta_path)

        fig = px.line(ta_df, x="window", y=["trs_mean", "mss_mean"],
                     labels={"value": "Score", "window": "Pseudotime Window"},
                     markers=True)
        fig.update_layout(height=400, legend_title="Axis")
        st.plotly_chart(fig, use_container_width=True)

        ta_summary_path = run_dir / "laneB" / "two_axis" / "two_axis_summary.json"
        if ta_summary_path.exists():
            with open(ta_summary_path) as f:
                ta_summary = json.load(f)
            lead_lag = ta_summary.get("lead_lag", "?")
            st.metric("Lead/Lag Pattern", lead_lag)

    # ── Download ──
    st.markdown("---")
    st.header("Export")
    if gene_scores_path.exists():
        gene_csv = pd.read_csv(gene_scores_path).to_csv(index=False)
        st.download_button(
            "Download gene phi scores (CSV)",
            gene_csv,
            file_name="gene_phi_scores.csv",
            mime="text/csv",
        )
