"""Streamlit interface for pairwise miRNA/transcript differential expression."""

import hmac

import streamlit as st


st.set_page_config(page_title="Pairwise Differential Expression", layout="wide")


def require_password() -> None:
    """Stop execution until the shared app password is verified."""
    try:
        expected_password = str(st.secrets["APP_PASSWORD"])
    except (KeyError, FileNotFoundError):
        st.error(
            "APP_PASSWORD is not configured. Add it to Streamlit Secrets "
            "before using this app."
        )
        st.stop()

    if st.session_state.get("authenticated", False):
        return

    st.title("Pairwise Differential Expression")
    st.caption("Enter the app password to continue.")
    with st.form("login_form", clear_on_submit=True):
        entered_password = st.text_input("Password", type="password")
        submitted = st.form_submit_button("Sign in", type="primary")

    if submitted:
        if hmac.compare_digest(entered_password, expected_password):
            st.session_state["authenticated"] = True
            st.rerun()
        else:
            st.error("Incorrect password.")
    st.stop()


require_password()

# Load the scientific stack only after authentication.
import pandas as pd

from deseq2_pairwise import (
    AnalysisConfig,
    available_comparisons,
    fit_dataset,
    infer_groups,
    read_count_matrix,
    results_zip,
    run_pairwise,
)


st.title("Pairwise Differential Expression")
st.caption("DESeq2-style analysis in Python using PyDESeq2")

with st.sidebar:
    st.header("Analysis settings")
    data_type = st.radio("Count type", ["miRNA", "Transcript"])
    alpha = st.number_input("Adjusted p-value threshold", 0.001, 0.20, 0.05, 0.01)
    min_count = st.number_input("Minimum raw count", 1, 1000, 10)
    if data_type == "miRNA":
        min_samples = st.number_input("Minimum samples meeting count", 1, 100, 2)
    else:
        min_samples = None
        st.caption("Transcript mode uses the smallest detected group size.")
    cpus = st.number_input("CPU workers", 1, 4, 2)

uploaded = st.file_uploader("Upload a raw-count CSV", type="csv")
st.info(
    "Expected format: feature IDs in the first column; sample names across the "
    "remaining columns. Names such as LPS_DMSO_1 and LPS_DMSO_2 are assigned "
    "to group LPS_DMSO."
)

if uploaded is not None:
    try:
        counts = read_count_matrix(uploaded)
        metadata = infer_groups(counts.columns)
    except Exception as exc:
        st.error(str(exc))
        st.stop()

    left, right = st.columns(2)
    with left:
        st.subheader("Uploaded matrix")
        st.metric("Features", f"{len(counts):,}")
        st.metric("Samples", f"{counts.shape[1]:,}")
        st.dataframe(counts.head(10), use_container_width=True)
    with right:
        st.subheader("Detected groups")
        group_table = metadata.value_counts().sort_index().rename("samples").to_frame()
        st.dataframe(group_table, use_container_width=True)

    pairs = available_comparisons(metadata)
    labels = [f"{a} vs {b}" for a, b in pairs]
    chosen_labels = st.multiselect(
        "Comparisons to run",
        labels,
        default=labels,
        help="A positive log2FoldChange means higher expression in the first group.",
    )
    chosen = [pair for pair, label in zip(pairs, labels) if label in chosen_labels]
    st.caption(f"{len(chosen)} comparison(s) selected.")

    if st.button("Run analysis", type="primary", disabled=not chosen):
        config = AnalysisConfig(
            feature_label="miRNA" if data_type == "miRNA" else "transcript",
            min_count=int(min_count),
            min_samples=int(min_samples) if min_samples is not None else None,
            alpha=float(alpha),
            n_cpus=int(cpus),
        )
        progress_bar = st.progress(0, text="Filtering counts and fitting model…")
        try:
            dds, retained = fit_dataset(counts, metadata, config)

            def update(done: int, total: int, name: str) -> None:
                progress_bar.progress(done / total, text=f"Completed {name}")

            tables, summary = run_pairwise(dds, chosen, config, update)
        except Exception as exc:
            progress_bar.empty()
            st.exception(exc)
            st.stop()

        progress_bar.progress(1.0, text="Analysis complete")
        st.success(f"Retained {retained:,} of {len(counts):,} features after filtering.")
        st.subheader("Comparison summary")
        st.dataframe(summary, use_container_width=True, hide_index=True)

        comparison = st.selectbox("Inspect one result", list(tables))
        result = tables[comparison]
        significant = result[result["padj"].lt(alpha).fillna(False)]
        st.write(f"**{len(significant):,} significant features** at adjusted p < {alpha:g}")
        st.dataframe(result.head(500), use_container_width=True, hide_index=True)
        st.download_button(
            "Download all results as ZIP",
            results_zip(tables, summary),
            file_name="pairwise_differential_expression_results.zip",
            mime="application/zip",
        )
