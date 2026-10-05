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

from gene_results import collapse_to_genes, fetch_gene_descriptions

@st.cache_data(ttl=86400, show_spinner=False)
def cached_descriptions(genes, species):
    return fetch_gene_descriptions(genes, species=species)

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
        st.session_state['de_results'] = (tables, summary, float(alpha))

if 'de_results' in st.session_state:
    tables, summary, result_alpha = st.session_state['de_results']
    # Keep full IDs internally for gene assignment; clean result outputs.
    output_tables = {}
    for name, table in tables.items():
        output = table.copy()
        if 'transcript' in output:
            output['transcript'] = output['transcript'].astype('string').str.replace(r'_.*$', '', regex=True)
        output_tables[name] = output
    show_gene_names = st.checkbox("Show gene annotation beside transcript", value=True)
    name_lookup = {}
    if show_gene_names and any('transcript' in t for t in output_tables.values()):
        name_species = st.selectbox("Species for full gene names", ['human', 'mouse', 'rat'])
        try:
            symbols = tuple(sorted(set(
                str(gene) for t in output_tables.values() if 'transcript' in t
                for gene in t['transcript'].dropna()
            )))
            with st.spinner("Looking up full gene names…"):
                names = cached_descriptions(symbols, name_species)
            if 'gene_name' in names:
                name_lookup = names.set_index('gene')['gene_name'].fillna('').to_dict()
            for output in output_tables.values():
                if 'transcript' in output:
                    output.insert(
                        output.columns.get_loc('transcript') + 1,
                        'gene annotation',
                        output['transcript'].map(name_lookup).fillna('Annotation unavailable'),
                    )
        except Exception as exc:
            st.warning(f"Full gene names could not be retrieved; gene symbols are still shown. {exc}")
    st.subheader("Comparison summary")
    st.dataframe(summary, use_container_width=True, hide_index=True)
    st.download_button("Export analysis results", results_zip(output_tables, summary),
                       file_name="pairwise_differential_expression_results.zip", mime="application/zip")
    st.subheader("Gene results")
    collapse = st.checkbox("Collapse probes / transcript IDs to one row per gene")
    gene_tables = None
    if collapse:
        st.caption("Keeps one existing result per gene; p-values and adjusted p-values are not recalculated gene-level tests.")
        id_col = st.selectbox("ID column", list(next(iter(tables.values())).columns))
        mode = st.selectbox("How to identify genes", ["Remove numeric probe suffix (CLEC4E_21846 → CLEC4E)", "IDs already contain gene names", "Upload transcript-to-gene mapping"])
        mapping = None
        suffix = None
        ready = True
        if mode.startswith("Remove"):
            suffix = st.text_input("Suffix pattern to remove", r"_[0-9]+$")
        elif mode.startswith("Upload"):
            map_file = st.file_uploader("Mapping CSV or TSV (header required)", type=['csv', 'tsv', 'txt'])
            ready = map_file is not None
            if ready:
                try:
                    m = pd.read_csv(map_file, sep=None, engine='python', dtype=str)
                    mi = st.selectbox("Mapping ID column", list(m.columns))
                    mg = st.selectbox("Mapping gene column", list(m.columns), index=min(1, len(m.columns)-1))
                    if mi == mg:
                        st.error("Choose different mapping columns.")
                        ready = False
                    else:
                        mapping = m.loc[:, [mi, mg]]
                except Exception as exc:
                    st.error(str(exc))
                    ready = False
        method_label = st.selectbox("Representative row", ['Lowest p-value', 'Largest absolute log2 fold change'])
        p_col = st.selectbox("P-value column", ['pvalue', 'padj'])
        strip_version = st.checkbox("Remove ID version suffixes before mapping")
        split_multi = st.checkbox("Split genes separated by ///, semicolon, or comma")
        if ready:
            try:
                gene_tables = {}
                for name, table in tables.items():
                    gene_tables[name], dropped = collapse_to_genes(table, id_col=id_col, mapping=mapping,
                        strip_suffix=suffix, method='min_p' if method_label.startswith('Lowest') else 'max_abs_fc',
                        p_col=p_col, strip_version=strip_version, split_multi=split_multi)
                    if 'transcript' in gene_tables[name]:
                        gene_tables[name]['transcript'] = (
                            gene_tables[name]['transcript'].astype('string')
                            .str.replace(r'_.*$', '', regex=True)
                        )
                        if show_gene_names:
                            gene_tables[name].insert(
                                gene_tables[name].columns.get_loc('transcript') + 1,
                                'gene annotation',
                                gene_tables[name]['transcript'].map(name_lookup).fillna('Annotation unavailable'),
                            )
                    if dropped:
                        st.warning(f"{name}: {dropped} entries without a gene were excluded.")
            except Exception as exc:
                st.error(str(exc))
                gene_tables = None
        if gene_tables is not None:
            species = st.selectbox("Species for gene descriptions", ['human', 'mouse', 'rat'])
            annotate = st.checkbox("Add gene-function descriptions from MyGene.info / NCBI")
            st.caption("Lookup sends gene identifiers and species only. Unmatched or ambiguous genes are labeled; descriptions require internet access.")
            if annotate:
                try:
                    genes = tuple(sorted(set(g for table in gene_tables.values() for g in table['gene'])))
                    with st.spinner("Looking up gene descriptions…"):
                        annotations = cached_descriptions(genes, species)
                    gene_tables = {name: table.merge(annotations, on='gene', how='left', validate='many_to_one') for name, table in gene_tables.items()}
                except Exception as exc:
                    st.warning(f"Descriptions could not be retrieved. Collapsed results are still available. {exc}")
            gene_summary = pd.DataFrame([{'comparison': name, 'n_genes': len(table),
                'n_selected_rows_padj_below_threshold': int(table['padj'].lt(result_alpha).sum())}
                for name, table in gene_tables.items()])
            st.dataframe(gene_summary, hide_index=True)
            st.download_button("Download gene results as ZIP", results_zip(
                {name + '_genes': table for name, table in gene_tables.items()}, gene_summary),
                file_name="collapsed_gene_results.zip", mime="application/zip")
    display_tables = gene_tables if gene_tables is not None else output_tables
    comparison = st.selectbox("Inspect one result", list(display_tables))
    result = display_tables[comparison]
    st.write(f"**{len(result):,} rows**")
    st.dataframe(result.head(500), use_container_width=True, hide_index=True)
