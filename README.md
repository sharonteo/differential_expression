# Pairwise Differential Expression (Python + Streamlit)

This project converts the supplied miRNA and transcript DESeq2 R workflows to
Python with PyDESeq2 and adds a Streamlit interface.

## Run the app

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
pip install -r requirements.txt
streamlit run streamlit_app.py
```

Before running locally, copy `.streamlit/secrets.toml.example` to
`.streamlit/secrets.toml` and replace the example value with a strong password.
The real `secrets.toml` file is excluded from Git and must never be committed.

## Deploy to Streamlit Community Cloud

1. Create a GitHub repository and upload the project files. Do not upload the
   `.venv` directory or `.streamlit/secrets.toml`.
2. Sign in at `share.streamlit.io` with GitHub and select **Create app**.
3. Select the repository, branch, and `streamlit_app.py` as the entrypoint.
4. In **Advanced settings**, select Python 3.11 and add this secret:

   ```toml
   APP_PASSWORD = "replace-with-a-strong-private-password"
   ```

5. Deploy the app and test the password screen before sharing its URL.

For additional access control, keep the GitHub repository private. A shared
password protects the app interface but does not provide separate user accounts
or an access audit trail.

## Run from the command line

```bash
python miRNA_pairwise_DESeq2.py SP0497_gene_counts.csv
python transcript_pairwise_DESeq2.py transcript_counts.csv
```

The input must be a raw-count CSV with feature IDs in the first column and one
sample per remaining column. A sample named `LPS_AFX2781_2` is assigned to the
group `LPS_AFX2781`.

The summary is saved as `_summary_all_comparisons.csv`. Each detailed pairwise
result is saved as a tab-delimited text file such as
`AFX220_vs_AFX220HCl.txt`.

## Important statistical note

PyDESeq2 reproduces the core DESeq2 workflow in Python, but results may differ
slightly from R DESeq2. This Python version reports unshrunk log2 fold changes
because PyDESeq2's coefficient-based apeGLM shrinkage does not support every
arbitrary pairwise contrast from one fitted multi-group model. For regulated
or publication-critical work requiring exact R DESeq2/ashr output, retain the
R engine and use Python/Streamlit only as the interface.

## Collapse transcript / probe results to genes

After analysis, enable **Collapse probes / transcript IDs to one row per gene**.
Choose embedded symbols (e.g. `CLEC4E_21846` becomes `CLEC4E`), existing gene
names, or a headered CSV/TSV mapping. Select its ID and gene columns. Ensembl
transcript IDs require a mapping; removing suffixes alone does not identify
their genes. Version removal and splitting multi-gene labels are optional.

The representative row defaults to the lowest p-value (`pvalue` or `padj`);
alternatively select the largest absolute log2 fold change. Missing scores rank
last and ties keep original input order. `n_probes` counts the original rows
represented by each gene. Unmapped / blank genes are excluded with a warning;
conflicting mappings are rejected. Original results remain downloadable.

Enable gene descriptions and select human, mouse, or rat. MyGene.info retrieves
published gene summaries; the first sentence is shown as `gene_description`,
with the full text in `gene_summary` and an NCBI link in `description_source`.
If no summary exists, the gene name is used. Missing / ambiguous matches are
labeled rather than guessed. Lookups send gene identifiers and species only,
require internet access, and are cached for 24 hours. Lookup failure does not
prevent downloading collapsed results. API reference: https://docs.mygene.info/en/latest/doc/quick_start.html

**Statistical interpretation:** this is representative-row selection after
transcript-level analysis, not a gene-level DESeq2 fit. Retained p-values and
adjusted p-values are the original transcript/probe statistics; choosing the
smallest p-value across probes can favor genes with more probes. For formal
gene-level inference, aggregate appropriate raw counts to genes before fitting.
