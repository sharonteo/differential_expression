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
