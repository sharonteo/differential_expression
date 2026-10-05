"""Post-process differential-expression rows; never refit gene-level statistics."""
import re
import json
from urllib.request import Request, urlopen
from urllib.parse import urlencode
import pandas as pd


def collapse_to_genes(table, id_col='transcript', mapping=None, strip_suffix=None,
                      method='min_p', p_col='pvalue', fc_col='log2FoldChange',
                      strip_version=False, split_multi=False):
    if method not in ('min_p', 'max_abs_fc'):
        raise ValueError('Unknown selection method.')
    score_col = p_col if method == 'min_p' else fc_col
    for col in (id_col, score_col):
        if col not in table:
            raise ValueError(f'Missing column: {col}')
    d = table.copy()
    ids = d[id_col].astype('string').str.strip()
    if strip_version:
        ids = ids.str.replace(r'\.\d+$', '', regex=True)
    if mapping is not None:
        if mapping.shape[1] < 2:
            raise ValueError('Mapping needs an ID column and a gene column.')
        m = mapping.iloc[:, :2].copy()
        m.columns = ['id', 'gene']
        m['id'] = m['id'].astype('string').str.strip()
        m['gene'] = m['gene'].astype('string').str.strip()
        if strip_version:
            m['id'] = m['id'].str.replace(r'\.\d+$', '', regex=True)
        if m.groupby('id')['gene'].nunique().gt(1).any():
            raise ValueError('Mapping contains conflicting genes for the same ID. Use one row with semicolon-separated genes and enable splitting.')
        d['gene'] = ids.map(m.drop_duplicates('id').set_index('id')['gene'])
    else:
        d['gene'] = ids.str.replace(strip_suffix, '', regex=True) if strip_suffix else ids
    if split_multi:
        d['gene'] = d['gene'].str.split(r'\s*(?:///|;|,)\s*')
        d = d.explode('gene')
    d['gene'] = d['gene'].astype('string').str.strip()
    valid = d['gene'].notna() & ~d['gene'].isin(['', '---', 'NA'])
    dropped = int((~valid).sum())
    d = d.loc[valid].copy()
    if d.empty:
        raise ValueError('No rows could be assigned to genes. Check the mapping and ID format.')
    # Count original rows once per gene, even if a multi-gene label repeats.
    d = d.reset_index().drop_duplicates(['index', 'gene']).drop(columns='index')
    d['n_probes'] = d.groupby('gene')['gene'].transform('size')
    d['_score'] = pd.to_numeric(d[score_col], errors='coerce')
    if method == 'max_abs_fc':
        d['_score'] = -d['_score'].abs()
    d = d.sort_values('_score', kind='stable', na_position='last')
    res = d.drop_duplicates('gene').drop(columns='_score')
    if p_col in res:
        res = res.assign(_p=pd.to_numeric(res[p_col], errors='coerce')).sort_values('_p', kind='stable', na_position='last').drop(columns='_p')
    return res.reset_index(drop=True), dropped


def fetch_gene_descriptions(genes, species='human'):
    """Retrieve exact-ID/symbol matches. Ambiguous matches are left unannotated."""
    genes = list(dict.fromkeys(str(g) for g in genes))
    records = []
    for start in range(0, len(genes), 200):
        batch = genes[start:start + 200]
        request = Request('https://mygene.info/v3/query', data=urlencode({
            'q': ','.join(batch), 'scopes': 'symbol,ensembl.gene,entrezgene',
            'fields': 'symbol,name,summary,entrezgene', 'species': species,
        }).encode(), headers={'Content-Type': 'application/x-www-form-urlencoded'})
        with urlopen(request, timeout=30) as response:
            hits = json.load(response)
        for gene in batch:
            matches = [h for h in hits if h.get('query') == gene and not h.get('notfound')]
            unique = {h.get('_id'): h for h in matches}
            if len(unique) != 1:
                records.append({'gene': gene, 'gene_description': 'Ambiguous match' if unique else 'Description unavailable', 'description_source': ''})
                continue
            hit = next(iter(unique.values()))
            summary = re.sub(r'\s+', ' ', hit.get('summary', '')).strip()
            brief = re.split(r'(?<=[.!?])\s+', summary)[0] if summary else hit.get('name', 'Description unavailable')
            records.append({'gene': gene, 'gene_name': hit.get('name', ''), 'gene_description': brief, 'gene_summary': summary,
                            'description_source': f"https://www.ncbi.nlm.nih.gov/gene/{hit['entrezgene']}" if hit.get('entrezgene') else 'https://mygene.info'})
    return pd.DataFrame(records)
