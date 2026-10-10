"""RCSB ligands present in structures with an exact gene/species match."""
import json
from urllib.parse import urlencode
from urllib.request import urlopen

TAXONOMY = {'human': 9606, 'mouse': 10090, 'rat': 10116}


def _get_json(url, parameters):
    with urlopen(url + '?' + urlencode(parameters), timeout=30) as response:
        if response.status == 204:
            return {}
        return json.load(response)


def fetch_ligands(gene, species='human'):
    """Return unique ligand IDs and source structures; failures propagate for retry."""
    query = {
        'query': {'type': 'terminal', 'service': 'text', 'parameters': {
            'attribute': 'rcsb_entity_source_organism.rcsb_gene_name.value',
            'operator': 'exact_match', 'value': str(gene)}},
        'return_type': 'polymer_entity', 'request_options': {'return_all_hits': True},
    }
    hits = _get_json('https://search.rcsb.org/rcsbsearch/v2/query',
                     {'json': json.dumps(query)}).get('result_set', [])
    entity_ids = sorted({hit['identifier'] for hit in hits})
    ligands, structures = set(), set()
    for start in range(0, len(entity_ids), 50):
        graphql = '''{ polymer_entities(entity_ids: %s) {
            rcsb_id rcsb_entity_source_organism {
                ncbi_taxonomy_id rcsb_gene_name { value }
            } entry { rcsb_id nonpolymer_entities {
                pdbx_entity_nonpoly { comp_id }
            } }
        } }''' % json.dumps(entity_ids[start:start + 50])
        response = _get_json('https://data.rcsb.org/graphql', {'query': graphql})
        if response.get('errors'):
            raise RuntimeError('RCSB data lookup failed')
        for entity in response['data']['polymer_entities'] or []:
            if not entity:
                raise RuntimeError('RCSB returned an incomplete entity')
            sources = entity.get('rcsb_entity_source_organism') or []
            if not any(source.get('ncbi_taxonomy_id') == TAXONOMY[species]
                       and str(gene).upper() in {
                           item['value'].upper() for item in source.get('rcsb_gene_name') or []
                       } for source in sources):
                continue
            entry = entity['entry']
            structures.add(entry['rcsb_id'])
            for molecule in entry.get('nonpolymer_entities') or []:
                component = molecule.get('pdbx_entity_nonpoly') or {}
                if component.get('comp_id'):
                    ligands.add(component['comp_id'])
    return {'ligands': sorted(ligands), 'structures': sorted(structures)}


def ligand_links(record):
    if not record['structures']:
        return 'No matching structure'
    if not record['ligands']:
        return 'No ligands listed'
    return ' · '.join(f'[{ligand}](https://www.rcsb.org/ligand/{ligand})'
                      for ligand in record['ligands'])


def _request_json(request):
    """Retry temporary API failures without caching failures as missing matches."""
    from urllib.error import HTTPError, URLError
    import time
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as exc:
            if exc.code not in (429, 500, 502, 503, 504) or attempt == 2:
                raise
            try:
                delay = min(10, max(1, float(exc.headers.get('Retry-After', 2 ** attempt))))
            except (TypeError, ValueError):
                delay = 2 ** attempt
            time.sleep(delay)
        except (URLError, TimeoutError):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def fetch_pubchem(ligand):
    """Resolve CCD OpenEye canonical SMILES using a form-encoded POST.

    SMILES may contain slashes and other reserved URL characters, so they
    belong in the request body rather than the URL path.
    """
    from urllib.parse import quote
    from urllib.request import Request
    from urllib.error import HTTPError
    component = _request_json(
        'https://data.rcsb.org/rest/v1/core/chemcomp/' + quote(ligand, safe=''))
    smiles = next((item['descriptor'] for item in
                   component.get('pdbx_chem_comp_descriptor', [])
                   if item.get('type') == 'SMILES_CANONICAL'
                   and item.get('program') == 'OpenEye OEToolkits'), None)
    if not smiles:
        return {'smiles': None, 'compounds': []}
    if '*' in smiles:
        return {'smiles': smiles, 'compounds': [],
                'reason': ('Unknown atom or ion — PubChem match not applicable'
                           if ligand.upper() == 'UNX' else
                           'Undefined atoms in SMILES — PubChem match not applicable')}
    request = Request(
        'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/'
        'property/Title,MolecularFormula,MolecularWeight/JSON',
        data=urlencode({'smiles': smiles}).encode(),
        headers={'Content-Type': 'application/x-www-form-urlencoded'})
    try:
        result = _request_json(request)
    except HTTPError as exc:
        if exc.code == 404:
            return {'smiles': smiles, 'compounds': []}
        raise
    return {'smiles': smiles,
            'compounds': result['PropertyTable']['Properties']}


def pubchem_links(ligand, record):
    if record.get('reason'):
        return f"{ligand}: {record['reason']}"
    if not record['smiles']:
        return f'{ligand}: OpenEye canonical SMILES unavailable'
    if not record['compounds']:
        return f'{ligand}: No PubChem match'
    return ' · '.join(
        f"{ligand}: [{item.get('Title', 'Compound')} (CID {item['CID']})]"
        f"(https://pubchem.ncbi.nlm.nih.gov/compound/{item['CID']})"
        for item in record['compounds'])
