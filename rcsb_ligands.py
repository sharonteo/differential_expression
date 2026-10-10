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
