# Database Integration Reference

## 40+ Database API Reference

### Protein & Structure Databases

#### UniProt
```python
import requests

def query_uniprot(gene_name, organism='human'):
    """Query UniProt for protein information."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    query = f"gene:{gene_name} AND organism_name:{organism}"
    params = {
        'query': query,
        'format': 'json',
        'fields': 'accession,id,gene_names,protein_name,sequence,ft_domain,ft_binding,cc_function'
    }

    response = requests.get(base_url, params=params)
    return response.json()

# Example usage
result = query_uniprot('VEGFR2')
# Returns: accession, domains, binding sites, functions
```

#### AlphaFold
```python
def get_alphafold_structure(uniprot_id):
    """Get predicted structure from AlphaFold DB."""
    base_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

    response = requests.get(base_url)
    data = response.json()[0]

    return {
        'pdb_url': data['pdbUrl'],
        'pae_url': data['paeImageUrl'],
        'confidence': data['confidenceVersion']
    }
```

#### PDB
```python
def search_pdb(gene_name, organism='Homo sapiens'):
    """Search PDB for protein structures."""
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                        "operator": "exact_match",
                        "value": organism
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": gene_name
                    }
                }
            ]
        },
        "return_type": "entry"
    }

    response = requests.post(search_url, json=query)
    return response.json()
```

### Pathway Databases

#### KEGG
```python
def kegg_get_pathway(pathway_id):
    """Get KEGG pathway information."""
    # Get pathway entry
    entry_url = f"https://rest.kegg.jp/get/{pathway_id}"
    entry = requests.get(entry_url).text

    # Get genes in pathway
    link_url = f"https://rest.kegg.jp/link/genes/{pathway_id}"
    genes = requests.get(link_url).text

    # Get pathway image
    image_url = f"https://rest.kegg.jp/get/{pathway_id}/image"

    return {
        'entry': entry,
        'genes': genes,
        'image_url': image_url
    }

def kegg_find_pathways(gene):
    """Find pathways containing a gene."""
    url = f"https://rest.kegg.jp/link/pathway/hsa:{gene}"
    return requests.get(url).text
```

#### Reactome
```python
def reactome_enrichment(gene_list):
    """Perform pathway enrichment analysis."""
    url = "https://reactome.org/AnalysisService/identifiers/"

    data = "\n".join(gene_list)
    headers = {"Content-Type": "text/plain"}

    response = requests.post(url, data=data, headers=headers)
    result = response.json()

    # Get top pathways
    pathways = []
    for p in result.get('pathways', []):
        pathways.append({
            'id': p['stId'],
            'name': p['name'],
            'pValue': p['entities']['pValue'],
            'fdr': p['entities']['fdr'],
            'genes_found': p['entities']['found']
        })

    return pathways

def reactome_get_interactors(gene):
    """Get interaction partners from Reactome."""
    url = f"https://reactome.org/ContentService/interactors/static/molecule/{gene}/details"
    return requests.get(url).json()
```

#### STRING
```python
def string_network(genes, species=9606, score=400):
    """Get protein interaction network from STRING."""
    url = "https://string-db.org/api/json/network"

    params = {
        'identifiers': '%0d'.join(genes),
        'species': species,
        'required_score': score,
        'network_type': 'functional'
    }

    response = requests.get(url, params=params)
    return response.json()

def string_enrichment(genes, species=9606):
    """Get functional enrichment from STRING."""
    url = "https://string-db.org/api/json/enrichment"

    params = {
        'identifiers': '%0d'.join(genes),
        'species': species
    }

    response = requests.get(url, params=params)
    return response.json()
```

### Genomics Databases

#### Ensembl
```python
def ensembl_gene_info(gene_symbol, species='human'):
    """Get gene information from Ensembl."""
    server = "https://rest.ensembl.org"

    # Search for gene
    search_url = f"{server}/xrefs/symbol/{species}/{gene_symbol}"
    params = {'content-type': 'application/json'}

    response = requests.get(search_url, params=params)
    gene_ids = response.json()

    if gene_ids:
        gene_id = gene_ids[0]['id']

        # Get full gene info
        gene_url = f"{server}/lookup/id/{gene_id}"
        gene_info = requests.get(gene_url, params=params).json()

        # Get sequence
        seq_url = f"{server}/sequence/id/{gene_id}"
        sequence = requests.get(seq_url, params=params).json()

        return {
            'gene_info': gene_info,
            'sequence': sequence
        }

    return None
```

#### NCBI Gene (Entrez)
```python
from Bio import Entrez

Entrez.email = "your@email.com"

def ncbi_gene_info(gene_symbol, organism='Homo sapiens'):
    """Get gene information from NCBI."""
    # Search for gene
    search_handle = Entrez.esearch(
        db="gene",
        term=f"{gene_symbol}[Gene Name] AND {organism}[Organism]"
    )
    search_results = Entrez.read(search_handle)

    if search_results['IdList']:
        gene_id = search_results['IdList'][0]

        # Fetch gene record
        fetch_handle = Entrez.efetch(
            db="gene",
            id=gene_id,
            rettype="xml"
        )
        gene_record = Entrez.read(fetch_handle)

        return gene_record

    return None
```

#### GEO (Gene Expression Omnibus)
```python
def search_geo(query, dataset_type='gds'):
    """Search GEO for expression datasets."""
    search_handle = Entrez.esearch(
        db="gds",
        term=query,
        retmax=50
    )
    results = Entrez.read(search_handle)

    datasets = []
    for gds_id in results['IdList']:
        summary = Entrez.esummary(db="gds", id=gds_id)
        datasets.append(Entrez.read(summary))

    return datasets
```

### Single-Cell Databases

#### CellxGene Census
```python
import cellxgene_census

def query_cellxgene(gene, cell_type=None, tissue=None):
    """Query CellxGene Census for single-cell data."""
    with cellxgene_census.open_soma() as census:
        # Get human data
        human = census["census_data"]["homo_sapiens"]

        # Build query
        query_filter = f"feature_name == '{gene}'"
        if cell_type:
            query_filter += f" and cell_type == '{cell_type}'"
        if tissue:
            query_filter += f" and tissue == '{tissue}'"

        # Get expression data
        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            var_value_filter=f"feature_name == '{gene}'",
            obs_value_filter=query_filter if cell_type or tissue else None
        )

        return adata
```

### Chemical & Drug Databases

#### ChEMBL
```python
from chembl_webresource_client.new_client import new_client

def chembl_target_activities(target_name):
    """Get bioactivities for a target from ChEMBL."""
    target = new_client.target
    activity = new_client.activity

    # Search for target
    target_query = target.search(target_name)
    targets = list(target_query)

    if targets:
        target_chembl_id = targets[0]['target_chembl_id']

        # Get activities
        activities = activity.filter(
            target_chembl_id=target_chembl_id
        ).filter(
            standard_type__in=['IC50', 'Ki', 'Kd', 'EC50']
        )

        return list(activities)

    return []

def chembl_compound_info(compound_name):
    """Get compound information from ChEMBL."""
    molecule = new_client.molecule

    result = molecule.search(compound_name)
    return list(result)
```

#### DrugBank
```python
def drugbank_search(drug_name):
    """Search DrugBank for drug information."""
    # Note: DrugBank requires API key for full access
    # Public data available through RDF downloads

    url = f"https://go.drugbank.com/drugs/search?query={drug_name}"
    # Returns HTML page - parse or use downloaded data

    # Alternative: Use local DrugBank XML/JSON download
    pass
```

#### PubChem
```python
import pubchempy as pcp

def pubchem_compound(name):
    """Get compound information from PubChem."""
    compounds = pcp.get_compounds(name, 'name')

    if compounds:
        c = compounds[0]
        return {
            'cid': c.cid,
            'iupac_name': c.iupac_name,
            'molecular_formula': c.molecular_formula,
            'molecular_weight': c.molecular_weight,
            'canonical_smiles': c.canonical_smiles,
            'isomeric_smiles': c.isomeric_smiles
        }

    return None

def pubchem_bioassays(target_gene):
    """Get bioassay data for a target."""
    # Search for gene target
    assays = pcp.get_assays(target_gene, 'name')
    return assays
```

### Clinical & Disease Databases

#### ClinicalTrials.gov
```python
def search_clinical_trials(condition=None, intervention=None, status='RECRUITING'):
    """Search ClinicalTrials.gov."""
    base_url = "https://clinicaltrials.gov/api/v2/studies"

    params = {
        'format': 'json',
        'pageSize': 50
    }

    query_parts = []
    if condition:
        query_parts.append(f'CONDITION:({condition})')
    if intervention:
        query_parts.append(f'INTERVENTION:({intervention})')

    params['query.cond'] = condition
    params['query.intr'] = intervention
    params['filter.overallStatus'] = status

    response = requests.get(base_url, params=params)
    return response.json()
```

#### ClinVar
```python
def clinvar_variants(gene):
    """Get clinical variants for a gene."""
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    params = {
        'db': 'clinvar',
        'term': f'{gene}[gene]',
        'retmode': 'json',
        'retmax': 100
    }

    response = requests.get(search_url, params=params)
    result = response.json()

    # Fetch variant details
    variant_ids = result['esearchresult']['idlist']

    if variant_ids:
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            'db': 'clinvar',
            'id': ','.join(variant_ids),
            'rettype': 'vcv',
            'retmode': 'json'
        }
        variants = requests.get(fetch_url, params=fetch_params)
        return variants.json()

    return []
```

#### COSMIC
```python
def cosmic_mutations(gene):
    """Get cancer mutations from COSMIC."""
    # Note: COSMIC requires registration for API access
    # Public data available through download

    # Alternative: Use COSMIC website search
    url = f"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln={gene}"

    # For programmatic access, download COSMIC data files
    # and query locally
    pass
```

#### DisGeNET
```python
def disgenet_associations(gene, api_key):
    """Get gene-disease associations from DisGeNET."""
    base_url = "https://www.disgenet.org/api"

    headers = {
        'Authorization': f'Bearer {api_key}',
        'accept': 'application/json'
    }

    # Get gene associations
    url = f"{base_url}/gda/gene/{gene}"
    response = requests.get(url, headers=headers)

    return response.json()
```

#### Open Targets
```python
def open_targets_associations(gene_symbol):
    """Get target-disease associations from Open Targets."""
    url = "https://api.platform.opentargets.org/api/v4/graphql"

    query = """
    query targetAssociations($ensemblId: String!) {
        target(ensemblId: $ensemblId) {
            id
            approvedSymbol
            associatedDiseases {
                count
                rows {
                    disease {
                        id
                        name
                    }
                    score
                    datasourceScores {
                        id
                        score
                    }
                }
            }
        }
    }
    """

    # First, get Ensembl ID
    search_query = """
    query searchTarget($q: String!) {
        search(queryString: $q, entityNames: ["target"]) {
            hits {
                id
            }
        }
    }
    """

    response = requests.post(url, json={
        'query': search_query,
        'variables': {'q': gene_symbol}
    })

    return response.json()
```

### Literature Databases

#### PubMed
```python
def pubmed_search(query, max_results=100):
    """Search PubMed for publications."""
    search_handle = Entrez.esearch(
        db="pubmed",
        term=query,
        retmax=max_results,
        sort="relevance"
    )
    results = Entrez.read(search_handle)

    # Fetch abstracts
    if results['IdList']:
        fetch_handle = Entrez.efetch(
            db="pubmed",
            id=','.join(results['IdList']),
            rettype="abstract",
            retmode="xml"
        )
        records = Entrez.read(fetch_handle)
        return records

    return []
```

#### bioRxiv/medRxiv
```python
def biorxiv_search(query, server='biorxiv'):
    """Search bioRxiv/medRxiv preprints."""
    base_url = f"https://api.biorxiv.org/details/{server}"

    # Get recent papers (API is date-based)
    from datetime import datetime, timedelta
    end_date = datetime.now().strftime('%Y-%m-%d')
    start_date = (datetime.now() - timedelta(days=365)).strftime('%Y-%m-%d')

    url = f"{base_url}/{start_date}/{end_date}"
    response = requests.get(url)

    results = response.json()

    # Filter by query
    filtered = [
        paper for paper in results.get('collection', [])
        if query.lower() in paper.get('title', '').lower() or
           query.lower() in paper.get('abstract', '').lower()
    ]

    return filtered
```

## Cross-Database ID Mapping

### BioMart (Ensembl)
```python
def biomart_id_mapping(gene_symbols, from_type='hgnc_symbol', to_type='ensembl_gene_id'):
    """Map gene IDs using BioMart."""
    import pybiomart

    dataset = pybiomart.Dataset(
        name='hsapiens_gene_ensembl',
        host='http://www.ensembl.org'
    )

    result = dataset.query(
        attributes=[from_type, to_type],
        filters={from_type: gene_symbols}
    )

    return result
```

### UniProt ID Mapping
```python
def uniprot_id_mapping(ids, from_db='Gene_Name', to_db='UniProtKB'):
    """Map IDs using UniProt ID Mapping service."""
    url = "https://rest.uniprot.org/idmapping/run"

    data = {
        'from': from_db,
        'to': to_db,
        'ids': ','.join(ids)
    }

    response = requests.post(url, data=data)
    job_id = response.json()['jobId']

    # Poll for results
    import time
    while True:
        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        status = requests.get(status_url).json()

        if 'results' in status or 'failedIds' in status:
            break
        time.sleep(1)

    # Get results
    results_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    results = requests.get(results_url).json()

    return results
```

## Unified Query Interface

```python
class DatabaseIntegrator:
    """Unified interface for querying multiple databases."""

    def __init__(self):
        self.cache = {}

    def analyze_target(self, gene, queries):
        """
        Comprehensive target analysis across databases.

        Parameters:
        -----------
        gene : str
            Gene symbol
        queries : dict
            Database-specific queries

        Example:
        --------
        db.analyze_target('VEGFR2', {
            'uniprot': ['sequence', 'domains'],
            'kegg': ['pathways'],
            'string': ['interactions'],
            'clinvar': ['variants']
        })
        """
        results = {}

        for database, query_types in queries.items():
            if database == 'uniprot':
                results['uniprot'] = query_uniprot(gene)
            elif database == 'kegg':
                results['kegg'] = kegg_find_pathways(gene)
            elif database == 'string':
                results['string'] = string_network([gene])
            elif database == 'reactome':
                results['reactome'] = reactome_enrichment([gene])
            elif database == 'clinvar':
                results['clinvar'] = clinvar_variants(gene)
            elif database == 'pubmed':
                results['pubmed'] = pubmed_search(gene)
            # Add more databases...

        return results

    def get_clinical_context(self, gene, queries):
        """Get clinical/disease context for a gene."""
        results = {}

        if 'clinicaltrials' in queries:
            results['trials'] = search_clinical_trials(intervention=gene)
        if 'clinvar' in queries:
            results['variants'] = clinvar_variants(gene)
        if 'disgenet' in queries:
            results['diseases'] = disgenet_associations(gene)

        return results
```
