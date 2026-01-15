# Universal Phenotype Database Reference

## Overview

범용 Phenotype 데이터베이스 - 어떤 유전자/pathway에도 적용 가능한 phenotype 정보 조회 시스템

## Phenotype Data Sources

### 1. Mouse Genome Informatics (MGI)
```python
def query_mgi_phenotypes(gene_symbol: str) -> dict:
    """
    Query MGI for mouse phenotype data.

    Returns:
    - Allele information
    - Phenotype annotations (MP ontology)
    - Genetic backgrounds
    - References
    """
    import requests

    # MGI MouseMine API
    base_url = "https://www.mousemine.org/mousemine/service"

    # Query for gene
    query = f'''
    <query model="genomic" view="Gene.symbol Gene.primaryIdentifier
           Gene.alleles.symbol Gene.alleles.alleleType
           Gene.alleles.phenotypeSummary">
        <constraint path="Gene.symbol" op="=" value="{gene_symbol}"/>
        <constraint path="Gene.organism.shortName" op="=" value="M. musculus"/>
    </query>
    '''

    response = requests.post(f"{base_url}/query/results",
                            data={'query': query, 'format': 'json'})

    return response.json()
```

### 2. International Mouse Phenotyping Consortium (IMPC)
```python
def query_impc_phenotypes(gene_symbol: str) -> dict:
    """
    Query IMPC for standardized knockout phenotypes.

    Returns systematic phenotyping data from:
    - Viability assessment
    - Fertility assessment
    - Body composition
    - Clinical chemistry
    - Behavioral analysis
    - And 20+ other phenotyping pipelines
    """
    import requests

    base_url = "https://www.ebi.ac.uk/mi/impc/solr"

    # Get gene information
    gene_url = f"{base_url}/gene/select"
    params = {
        'q': f'marker_symbol:{gene_symbol}',
        'wt': 'json',
        'rows': 1
    }
    gene_data = requests.get(gene_url, params=params).json()

    if gene_data['response']['numFound'] == 0:
        return None

    gene_info = gene_data['response']['docs'][0]
    mgi_id = gene_info.get('mgi_accession_id')

    # Get phenotype data
    phenotype_url = f"{base_url}/genotype-phenotype/select"
    params = {
        'q': f'marker_accession_id:"{mgi_id}"',
        'wt': 'json',
        'rows': 1000
    }
    phenotype_data = requests.get(phenotype_url, params=params).json()

    return {
        'gene': gene_info,
        'phenotypes': phenotype_data['response']['docs']
    }
```

### 3. Human Phenotype Ontology (HPO)
```python
def query_hpo_phenotypes(gene_symbol: str) -> dict:
    """
    Query HPO for human disease phenotypes.

    Returns:
    - Associated HPO terms
    - Disease associations
    - Phenotype frequencies
    """
    import requests

    base_url = "https://hpo.jax.org/api/hpo"

    # Search for gene
    search_url = f"{base_url}/search"
    params = {'q': gene_symbol, 'category': 'genes'}
    search_result = requests.get(search_url, params=params).json()

    if not search_result.get('genes'):
        return None

    gene_id = search_result['genes'][0]['entrezGeneId']

    # Get gene-phenotype associations
    gene_url = f"{base_url}/gene/{gene_id}"
    gene_data = requests.get(gene_url).json()

    return gene_data
```

### 4. GWAS Catalog
```python
def query_gwas_phenotypes(gene_symbol: str) -> list:
    """
    Query GWAS Catalog for trait associations.

    Returns SNPs and associated traits near the gene.
    """
    import requests

    base_url = "https://www.ebi.ac.uk/gwas/rest/api"

    # Get associations for gene
    url = f"{base_url}/singleNucleotidePolymorphisms/search/findByGene"
    params = {'geneName': gene_symbol}

    response = requests.get(url, params=params)
    return response.json() if response.status_code == 200 else []
```

## Mammalian Phenotype Ontology (MP) Categories

```python
MP_ONTOLOGY_CATEGORIES = {
    'MP:0005386': {
        'name': 'behavior/neurological phenotype',
        'subcategories': [
            'MP:0001392',  # abnormal locomotor behavior
            'MP:0001406',  # abnormal sensory capabilities
            'MP:0002063',  # abnormal learning/memory
        ]
    },
    'MP:0005385': {
        'name': 'cardiovascular system phenotype',
        'subcategories': [
            'MP:0000266',  # abnormal heart morphology
            'MP:0001544',  # abnormal cardiovascular physiology
            'MP:0002127',  # abnormal angiogenesis
            'MP:0005594',  # abnormal vascular development
        ]
    },
    'MP:0005376': {
        'name': 'homeostasis/metabolism phenotype',
        'subcategories': [
            'MP:0001764',  # abnormal blood chemistry
            'MP:0005559',  # abnormal glucose homeostasis
        ]
    },
    'MP:0005378': {
        'name': 'growth/size/body region phenotype',
        'subcategories': [
            'MP:0001265',  # decreased body size
            'MP:0001262',  # increased body size
        ]
    },
    'MP:0005389': {
        'name': 'reproductive system phenotype',
        'subcategories': [
            'MP:0001919',  # abnormal fertility/fecundity
            'MP:0001920',  # abnormal reproductive system morphology
        ]
    },
    'MP:0005381': {
        'name': 'digestive/alimentary phenotype',
        'subcategories': []
    },
    'MP:0005382': {
        'name': 'craniofacial phenotype',
        'subcategories': []
    },
    'MP:0005384': {
        'name': 'cellular phenotype',
        'subcategories': [
            'MP:0000607',  # abnormal cell death
            'MP:0000350',  # abnormal cell proliferation
            'MP:0000353',  # abnormal cell differentiation
        ]
    },
    'MP:0005387': {
        'name': 'immune system phenotype',
        'subcategories': [
            'MP:0002419',  # abnormal innate immunity
            'MP:0002420',  # abnormal adaptive immunity
        ]
    },
    'MP:0010771': {
        'name': 'integument phenotype',
        'subcategories': []
    },
    'MP:0005388': {
        'name': 'respiratory system phenotype',
        'subcategories': []
    },
    'MP:0005367': {
        'name': 'renal/urinary system phenotype',
        'subcategories': []
    },
    'MP:0005390': {
        'name': 'skeleton phenotype',
        'subcategories': []
    },
    'MP:0005377': {
        'name': 'hearing/vestibular/ear phenotype',
        'subcategories': []
    },
    'MP:0005391': {
        'name': 'vision/eye phenotype',
        'subcategories': []
    },
    'MP:0005380': {
        'name': 'embryo phenotype',
        'subcategories': [
            'MP:0001698',  # embryonic lethality
            'MP:0002086',  # abnormal embryonic tissue morphology
        ]
    },
    'MP:0002873': {
        'name': 'normal phenotype',
        'subcategories': []
    },
    'MP:0010768': {
        'name': 'mortality/aging',
        'subcategories': [
            'MP:0011100',  # preweaning lethality
            'MP:0011400',  # premature death
        ]
    },
    'MP:0001186': {
        'name': 'pigmentation phenotype',
        'subcategories': []
    },
    'MP:0002006': {
        'name': 'neoplasm',
        'subcategories': []
    },
    'MP:0003631': {
        'name': 'nervous system phenotype',
        'subcategories': []
    },
    'MP:0005375': {
        'name': 'adipose tissue phenotype',
        'subcategories': []
    },
    'MP:0005369': {
        'name': 'muscle phenotype',
        'subcategories': []
    },
    'MP:0005371': {
        'name': 'limbs/digits/tail phenotype',
        'subcategories': []
    },
    'MP:0005370': {
        'name': 'liver/biliary system phenotype',
        'subcategories': []
    },
    'MP:0005379': {
        'name': 'endocrine/exocrine gland phenotype',
        'subcategories': []
    },
    'MP:0005394': {
        'name': 'taste/olfaction phenotype',
        'subcategories': []
    },
    'MP:0005397': {
        'name': 'hematopoietic system phenotype',
        'subcategories': []
    }
}
```

## Unified Phenotype Query System

```python
class UniversalPhenotypeDB:
    """
    Unified interface for querying phenotype data across multiple databases.
    """

    def __init__(self):
        self.cache = {}

    def get_phenotypes(
        self,
        gene: str,
        species: str = 'all',
        phenotype_categories: list = None
    ) -> dict:
        """
        Get comprehensive phenotype data for any gene.

        Parameters
        ----------
        gene : str
            Gene symbol
        species : str
            'human', 'mouse', 'rat', or 'all'
        phenotype_categories : list, optional
            Filter by MP/HPO categories

        Returns
        -------
        dict
            Unified phenotype data from all sources
        """
        results = {
            'gene': gene,
            'mouse': {},
            'human': {},
            'summary': {}
        }

        # Mouse phenotypes
        if species in ['mouse', 'all']:
            results['mouse']['mgi'] = query_mgi_phenotypes(gene)
            results['mouse']['impc'] = query_impc_phenotypes(gene)

        # Human phenotypes
        if species in ['human', 'all']:
            results['human']['hpo'] = query_hpo_phenotypes(gene)
            results['human']['gwas'] = query_gwas_phenotypes(gene)

        # Generate summary
        results['summary'] = self._summarize_phenotypes(results)

        return results

    def _summarize_phenotypes(self, data: dict) -> dict:
        """Generate phenotype summary across sources."""
        summary = {
            'total_phenotypes': 0,
            'lethality': None,
            'major_systems_affected': [],
            'key_phenotypes': []
        }

        # Check for lethality
        if data['mouse'].get('impc'):
            phenotypes = data['mouse']['impc'].get('phenotypes', [])
            for p in phenotypes:
                if 'lethal' in p.get('mp_term_name', '').lower():
                    summary['lethality'] = {
                        'type': p.get('mp_term_name'),
                        'source': 'IMPC'
                    }
                    break

        return summary

    def compare_ko_phenotypes(
        self,
        genes: list,
        phenotype_focus: str = None
    ) -> dict:
        """
        Compare phenotypes across multiple gene knockouts.

        Parameters
        ----------
        genes : list
            List of gene symbols to compare
        phenotype_focus : str, optional
            Focus on specific phenotype category

        Returns
        -------
        dict
            Comparative phenotype analysis
        """
        comparison = {
            'genes': genes,
            'individual': {},
            'comparison': {
                'shared_phenotypes': [],
                'unique_phenotypes': {},
                'severity_ranking': []
            }
        }

        # Get phenotypes for each gene
        all_phenotypes = {}
        for gene in genes:
            pheno_data = self.get_phenotypes(gene, species='mouse')
            comparison['individual'][gene] = pheno_data
            all_phenotypes[gene] = set()

            # Extract phenotype terms
            if pheno_data['mouse'].get('impc'):
                for p in pheno_data['mouse']['impc'].get('phenotypes', []):
                    all_phenotypes[gene].add(p.get('mp_term_name', ''))

        # Find shared phenotypes
        if len(genes) >= 2:
            shared = all_phenotypes[genes[0]]
            for gene in genes[1:]:
                shared = shared & all_phenotypes[gene]
            comparison['comparison']['shared_phenotypes'] = list(shared)

        # Find unique phenotypes
        for gene in genes:
            others = set()
            for other_gene in genes:
                if other_gene != gene:
                    others = others | all_phenotypes[other_gene]
            comparison['comparison']['unique_phenotypes'][gene] = list(
                all_phenotypes[gene] - others
            )

        return comparison

    def get_phenotype_by_system(
        self,
        gene: str,
        system: str
    ) -> dict:
        """
        Get phenotypes filtered by biological system.

        Parameters
        ----------
        gene : str
            Gene symbol
        system : str
            System name: 'cardiovascular', 'immune', 'nervous', etc.

        Returns
        -------
        dict
            System-specific phenotypes
        """
        system_mp_terms = {
            'cardiovascular': 'MP:0005385',
            'immune': 'MP:0005387',
            'nervous': 'MP:0003631',
            'reproductive': 'MP:0005389',
            'skeletal': 'MP:0005390',
            'metabolic': 'MP:0005376',
            'embryonic': 'MP:0005380',
            'hematopoietic': 'MP:0005397'
        }

        mp_root = system_mp_terms.get(system.lower())
        if not mp_root:
            raise ValueError(f"Unknown system: {system}")

        pheno_data = self.get_phenotypes(gene, species='mouse')

        # Filter by system
        filtered = []
        if pheno_data['mouse'].get('impc'):
            for p in pheno_data['mouse']['impc'].get('phenotypes', []):
                mp_id = p.get('mp_term_id', '')
                # Check if belongs to system (simplified - full implementation would use ontology)
                if system.lower() in p.get('mp_term_name', '').lower():
                    filtered.append(p)

        return {
            'gene': gene,
            'system': system,
            'phenotypes': filtered
        }
```

## Phenotype Template Database

### Common Phenotype Patterns by Pathway

```python
PATHWAY_PHENOTYPE_TEMPLATES = {
    'angiogenesis': {
        'expected_phenotypes': [
            'abnormal angiogenesis',
            'abnormal vascular development',
            'abnormal blood vessel morphology',
            'embryonic lethality during organogenesis',
            'abnormal vasculogenesis',
            'abnormal sprouting angiogenesis'
        ],
        'cell_type_phenotypes': {
            'tip_cell': [
                'absent tip cells',
                'reduced filopodia',
                'impaired migration',
                'abnormal vessel sprouting'
            ],
            'stalk_cell': [
                'reduced proliferation',
                'abnormal lumen formation',
                'vessel hypersprouting'
            ],
            'pericyte': [
                'vascular instability',
                'hemorrhage',
                'increased permeability'
            ]
        }
    },

    'notch_signaling': {
        'expected_phenotypes': [
            'somite development defects',
            'arteriovenous malformation',
            'abnormal lateral inhibition',
            'abnormal cell fate determination'
        ],
        'cell_type_phenotypes': {
            'arterial_ec': ['loss of arterial identity', 'venous marker expression'],
            'venous_ec': ['arterial marker expression']
        }
    },

    'wnt_signaling': {
        'expected_phenotypes': [
            'axis duplication',
            'abnormal limb development',
            'abnormal hair follicle',
            'abnormal stem cell maintenance'
        ]
    },

    'tgfb_signaling': {
        'expected_phenotypes': [
            'embryonic lethality',
            'abnormal immune tolerance',
            'fibrosis',
            'abnormal wound healing'
        ]
    },

    'hedgehog_signaling': {
        'expected_phenotypes': [
            'holoprosencephaly',
            'polydactyly',
            'abnormal neural tube development'
        ]
    },

    'mapk_signaling': {
        'expected_phenotypes': [
            'abnormal cell proliferation',
            'abnormal cell survival',
            'abnormal immune response'
        ]
    },

    'pi3k_akt_signaling': {
        'expected_phenotypes': [
            'abnormal glucose homeostasis',
            'abnormal cell growth',
            'increased susceptibility to cancer'
        ]
    },

    'jak_stat_signaling': {
        'expected_phenotypes': [
            'abnormal cytokine signaling',
            'immunodeficiency',
            'abnormal hematopoiesis'
        ]
    },

    'nfkb_signaling': {
        'expected_phenotypes': [
            'abnormal inflammatory response',
            'immunodeficiency',
            'abnormal B cell development'
        ]
    },

    'hippo_signaling': {
        'expected_phenotypes': [
            'organ overgrowth',
            'abnormal contact inhibition',
            'increased susceptibility to cancer'
        ]
    }
}
```

## Auto-Generated Phenotype Report

```python
def generate_phenotype_report(gene: str, context: str = None) -> str:
    """
    Generate comprehensive phenotype report for any gene.

    Parameters
    ----------
    gene : str
        Gene symbol
    context : str, optional
        Biological context (e.g., 'angiogenesis', 'cancer')

    Returns
    -------
    str
        Markdown-formatted report
    """
    db = UniversalPhenotypeDB()
    data = db.get_phenotypes(gene)

    report = f"""# Phenotype Report: {gene}

## Summary
- **Gene Symbol**: {gene}
- **Data Sources**: MGI, IMPC, HPO, GWAS

## Mouse Phenotypes

### IMPC Systematic Phenotyping
"""

    if data['mouse'].get('impc') and data['mouse']['impc'].get('phenotypes'):
        phenotypes = data['mouse']['impc']['phenotypes']

        # Group by system
        by_system = {}
        for p in phenotypes:
            system = p.get('top_level_mp_term_name', 'Other')
            if system not in by_system:
                by_system[system] = []
            by_system[system].append(p)

        for system, phenos in by_system.items():
            report += f"\n#### {system}\n"
            for p in phenos[:5]:  # Top 5 per system
                report += f"- {p.get('mp_term_name', 'Unknown')}\n"

    report += f"""
## Human Phenotypes

### HPO Disease Associations
"""

    if data['human'].get('hpo'):
        hpo_data = data['human']['hpo']
        diseases = hpo_data.get('diseaseAssoc', [])
        for d in diseases[:10]:
            report += f"- {d.get('diseaseName', 'Unknown')}\n"

    report += f"""
## Clinical Relevance

### GWAS Associations
"""

    return report
```

## Usage Examples

### Example 1: Compare VEGFR2 vs VEGFR3 Knockouts
```python
db = UniversalPhenotypeDB()

comparison = db.compare_ko_phenotypes(
    genes=['KDR', 'FLT4'],  # VEGFR2, VEGFR3
    phenotype_focus='cardiovascular'
)

print("Shared phenotypes:", comparison['comparison']['shared_phenotypes'])
print("VEGFR2-specific:", comparison['comparison']['unique_phenotypes']['KDR'])
print("VEGFR3-specific:", comparison['comparison']['unique_phenotypes']['FLT4'])
```

### Example 2: Get System-Specific Phenotypes
```python
# Get cardiovascular phenotypes for any gene
cv_pheno = db.get_phenotype_by_system('NOTCH1', 'cardiovascular')

# Get immune phenotypes
immune_pheno = db.get_phenotype_by_system('JAK2', 'immune')
```

### Example 3: Generate Report for Novel Target
```python
report = generate_phenotype_report('DLL4', context='angiogenesis')
print(report)
```
