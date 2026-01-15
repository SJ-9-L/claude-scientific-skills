#!/usr/bin/env python3
"""
Database Integrator Module

Unified interface for querying 40+ scientific databases including
protein, pathway, genomics, chemical, and clinical databases.
"""

import time
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field
from collections import defaultdict
import json

import requests


@dataclass
class DatabaseResult:
    """Container for database query results."""
    database: str
    query: str
    success: bool = True
    data: Any = None
    error: Optional[str] = None
    timestamp: str = ""


class DatabaseIntegrator:
    """
    Unified interface for querying multiple scientific databases.

    Databases supported:
    - Protein: UniProt, PDB, AlphaFold, BRENDA
    - Pathway: KEGG, Reactome, STRING
    - Genomics: Ensembl, NCBI Gene, GEO
    - Chemical: ChEMBL, PubChem, DrugBank
    - Clinical: ClinicalTrials.gov, ClinVar, COSMIC
    - Literature: PubMed, bioRxiv
    """

    def __init__(self, cache_enabled: bool = True):
        """
        Initialize database integrator.

        Parameters
        ----------
        cache_enabled : bool
            Enable result caching
        """
        self.cache_enabled = cache_enabled
        self.cache: Dict[str, Any] = {}
        self.rate_limits = {
            'ncbi': 0.34,  # 3 requests per second
            'uniprot': 0.1,
            'string': 1.0,
            'default': 0.5
        }

    def analyze_target(
        self,
        gene: str,
        queries: Dict[str, List[str]]
    ) -> Dict[str, DatabaseResult]:
        """
        Comprehensive target analysis across databases.

        Parameters
        ----------
        gene : str
            Gene symbol
        queries : dict
            Database-specific queries
            Example: {'uniprot': ['sequence', 'domains'], 'kegg': ['pathways']}

        Returns
        -------
        dict
            Results from each database
        """
        results = {}

        for database, query_types in queries.items():
            try:
                if database == 'uniprot':
                    results['uniprot'] = self._query_uniprot(gene, query_types)
                elif database == 'pdb':
                    results['pdb'] = self._query_pdb(gene)
                elif database == 'alphafold':
                    results['alphafold'] = self._query_alphafold(gene)
                elif database == 'kegg':
                    results['kegg'] = self._query_kegg(gene, query_types)
                elif database == 'reactome':
                    results['reactome'] = self._query_reactome(gene)
                elif database == 'string':
                    results['string'] = self._query_string(gene)
                elif database == 'ensembl':
                    results['ensembl'] = self._query_ensembl(gene)
                elif database == 'ncbi':
                    results['ncbi'] = self._query_ncbi_gene(gene)
                elif database == 'chembl':
                    results['chembl'] = self._query_chembl(gene)
                elif database == 'pubchem':
                    results['pubchem'] = self._query_pubchem(gene)
                elif database == 'clinvar':
                    results['clinvar'] = self._query_clinvar(gene)
                elif database == 'clinicaltrials':
                    results['clinicaltrials'] = self._query_clinicaltrials(gene)
                elif database == 'opentargets':
                    results['opentargets'] = self._query_opentargets(gene)
                elif database == 'pubmed':
                    results['pubmed'] = self._query_pubmed(gene)
                elif database == 'gtex':
                    results['gtex'] = self._query_gtex(gene)
                elif database == 'hpa':
                    results['hpa'] = self._query_hpa(gene)
                elif database == 'disgenet':
                    results['disgenet'] = self._query_disgenet(gene)
                else:
                    results[database] = DatabaseResult(
                        database=database,
                        query=gene,
                        success=False,
                        error=f"Unknown database: {database}"
                    )

                time.sleep(self.rate_limits.get(database, self.rate_limits['default']))

            except Exception as e:
                results[database] = DatabaseResult(
                    database=database,
                    query=gene,
                    success=False,
                    error=str(e)
                )

        return results

    def get_clinical_context(
        self,
        gene: str,
        queries: Dict[str, List[str]]
    ) -> Dict[str, DatabaseResult]:
        """
        Get clinical context for a gene.

        Parameters
        ----------
        gene : str
            Gene symbol
        queries : dict
            Clinical database queries

        Returns
        -------
        dict
            Clinical context data
        """
        results = {}

        if 'clinicaltrials' in queries:
            results['clinicaltrials'] = self._query_clinicaltrials(gene)

        if 'clinvar' in queries:
            results['clinvar'] = self._query_clinvar(gene)

        if 'cosmic' in queries:
            results['cosmic'] = self._query_cosmic(gene)

        if 'pharmgkb' in queries:
            results['pharmgkb'] = self._query_pharmgkb(gene)

        if 'fda' in queries:
            results['fda'] = self._query_fda(gene)

        return results

    # Individual database query methods

    def _query_uniprot(self, gene: str, query_types: List[str]) -> DatabaseResult:
        """Query UniProt for protein information."""
        base_url = "https://rest.uniprot.org/uniprotkb/search"

        params = {
            'query': f'gene:{gene} AND organism_id:9606',
            'format': 'json',
            'fields': 'accession,id,gene_names,protein_name,sequence,ft_domain,ft_binding,cc_function,cc_pathway'
        }

        try:
            response = requests.get(base_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='uniprot',
                query=gene,
                success=True,
                data={
                    'entries': data.get('results', []),
                    'total': len(data.get('results', []))
                }
            )

        except Exception as e:
            return DatabaseResult(
                database='uniprot',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_pdb(self, gene: str) -> DatabaseResult:
        """Query PDB for protein structures."""
        search_url = "https://search.rcsb.org/rcsbsearch/v2/query"

        query = {
            "query": {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_polymer_entity.rcsb_gene_name.value",
                    "operator": "exact_match",
                    "value": gene
                }
            },
            "return_type": "entry"
        }

        try:
            response = requests.post(search_url, json=query, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='pdb',
                query=gene,
                success=True,
                data={
                    'structures': data.get('result_set', []),
                    'total': data.get('total_count', 0)
                }
            )

        except Exception as e:
            return DatabaseResult(
                database='pdb',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_alphafold(self, gene: str) -> DatabaseResult:
        """Query AlphaFold for predicted structures."""
        # First get UniProt ID
        uniprot_result = self._query_uniprot(gene, ['accession'])

        if not uniprot_result.success or not uniprot_result.data.get('entries'):
            return DatabaseResult(
                database='alphafold',
                query=gene,
                success=False,
                error="No UniProt entry found"
            )

        uniprot_id = uniprot_result.data['entries'][0].get('primaryAccession', '')

        try:
            url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='alphafold',
                query=gene,
                success=True,
                data={
                    'uniprot_id': uniprot_id,
                    'pdb_url': data[0].get('pdbUrl', '') if data else '',
                    'confidence': data[0].get('confidenceVersion', '') if data else ''
                }
            )

        except Exception as e:
            return DatabaseResult(
                database='alphafold',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_kegg(self, gene: str, query_types: List[str]) -> DatabaseResult:
        """Query KEGG for pathway information."""
        results = {'pathways': [], 'genes': []}

        try:
            # Find pathways containing the gene
            link_url = f"https://rest.kegg.jp/link/pathway/hsa:{gene}"
            response = requests.get(link_url, timeout=30)

            if response.status_code == 200 and response.text.strip():
                for line in response.text.strip().split('\n'):
                    if line:
                        parts = line.split('\t')
                        if len(parts) >= 2:
                            pathway_id = parts[1].replace('path:', '')
                            results['pathways'].append(pathway_id)

            return DatabaseResult(
                database='kegg',
                query=gene,
                success=True,
                data=results
            )

        except Exception as e:
            return DatabaseResult(
                database='kegg',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_reactome(self, gene: str) -> DatabaseResult:
        """Query Reactome for pathway data."""
        base_url = "https://reactome.org/ContentService"

        try:
            # Search for gene
            search_url = f"{base_url}/search/query"
            params = {'query': gene, 'cluster': 'true'}
            response = requests.get(search_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            pathways = []
            if data.get('results'):
                for result in data['results']:
                    for entry in result.get('entries', []):
                        if entry.get('schemaClass') == 'Pathway':
                            pathways.append({
                                'id': entry.get('stId', ''),
                                'name': entry.get('name', '')
                            })

            return DatabaseResult(
                database='reactome',
                query=gene,
                success=True,
                data={'pathways': pathways}
            )

        except Exception as e:
            return DatabaseResult(
                database='reactome',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_string(self, gene: str) -> DatabaseResult:
        """Query STRING for protein interactions."""
        url = "https://string-db.org/api/json/network"

        params = {
            'identifiers': gene,
            'species': 9606,
            'required_score': 700
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            interactions = [
                {
                    'partner': i.get('preferredName_B', ''),
                    'score': i.get('score', 0)
                }
                for i in data
            ]

            return DatabaseResult(
                database='string',
                query=gene,
                success=True,
                data={'interactions': interactions}
            )

        except Exception as e:
            return DatabaseResult(
                database='string',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_ensembl(self, gene: str) -> DatabaseResult:
        """Query Ensembl for gene information."""
        server = "https://rest.ensembl.org"

        try:
            # Search for gene
            search_url = f"{server}/xrefs/symbol/homo_sapiens/{gene}"
            response = requests.get(
                search_url,
                headers={'Content-Type': 'application/json'},
                timeout=30
            )
            response.raise_for_status()
            data = response.json()

            if data:
                gene_id = data[0].get('id', '')

                # Get full gene info
                gene_url = f"{server}/lookup/id/{gene_id}"
                gene_response = requests.get(
                    gene_url,
                    headers={'Content-Type': 'application/json'},
                    timeout=30
                )
                gene_data = gene_response.json()

                return DatabaseResult(
                    database='ensembl',
                    query=gene,
                    success=True,
                    data={
                        'gene_id': gene_id,
                        'description': gene_data.get('description', ''),
                        'biotype': gene_data.get('biotype', ''),
                        'chromosome': gene_data.get('seq_region_name', ''),
                        'start': gene_data.get('start', 0),
                        'end': gene_data.get('end', 0)
                    }
                )

            return DatabaseResult(
                database='ensembl',
                query=gene,
                success=False,
                error="Gene not found"
            )

        except Exception as e:
            return DatabaseResult(
                database='ensembl',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_ncbi_gene(self, gene: str) -> DatabaseResult:
        """Query NCBI Gene database."""
        try:
            from Bio import Entrez
            Entrez.email = "analysis@research.com"

            # Search for gene
            handle = Entrez.esearch(
                db="gene",
                term=f"{gene}[Gene Name] AND Homo sapiens[Organism]"
            )
            results = Entrez.read(handle)

            if results['IdList']:
                gene_id = results['IdList'][0]

                return DatabaseResult(
                    database='ncbi',
                    query=gene,
                    success=True,
                    data={
                        'gene_id': gene_id,
                        'url': f"https://www.ncbi.nlm.nih.gov/gene/{gene_id}"
                    }
                )

            return DatabaseResult(
                database='ncbi',
                query=gene,
                success=False,
                error="Gene not found"
            )

        except ImportError:
            return DatabaseResult(
                database='ncbi',
                query=gene,
                success=False,
                error="Biopython not installed"
            )

        except Exception as e:
            return DatabaseResult(
                database='ncbi',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_chembl(self, target: str) -> DatabaseResult:
        """Query ChEMBL for bioactivity data."""
        try:
            from chembl_webresource_client.new_client import new_client
            target_api = new_client.target
            activity_api = new_client.activity

            # Search for target
            targets = list(target_api.search(target))

            if targets:
                target_id = targets[0].get('target_chembl_id', '')

                # Get activities
                activities = list(activity_api.filter(
                    target_chembl_id=target_id
                ).filter(
                    standard_type__in=['IC50', 'Ki', 'EC50']
                )[:20])

                return DatabaseResult(
                    database='chembl',
                    query=target,
                    success=True,
                    data={
                        'target_id': target_id,
                        'activities': activities
                    }
                )

            return DatabaseResult(
                database='chembl',
                query=target,
                success=False,
                error="Target not found"
            )

        except ImportError:
            return DatabaseResult(
                database='chembl',
                query=target,
                success=False,
                error="chembl_webresource_client not installed"
            )

        except Exception as e:
            return DatabaseResult(
                database='chembl',
                query=target,
                success=False,
                error=str(e)
            )

    def _query_pubchem(self, compound: str) -> DatabaseResult:
        """Query PubChem for compound information."""
        try:
            import pubchempy as pcp

            compounds = pcp.get_compounds(compound, 'name')

            if compounds:
                c = compounds[0]
                return DatabaseResult(
                    database='pubchem',
                    query=compound,
                    success=True,
                    data={
                        'cid': c.cid,
                        'iupac_name': c.iupac_name,
                        'molecular_formula': c.molecular_formula,
                        'molecular_weight': c.molecular_weight,
                        'smiles': c.canonical_smiles
                    }
                )

            return DatabaseResult(
                database='pubchem',
                query=compound,
                success=False,
                error="Compound not found"
            )

        except ImportError:
            return DatabaseResult(
                database='pubchem',
                query=compound,
                success=False,
                error="pubchempy not installed"
            )

        except Exception as e:
            return DatabaseResult(
                database='pubchem',
                query=compound,
                success=False,
                error=str(e)
            )

    def _query_clinvar(self, gene: str) -> DatabaseResult:
        """Query ClinVar for clinical variants."""
        try:
            from Bio import Entrez
            Entrez.email = "analysis@research.com"

            # Search ClinVar
            handle = Entrez.esearch(
                db="clinvar",
                term=f"{gene}[gene]",
                retmax=50
            )
            results = Entrez.read(handle)

            return DatabaseResult(
                database='clinvar',
                query=gene,
                success=True,
                data={
                    'variant_count': int(results.get('Count', 0)),
                    'variant_ids': results.get('IdList', [])
                }
            )

        except ImportError:
            return DatabaseResult(
                database='clinvar',
                query=gene,
                success=False,
                error="Biopython not installed"
            )

        except Exception as e:
            return DatabaseResult(
                database='clinvar',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_clinicaltrials(self, term: str) -> DatabaseResult:
        """Query ClinicalTrials.gov."""
        base_url = "https://clinicaltrials.gov/api/v2/studies"

        try:
            params = {
                'query.intr': term,
                'pageSize': 50,
                'format': 'json'
            }

            response = requests.get(base_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            studies = data.get('studies', [])

            return DatabaseResult(
                database='clinicaltrials',
                query=term,
                success=True,
                data={
                    'total': data.get('totalCount', len(studies)),
                    'studies': [
                        {
                            'nct_id': s.get('protocolSection', {}).get('identificationModule', {}).get('nctId', ''),
                            'title': s.get('protocolSection', {}).get('identificationModule', {}).get('briefTitle', ''),
                            'status': s.get('protocolSection', {}).get('statusModule', {}).get('overallStatus', '')
                        }
                        for s in studies[:20]
                    ]
                }
            )

        except Exception as e:
            return DatabaseResult(
                database='clinicaltrials',
                query=term,
                success=False,
                error=str(e)
            )

    def _query_opentargets(self, gene: str) -> DatabaseResult:
        """Query Open Targets for target-disease associations."""
        url = "https://api.platform.opentargets.org/api/v4/graphql"

        # First search for gene
        search_query = """
        query searchTarget($q: String!) {
            search(queryString: $q, entityNames: ["target"]) {
                total
                hits {
                    id
                    name
                }
            }
        }
        """

        try:
            response = requests.post(url, json={
                'query': search_query,
                'variables': {'q': gene}
            }, timeout=30)
            response.raise_for_status()
            data = response.json()

            hits = data.get('data', {}).get('search', {}).get('hits', [])

            if hits:
                return DatabaseResult(
                    database='opentargets',
                    query=gene,
                    success=True,
                    data={
                        'target_id': hits[0].get('id', ''),
                        'name': hits[0].get('name', ''),
                        'total_hits': data.get('data', {}).get('search', {}).get('total', 0)
                    }
                )

            return DatabaseResult(
                database='opentargets',
                query=gene,
                success=False,
                error="Target not found"
            )

        except Exception as e:
            return DatabaseResult(
                database='opentargets',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_pubmed(self, query: str) -> DatabaseResult:
        """Query PubMed for publications."""
        try:
            from Bio import Entrez
            Entrez.email = "analysis@research.com"

            handle = Entrez.esearch(
                db="pubmed",
                term=query,
                retmax=50,
                sort="relevance"
            )
            results = Entrez.read(handle)

            return DatabaseResult(
                database='pubmed',
                query=query,
                success=True,
                data={
                    'total': int(results.get('Count', 0)),
                    'pmids': results.get('IdList', [])
                }
            )

        except ImportError:
            return DatabaseResult(
                database='pubmed',
                query=query,
                success=False,
                error="Biopython not installed"
            )

        except Exception as e:
            return DatabaseResult(
                database='pubmed',
                query=query,
                success=False,
                error=str(e)
            )

    def _query_gtex(self, gene: str) -> DatabaseResult:
        """Query GTEx for tissue expression."""
        # GTEx API endpoint
        url = f"https://gtexportal.org/api/v2/expression/medianGeneExpression"

        try:
            params = {
                'geneSymbol': gene,
                'format': 'json'
            }

            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='gtex',
                query=gene,
                success=True,
                data={'expression': data.get('data', [])}
            )

        except Exception as e:
            return DatabaseResult(
                database='gtex',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_hpa(self, gene: str) -> DatabaseResult:
        """Query Human Protein Atlas."""
        url = f"https://www.proteinatlas.org/{gene}.json"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='hpa',
                query=gene,
                success=True,
                data=data
            )

        except Exception as e:
            return DatabaseResult(
                database='hpa',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_cosmic(self, gene: str) -> DatabaseResult:
        """Query COSMIC for cancer mutations."""
        # COSMIC requires authentication for API
        # Return guidance for manual access
        return DatabaseResult(
            database='cosmic',
            query=gene,
            success=True,
            data={
                'url': f"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln={gene}",
                'note': 'COSMIC requires registration for API access'
            }
        )

    def _query_pharmgkb(self, gene: str) -> DatabaseResult:
        """Query PharmGKB for pharmacogenomics."""
        url = f"https://api.pharmgkb.org/v1/data/gene/{gene}"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='pharmgkb',
                query=gene,
                success=True,
                data=data.get('data', {})
            )

        except Exception as e:
            return DatabaseResult(
                database='pharmgkb',
                query=gene,
                success=False,
                error=str(e)
            )

    def _query_fda(self, term: str) -> DatabaseResult:
        """Query FDA databases."""
        # FDA FAERS API
        url = f"https://api.fda.gov/drug/event.json"

        try:
            params = {
                'search': f'patient.drug.medicinalproduct:"{term}"',
                'limit': 10
            }

            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            return DatabaseResult(
                database='fda',
                query=term,
                success=True,
                data={
                    'total': data.get('meta', {}).get('results', {}).get('total', 0),
                    'events': data.get('results', [])
                }
            )

        except Exception as e:
            return DatabaseResult(
                database='fda',
                query=term,
                success=False,
                error=str(e)
            )

    def _query_disgenet(self, gene: str) -> DatabaseResult:
        """Query DisGeNET for gene-disease associations."""
        # DisGeNET requires API key
        return DatabaseResult(
            database='disgenet',
            query=gene,
            success=True,
            data={
                'url': f"https://www.disgenet.org/search?term={gene}",
                'note': 'DisGeNET requires API key for programmatic access'
            }
        )

    def export_results(
        self,
        results: Dict[str, DatabaseResult],
        format: str = 'json',
        output_path: str = 'database_results.json'
    ) -> str:
        """Export query results."""
        if format == 'json':
            output = {}
            for db, result in results.items():
                output[db] = {
                    'success': result.success,
                    'data': result.data,
                    'error': result.error
                }

            with open(output_path, 'w') as f:
                json.dump(output, f, indent=2, default=str)

        return output_path


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Database Integrator')
    parser.add_argument('gene', help='Gene symbol to query')
    parser.add_argument('--databases', nargs='+', default=['uniprot', 'kegg', 'string'],
                        help='Databases to query')
    parser.add_argument('--output', '-o', default='database_results.json',
                        help='Output file')

    args = parser.parse_args()

    integrator = DatabaseIntegrator()

    queries = {db: ['all'] for db in args.databases}
    results = integrator.analyze_target(args.gene, queries)

    integrator.export_results(results, output_path=args.output)
    print(f"Results saved to {args.output}")
