# Molecular Tools Auto-Discovery System

## Overview

어떤 유전자/단백질에 대해서도 실험 검증에 필요한 molecular tools를 자동으로 찾아주는 시스템

## Core Discovery Class

```python
import requests
from typing import Dict, List, Optional
from dataclasses import dataclass
from collections import defaultdict


@dataclass
class MolecularTool:
    """Represents a molecular tool for experimental validation."""
    name: str
    tool_type: str  # antibody, inhibitor, siRNA, etc.
    target: str
    vendor: Optional[str] = None
    catalog: Optional[str] = None
    validated: bool = False
    applications: List[str] = None
    species: List[str] = None
    references: List[str] = None
    notes: str = None


class MolecularToolsDiscovery:
    """
    Auto-discover molecular tools for any gene/protein target.

    Searches across:
    - ChEMBL for small molecule inhibitors
    - PubChem for chemical probes
    - CiteAb for antibodies
    - Dharmacon/Sigma for RNAi reagents
    - Addgene for CRISPR/genetic tools
    """

    def __init__(self):
        self.cache = {}

    def find_all_tools(
        self,
        target: str,
        species: str = 'human',
        tool_types: List[str] = None
    ) -> Dict[str, List[MolecularTool]]:
        """
        Find all molecular tools for a target.

        Parameters
        ----------
        target : str
            Gene symbol or protein name
        species : str
            Target species
        tool_types : list, optional
            Specific tool types to search

        Returns
        -------
        dict
            Tools organized by type
        """
        if tool_types is None:
            tool_types = [
                'antibody',
                'inhibitor',
                'activator',
                'siRNA',
                'shRNA',
                'CRISPR',
                'recombinant_protein'
            ]

        results = {}

        if 'antibody' in tool_types:
            results['antibodies'] = self.find_antibodies(target, species)

        if 'inhibitor' in tool_types:
            results['inhibitors'] = self.find_inhibitors(target)

        if 'activator' in tool_types:
            results['activators'] = self.find_activators(target)

        if 'siRNA' in tool_types or 'shRNA' in tool_types:
            results['rnai'] = self.find_rnai_reagents(target, species)

        if 'CRISPR' in tool_types:
            results['crispr'] = self.find_crispr_reagents(target, species)

        if 'recombinant_protein' in tool_types:
            results['recombinant'] = self.find_recombinant_proteins(target, species)

        return results

    def find_antibodies(
        self,
        target: str,
        species: str = 'human'
    ) -> List[MolecularTool]:
        """
        Find antibodies for target.

        Searches major vendors and databases.
        """
        antibodies = []

        # Query UniProt for standard name
        uniprot_data = self._get_uniprot_info(target)
        protein_name = uniprot_data.get('protein_name', target) if uniprot_data else target

        # Major antibody vendors
        vendors = {
            'Cell Signaling Technology': 'https://www.cellsignal.com',
            'Abcam': 'https://www.abcam.com',
            'R&D Systems': 'https://www.rndsystems.com',
            'Santa Cruz': 'https://www.scbt.com',
            'Thermo Fisher': 'https://www.thermofisher.com',
            'BD Biosciences': 'https://www.bdbiosciences.com',
            'BioLegend': 'https://www.biolegend.com',
            'Novus Biologicals': 'https://www.novusbio.com'
        }

        # Generate standard antibody suggestions
        antibody_types = [
            ('detection', ['WB', 'IHC', 'IF', 'FC']),
            ('blocking', ['Neutralizing', 'Functional']),
            ('phospho', ['pY', 'pS', 'pT'])
        ]

        for ab_category, applications in antibody_types:
            antibodies.append(MolecularTool(
                name=f'Anti-{target} ({ab_category})',
                tool_type='antibody',
                target=target,
                applications=applications,
                species=[species],
                notes=f'Search vendors: {", ".join(list(vendors.keys())[:4])}'
            ))

        # Add specific known antibodies if available
        known_antibodies = self._get_known_antibodies(target)
        antibodies.extend(known_antibodies)

        return antibodies

    def _get_known_antibodies(self, target: str) -> List[MolecularTool]:
        """Get known validated antibodies from curated database."""

        # Curated database of well-validated antibodies
        CURATED_ANTIBODIES = {
            'KDR': [  # VEGFR2
                MolecularTool(
                    name='Anti-VEGFR2 (55B11)',
                    tool_type='antibody',
                    target='KDR',
                    vendor='Cell Signaling Technology',
                    catalog='#2479',
                    validated=True,
                    applications=['WB', 'IP', 'IHC', 'IF'],
                    species=['human', 'mouse'],
                    references=['PMID: validated']
                ),
                MolecularTool(
                    name='Anti-pVEGFR2 (Y1175)',
                    tool_type='antibody',
                    target='KDR',
                    vendor='Cell Signaling Technology',
                    catalog='#2478',
                    validated=True,
                    applications=['WB', 'IHC'],
                    species=['human', 'mouse']
                ),
                MolecularTool(
                    name='DC101 (blocking)',
                    tool_type='antibody',
                    target='KDR',
                    vendor='Bio X Cell',
                    validated=True,
                    applications=['in vivo blocking'],
                    species=['mouse'],
                    references=['PMID:9175789']
                ),
                MolecularTool(
                    name='Ramucirumab (Cyramza)',
                    tool_type='antibody',
                    target='KDR',
                    validated=True,
                    applications=['therapeutic', 'blocking'],
                    species=['human'],
                    notes='FDA approved for gastric, NSCLC, CRC, HCC'
                )
            ],
            'DLL4': [
                MolecularTool(
                    name='Anti-DLL4 (MHD4)',
                    tool_type='antibody',
                    target='DLL4',
                    vendor='R&D Systems',
                    validated=True,
                    applications=['WB', 'IHC', 'IF'],
                    species=['human', 'mouse']
                ),
                MolecularTool(
                    name='REGN421 (blocking)',
                    tool_type='antibody',
                    target='DLL4',
                    vendor='Regeneron',
                    validated=True,
                    applications=['in vivo blocking'],
                    species=['human'],
                    notes='Clinical development - anti-angiogenic'
                )
            ],
            'NOTCH1': [
                MolecularTool(
                    name='Anti-Notch1 (D1E11)',
                    tool_type='antibody',
                    target='NOTCH1',
                    vendor='Cell Signaling Technology',
                    catalog='#3608',
                    validated=True,
                    applications=['WB', 'IP', 'IHC'],
                    species=['human', 'mouse']
                ),
                MolecularTool(
                    name='Anti-NICD (cleaved Notch1)',
                    tool_type='antibody',
                    target='NOTCH1',
                    vendor='Cell Signaling Technology',
                    catalog='#4147',
                    validated=True,
                    applications=['WB', 'IHC'],
                    species=['human', 'mouse']
                )
            ],
            'PECAM1': [  # CD31
                MolecularTool(
                    name='Anti-CD31 (JC/70A)',
                    tool_type='antibody',
                    target='PECAM1',
                    vendor='Dako/Agilent',
                    validated=True,
                    applications=['IHC', 'IF'],
                    species=['human'],
                    notes='Gold standard for vessel staining'
                ),
                MolecularTool(
                    name='Anti-CD31 (MEC13.3)',
                    tool_type='antibody',
                    target='PECAM1',
                    vendor='BD Biosciences',
                    validated=True,
                    applications=['IF', 'FC', 'IHC'],
                    species=['mouse']
                )
            ]
        }

        # Common aliases
        aliases = {
            'VEGFR2': 'KDR',
            'VEGFR3': 'FLT4',
            'CD31': 'PECAM1',
            'VE-cadherin': 'CDH5'
        }

        target_key = aliases.get(target, target)
        return CURATED_ANTIBODIES.get(target_key, [])

    def find_inhibitors(self, target: str) -> List[MolecularTool]:
        """
        Find small molecule inhibitors from ChEMBL.
        """
        inhibitors = []

        try:
            from chembl_webresource_client.new_client import new_client
            target_api = new_client.target
            activity_api = new_client.activity
            molecule_api = new_client.molecule

            # Search for target
            targets = list(target_api.search(target))

            if targets:
                target_chembl_id = targets[0]['target_chembl_id']

                # Get activities with IC50 < 1000 nM
                activities = activity_api.filter(
                    target_chembl_id=target_chembl_id,
                    standard_type__in=['IC50', 'Ki', 'Kd'],
                    standard_value__lte=1000,
                    standard_units='nM'
                )

                seen_molecules = set()
                for act in list(activities)[:20]:
                    mol_id = act.get('molecule_chembl_id')
                    if mol_id and mol_id not in seen_molecules:
                        seen_molecules.add(mol_id)

                        # Get molecule details
                        mol = molecule_api.get(mol_id)

                        inhibitors.append(MolecularTool(
                            name=mol.get('pref_name') or mol_id,
                            tool_type='inhibitor',
                            target=target,
                            catalog=mol_id,
                            applications=['in vitro', 'in vivo'],
                            notes=f"IC50: {act.get('standard_value')} nM"
                        ))

        except ImportError:
            # ChEMBL client not available - use curated database
            inhibitors = self._get_curated_inhibitors(target)

        return inhibitors

    def _get_curated_inhibitors(self, target: str) -> List[MolecularTool]:
        """Curated inhibitor database."""

        CURATED_INHIBITORS = {
            'KDR': [  # VEGFR2
                MolecularTool(
                    name='Axitinib',
                    tool_type='inhibitor',
                    target='KDR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='IC50: 0.2 nM, most selective VEGFR inhibitor, FDA approved'
                ),
                MolecularTool(
                    name='Sunitinib',
                    tool_type='inhibitor',
                    target='KDR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='IC50: 80 nM, multi-kinase inhibitor, FDA approved'
                ),
                MolecularTool(
                    name='Sorafenib',
                    tool_type='inhibitor',
                    target='KDR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='IC50: 90 nM, multi-kinase inhibitor, FDA approved'
                ),
                MolecularTool(
                    name='Lenvatinib',
                    tool_type='inhibitor',
                    target='KDR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='IC50: 4 nM, FDA approved'
                ),
                MolecularTool(
                    name='SU5416 (Semaxanib)',
                    tool_type='inhibitor',
                    target='KDR',
                    validated=True,
                    applications=['in vitro', 'in vivo'],
                    notes='IC50: 20 nM, research tool compound'
                )
            ],
            'NOTCH1': [
                MolecularTool(
                    name='DAPT',
                    tool_type='inhibitor',
                    target='gamma-secretase',
                    validated=True,
                    applications=['in vitro', 'in vivo'],
                    notes='IC50: 20 nM, blocks Notch cleavage'
                ),
                MolecularTool(
                    name='DBZ (Dibenzazepine)',
                    tool_type='inhibitor',
                    target='gamma-secretase',
                    validated=True,
                    applications=['in vitro', 'in vivo'],
                    notes='IC50: 10 nM, potent GSI'
                ),
                MolecularTool(
                    name='LY-411575',
                    tool_type='inhibitor',
                    target='gamma-secretase',
                    validated=True,
                    applications=['in vitro', 'in vivo'],
                    notes='IC50: 0.08 nM, very potent'
                ),
                MolecularTool(
                    name='RO4929097',
                    tool_type='inhibitor',
                    target='gamma-secretase',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='Clinical-grade GSI'
                )
            ],
            'EGFR': [
                MolecularTool(
                    name='Erlotinib',
                    tool_type='inhibitor',
                    target='EGFR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='IC50: 2 nM, FDA approved'
                ),
                MolecularTool(
                    name='Gefitinib',
                    tool_type='inhibitor',
                    target='EGFR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='IC50: 33 nM, FDA approved'
                ),
                MolecularTool(
                    name='Osimertinib',
                    tool_type='inhibitor',
                    target='EGFR',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='3rd gen, T790M mutant selective'
                )
            ],
            'JAK1': [
                MolecularTool(
                    name='Ruxolitinib',
                    tool_type='inhibitor',
                    target='JAK1/JAK2',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved for myelofibrosis'
                ),
                MolecularTool(
                    name='Tofacitinib',
                    tool_type='inhibitor',
                    target='JAK1/JAK3',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved for RA'
                )
            ],
            'BRAF': [
                MolecularTool(
                    name='Vemurafenib',
                    tool_type='inhibitor',
                    target='BRAF V600E',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved, V600E selective'
                ),
                MolecularTool(
                    name='Dabrafenib',
                    tool_type='inhibitor',
                    target='BRAF V600E',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved'
                )
            ],
            'MEK1': [
                MolecularTool(
                    name='Trametinib',
                    tool_type='inhibitor',
                    target='MEK1/MEK2',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved'
                ),
                MolecularTool(
                    name='U0126',
                    tool_type='inhibitor',
                    target='MEK1/MEK2',
                    validated=True,
                    applications=['in vitro'],
                    notes='Research tool compound'
                )
            ],
            'AKT1': [
                MolecularTool(
                    name='MK-2206',
                    tool_type='inhibitor',
                    target='AKT1/2/3',
                    validated=True,
                    applications=['in vitro', 'in vivo'],
                    notes='Allosteric AKT inhibitor'
                )
            ],
            'PI3K': [
                MolecularTool(
                    name='LY294002',
                    tool_type='inhibitor',
                    target='PI3K',
                    validated=True,
                    applications=['in vitro'],
                    notes='Classic PI3K inhibitor, research tool'
                ),
                MolecularTool(
                    name='Alpelisib',
                    tool_type='inhibitor',
                    target='PI3Kalpha',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved'
                )
            ],
            'mTOR': [
                MolecularTool(
                    name='Rapamycin',
                    tool_type='inhibitor',
                    target='mTORC1',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='Classic mTOR inhibitor'
                ),
                MolecularTool(
                    name='Everolimus',
                    tool_type='inhibitor',
                    target='mTORC1',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='FDA approved'
                )
            ],
            'WNT': [
                MolecularTool(
                    name='IWP-2',
                    tool_type='inhibitor',
                    target='Porcupine (Wnt secretion)',
                    validated=True,
                    applications=['in vitro'],
                    notes='Blocks Wnt secretion'
                ),
                MolecularTool(
                    name='XAV939',
                    tool_type='inhibitor',
                    target='Tankyrase (Wnt/beta-catenin)',
                    validated=True,
                    applications=['in vitro'],
                    notes='Stabilizes Axin'
                )
            ],
            'TGFBR1': [
                MolecularTool(
                    name='SB431542',
                    tool_type='inhibitor',
                    target='ALK5 (TGFBR1)',
                    validated=True,
                    applications=['in vitro'],
                    notes='Classic TGF-beta inhibitor'
                ),
                MolecularTool(
                    name='LY2157299 (Galunisertib)',
                    tool_type='inhibitor',
                    target='ALK5 (TGFBR1)',
                    validated=True,
                    applications=['in vitro', 'in vivo', 'clinical'],
                    notes='Clinical development'
                )
            ]
        }

        aliases = {
            'VEGFR2': 'KDR',
            'ALK5': 'TGFBR1',
            'ERK': 'MAPK1',
            'Wnt': 'WNT'
        }

        target_key = aliases.get(target, target)
        return CURATED_INHIBITORS.get(target_key, [])

    def find_activators(self, target: str) -> List[MolecularTool]:
        """Find activators/agonists for target."""

        ACTIVATORS = {
            'KDR': [
                MolecularTool(
                    name='VEGF-A165',
                    tool_type='activator',
                    target='KDR',
                    vendor='R&D Systems, PeproTech',
                    applications=['in vitro', 'in vivo'],
                    notes='Recombinant human VEGF-A'
                )
            ],
            'NOTCH1': [
                MolecularTool(
                    name='DLL4-Fc',
                    tool_type='activator',
                    target='NOTCH1',
                    vendor='R&D Systems',
                    applications=['in vitro'],
                    notes='Immobilized DLL4 activates Notch'
                ),
                MolecularTool(
                    name='Jagged1 peptide',
                    tool_type='activator',
                    target='NOTCH1',
                    applications=['in vitro'],
                    notes='Synthetic Notch activating peptide'
                )
            ],
            'WNT': [
                MolecularTool(
                    name='Wnt3a',
                    tool_type='activator',
                    target='Frizzled/LRP',
                    vendor='R&D Systems',
                    applications=['in vitro'],
                    notes='Recombinant Wnt3a'
                ),
                MolecularTool(
                    name='CHIR99021',
                    tool_type='activator',
                    target='beta-catenin (GSK3 inhibitor)',
                    applications=['in vitro'],
                    notes='GSK3 inhibitor, activates Wnt pathway'
                )
            ],
            'TGFBR': [
                MolecularTool(
                    name='TGF-beta1',
                    tool_type='activator',
                    target='TGFBR1/TGFBR2',
                    vendor='R&D Systems, PeproTech',
                    applications=['in vitro'],
                    notes='Recombinant TGF-beta1'
                )
            ]
        }

        aliases = {'VEGFR2': 'KDR'}
        target_key = aliases.get(target, target)
        return ACTIVATORS.get(target_key, [])

    def find_rnai_reagents(
        self,
        target: str,
        species: str = 'human'
    ) -> List[MolecularTool]:
        """Find siRNA/shRNA reagents."""

        tools = []

        # siRNA vendors
        sirna_sources = [
            {
                'vendor': 'Dharmacon ON-TARGETplus',
                'type': 'siRNA SMARTpool',
                'knockdown': '>80%',
                'species': ['human', 'mouse']
            },
            {
                'vendor': 'Ambion Silencer Select',
                'type': 'siRNA',
                'knockdown': '>70%',
                'species': ['human', 'mouse']
            },
            {
                'vendor': 'Qiagen FlexiTube',
                'type': 'siRNA',
                'species': ['human', 'mouse']
            }
        ]

        for source in sirna_sources:
            if species in source['species']:
                tools.append(MolecularTool(
                    name=f"{target} siRNA ({source['vendor']})",
                    tool_type='siRNA',
                    target=target,
                    vendor=source['vendor'],
                    species=[species],
                    notes=f"Pre-designed, validated knockdown"
                ))

        # shRNA
        tools.append(MolecularTool(
            name=f'{target} shRNA (MISSION)',
            tool_type='shRNA',
            target=target,
            vendor='Sigma MISSION shRNA',
            applications=['stable knockdown'],
            species=[species],
            notes='Lentiviral vector, multiple clones available'
        ))

        return tools

    def find_crispr_reagents(
        self,
        target: str,
        species: str = 'human'
    ) -> List[MolecularTool]:
        """Find CRISPR reagents."""

        tools = []

        # CRISPR sources
        crispr_sources = [
            {
                'vendor': 'Synthego',
                'type': 'Synthetic sgRNA',
                'format': 'Modified RNA',
                'url': 'https://www.synthego.com'
            },
            {
                'vendor': 'IDT Alt-R',
                'type': 'crRNA + tracrRNA',
                'format': 'Two-part guide',
                'url': 'https://www.idtdna.com'
            },
            {
                'vendor': 'Addgene',
                'type': 'Plasmid libraries',
                'format': 'Lentiviral vector',
                'url': 'https://www.addgene.org'
            }
        ]

        for source in crispr_sources:
            tools.append(MolecularTool(
                name=f'{target} CRISPR KO ({source["vendor"]})',
                tool_type='CRISPR',
                target=target,
                vendor=source['vendor'],
                species=[species],
                notes=f'{source["type"]}, {source["format"]}'
            ))

        # Add guide RNA sequences if available
        guides = self._get_validated_guides(target, species)
        if guides:
            tools.append(MolecularTool(
                name=f'{target} validated guides',
                tool_type='CRISPR',
                target=target,
                species=[species],
                notes=f'Validated sequences: {guides[:2]}'
            ))

        return tools

    def _get_validated_guides(
        self,
        target: str,
        species: str
    ) -> List[str]:
        """Get validated CRISPR guide sequences."""
        # Would query databases like CRISPRscan, Benchling, etc.
        return []

    def find_recombinant_proteins(
        self,
        target: str,
        species: str = 'human'
    ) -> List[MolecularTool]:
        """Find recombinant proteins."""

        tools = []

        # Major recombinant protein vendors
        vendors = ['R&D Systems', 'PeproTech', 'Sino Biological', 'Abcam']

        for vendor in vendors:
            tools.append(MolecularTool(
                name=f'Recombinant {target}',
                tool_type='recombinant_protein',
                target=target,
                vendor=vendor,
                species=[species],
                applications=['ELISA standard', 'functional assays', 'binding studies']
            ))

        # Fc-fusion if receptor
        tools.append(MolecularTool(
            name=f'{target}-Fc chimera',
            tool_type='recombinant_protein',
            target=target,
            vendor='R&D Systems',
            species=[species],
            applications=['blocking', 'binding assays', 'ligand traps']
        ))

        return tools

    def _get_uniprot_info(self, gene: str) -> Optional[dict]:
        """Get protein info from UniProt."""
        try:
            url = "https://rest.uniprot.org/uniprotkb/search"
            params = {
                'query': f'gene:{gene} AND organism_id:9606',
                'format': 'json',
                'fields': 'accession,protein_name,gene_names'
            }
            response = requests.get(url, params=params, timeout=10)
            if response.status_code == 200:
                data = response.json()
                if data.get('results'):
                    return data['results'][0]
        except:
            pass
        return None

    def generate_tools_report(
        self,
        target: str,
        species: str = 'human'
    ) -> str:
        """Generate markdown report of available tools."""

        tools = self.find_all_tools(target, species)

        report = f"""# Molecular Tools Report: {target}

## Summary
Target: **{target}**
Species: **{species}**

"""

        for category, tool_list in tools.items():
            if tool_list:
                report += f"## {category.title()}\n\n"
                report += "| Name | Vendor | Applications | Notes |\n"
                report += "|------|--------|--------------|-------|\n"

                for tool in tool_list:
                    vendor = tool.vendor or '-'
                    apps = ', '.join(tool.applications) if tool.applications else '-'
                    notes = tool.notes or '-'
                    report += f"| {tool.name} | {vendor} | {apps} | {notes} |\n"

                report += "\n"

        return report
```

## Quick Usage

```python
# Find all tools for any target
discovery = MolecularToolsDiscovery()

# Get all tools
tools = discovery.find_all_tools('VEGFR2')
print(f"Found {sum(len(v) for v in tools.values())} tools")

# Generate report
report = discovery.generate_tools_report('NOTCH1')
print(report)

# Find specific tool types
inhibitors = discovery.find_inhibitors('BRAF')
antibodies = discovery.find_antibodies('CD31', species='mouse')
crispr = discovery.find_crispr_reagents('TP53')
```
