# Pathway Deep Dive Reference

## Experimental Evidence Categories

### 1. Knockout/Knockdown Studies

#### Query Structure for KO Analysis
```python
def analyze_knockout_effects(gene, species='mouse', context=None):
    """
    Retrieve comprehensive knockout phenotype data.

    Returns:
    - Direct phenotypes
    - Compensatory mechanisms
    - Cell-type specific effects
    - Developmental timing effects
    """
    results = {
        'phenotypes': [],
        'compensatory': [],
        'cell_specific': [],
        'temporal': []
    }

    # Query databases
    # 1. Mouse Genome Informatics (MGI)
    mgi_data = query_mgi(gene)

    # 2. International Mouse Phenotyping Consortium (IMPC)
    impc_data = query_impc(gene)

    # 3. Literature (PubMed)
    lit_data = search_pubmed(f'{gene} knockout OR knockdown {context}')

    return results
```

#### Example: VEGFR2 Knockout Analysis
```yaml
Gene: KDR (VEGFR2)
Species: Mouse (Kdr)

Knockout Phenotypes:
  Embryonic:
    - Lethal at E8.5-E9.5
    - No blood island formation
    - Absent endothelial cells
    - No vascular development

  Conditional (Endothelial-specific):
    - Reduced tumor angiogenesis
    - Impaired wound healing
    - Decreased vascular permeability

  Heterozygous:
    - Viable, fertile
    - Reduced angiogenic response
    - Mild vascular defects

Compensatory Mechanisms:
  - VEGFR1 upregulation
  - VEGFR3 partial compensation in lymphatics
  - Alternative angiogenic pathways (FGF, PDGF)

Cell-Type Specific Effects:
  Tip Cells:
    - Complete loss of tip cell identity
    - No filopodia formation
    - Migration abolished

  Stalk Cells:
    - Reduced proliferation
    - Lumen formation defects

  Phalanx Cells:
    - Decreased vessel stability

References:
  - PMID: 9697695 (Shalaby et al., 1995)
  - PMID: 15314025 (Ferrara et al., 2004)
```

### 2. Treatment/Perturbation Studies

#### Chemical Perturbation Database
```python
treatment_database = {
    'DAPT': {
        'target': 'gamma-secretase',
        'mechanism': 'Notch cleavage inhibition',
        'pathway_effects': {
            'Notch': {
                'HES1': 'downregulated',
                'HEY1': 'downregulated',
                'HEY2': 'downregulated',
                'DLL4': 'unchanged (ligand)',
                'Jagged1': 'unchanged (ligand)'
            },
            'VEGF_crosstalk': {
                'tip_cell_formation': 'enhanced',
                'sprouting': 'increased',
                'VEGFR2': 'upregulated in stalk cells'
            }
        },
        'cell_type_responses': {
            'tip_cell': 'expansion',
            'stalk_cell': 'conversion to tip-like',
            'phalanx_cell': 'disrupted quiescence'
        },
        'dosage_effects': {
            '1uM': 'partial inhibition',
            '10uM': 'complete inhibition',
            '50uM': 'off-target effects'
        },
        'references': ['PMID:19033656', 'PMID:17855499']
    },

    'anti-VEGF (Bevacizumab)': {
        'target': 'VEGF-A',
        'mechanism': 'Ligand sequestration',
        'pathway_effects': {
            'VEGF': {
                'VEGFR2_activation': 'blocked',
                'downstream_ERK': 'reduced',
                'downstream_AKT': 'reduced'
            }
        },
        'resistance_mechanisms': [
            'VEGF-independent angiogenesis',
            'Vessel co-option',
            'Alternative growth factors (FGF, PDGF)',
            'Bone marrow-derived cells'
        ],
        'clinical_context': {
            'colorectal_cancer': 'first-line with chemo',
            'nsclc': 'with carboplatin/paclitaxel',
            'glioblastoma': 'second-line'
        }
    }
}
```

### 3. Cell-Type Specific Pathway Analysis

#### Endothelial Cell Types in Angiogenesis
```python
cell_type_pathways = {
    'tip_cell': {
        'markers': ['DLL4', 'VEGFR2', 'CXCR4', 'UNC5B', 'PDGFB'],
        'high_activity': ['VEGF signaling', 'Semaphorin signaling'],
        'low_activity': ['Notch signaling (receiving)'],
        'functions': [
            'Guide vessel sprouting',
            'Filopodia extension',
            'Matrix degradation',
            'Chemotaxis toward VEGF'
        ],
        'regulatory_circuits': {
            'positive': ['VEGF-A → VEGFR2 → DLL4'],
            'negative': ['DLL4 → Notch1 (in neighbors)']
        }
    },

    'stalk_cell': {
        'markers': ['Jagged1', 'VEGFR1', 'Notch1', 'HES1', 'HEY1'],
        'high_activity': ['Notch signaling', 'Proliferation'],
        'low_activity': ['VEGF-induced migration'],
        'functions': [
            'Proliferation',
            'Lumen formation',
            'Junctional stability',
            'Basement membrane deposition'
        ],
        'regulatory_circuits': {
            'positive': ['DLL4 → Notch1 → HES1/HEY1'],
            'negative': ['Notch → suppresses VEGFR2']
        }
    },

    'phalanx_cell': {
        'markers': ['VE-cadherin', 'CD31', 'Claudin-5'],
        'high_activity': ['Junctional signaling', 'Quiescence'],
        'low_activity': ['Proliferation', 'Migration'],
        'functions': [
            'Vessel stability',
            'Barrier function',
            'Perfusion'
        ]
    }
}
```

### 4. Crosstalk Analysis

#### VEGF-Notch Crosstalk
```python
crosstalk_map = {
    'VEGF_to_Notch': {
        'mechanism': 'VEGF induces DLL4 in tip cells',
        'pathway': 'VEGF-A → VEGFR2 → ERK → ETS → DLL4 transcription',
        'effect': 'Lateral inhibition of neighbors',
        'evidence': 'experimental',
        'references': ['PMID:17855499']
    },

    'Notch_to_VEGF': {
        'mechanism': 'Notch suppresses VEGFR2/VEGFR3',
        'pathway': 'Notch ICD → RBPJ → HEY1/HEY2 → VEGFR2 repression',
        'effect': 'Reduced VEGF responsiveness',
        'evidence': 'experimental',
        'references': ['PMID:19033656']
    },

    'balance_outcome': {
        'too_much_VEGF': 'Hypervascular, disorganized vessels',
        'too_much_Notch': 'Hypovascular, sparse vessels',
        'balanced': 'Organized, functional vasculature'
    }
}
```

#### Multi-pathway Crosstalk Network
```yaml
Crosstalk Map:
  VEGF ↔ Notch:
    Type: Mutual regulation
    Effect: Tip/stalk cell balance

  VEGF ↔ Wnt:
    Type: Synergistic
    Effect: Enhanced angiogenesis
    Mechanism: Wnt/β-catenin → VEGF transcription

  VEGF ↔ TGF-β:
    Type: Context-dependent
    Effect:
      - Low TGF-β: Pro-angiogenic
      - High TGF-β: Vessel stabilization

  VEGF ↔ HIF:
    Type: Upstream regulation
    Effect: Hypoxia → HIF-1α → VEGF transcription

  Notch ↔ Wnt:
    Type: Mutual inhibition
    Effect: Cell fate decisions

  VEGF ↔ Semaphorin:
    Type: Guidance cue integration
    Effect: Directional migration
```

### 5. Molecular Interaction Database

#### Ligand-Receptor Pairs
```python
ligand_receptor_db = {
    'VEGF_family': {
        'VEGF-A': {
            'receptors': ['VEGFR1', 'VEGFR2', 'NRP1', 'NRP2'],
            'affinity': {'VEGFR1': 'high (10pM)', 'VEGFR2': 'moderate (75pM)'},
            'biological_effect': {
                'VEGFR2': 'Angiogenesis, permeability',
                'VEGFR1': 'Decoy receptor, monocyte chemotaxis',
                'NRP1': 'Co-receptor, enhances VEGFR2 signaling'
            }
        },
        'VEGF-B': {
            'receptors': ['VEGFR1', 'NRP1'],
            'biological_effect': 'Fatty acid transport, neuronal survival'
        },
        'VEGF-C': {
            'receptors': ['VEGFR2', 'VEGFR3'],
            'biological_effect': 'Lymphangiogenesis, developmental angiogenesis'
        },
        'VEGF-D': {
            'receptors': ['VEGFR2', 'VEGFR3'],
            'biological_effect': 'Lymphangiogenesis'
        },
        'PlGF': {
            'receptors': ['VEGFR1', 'NRP1', 'NRP2'],
            'biological_effect': 'Pathological angiogenesis, monocyte recruitment'
        }
    },

    'Notch_family': {
        'DLL1': {'receptors': ['Notch1', 'Notch2', 'Notch3', 'Notch4']},
        'DLL3': {'receptors': ['Notch1']},
        'DLL4': {
            'receptors': ['Notch1', 'Notch4'],
            'expression': 'Arterial endothelium, tip cells',
            'biological_effect': 'Tip/stalk selection, arterial identity'
        },
        'Jagged1': {
            'receptors': ['Notch1', 'Notch2'],
            'expression': 'Stalk cells, smooth muscle',
            'biological_effect': 'Pro-angiogenic in some contexts'
        },
        'Jagged2': {'receptors': ['Notch1', 'Notch2']}
    }
}
```

#### Kinase-Substrate Relationships
```python
kinase_substrate_db = {
    'VEGFR2': {
        'substrates': {
            'PLCgamma1': {
                'site': 'Y1175',
                'effect': 'PKC activation → ERK pathway',
                'phenotype': 'Proliferation, permeability'
            },
            'Src': {
                'site': 'Y951',
                'effect': 'Adhesion, migration',
                'phenotype': 'Vascular permeability'
            },
            'PI3K': {
                'site': 'Y1175 (indirect via adapters)',
                'effect': 'AKT activation',
                'phenotype': 'Survival, migration'
            },
            'Shb': {
                'site': 'Y1175',
                'effect': 'PI3K/AKT pathway',
                'phenotype': 'Angiogenesis'
            }
        },
        'autophosphorylation': ['Y951', 'Y1054', 'Y1059', 'Y1175', 'Y1214']
    }
}
```

### 6. Validation Reagents Database

#### Antibodies
```python
antibody_db = {
    'VEGFR2': {
        'blocking': [
            {
                'name': 'DC101',
                'species': 'mouse',
                'type': 'rat anti-mouse',
                'application': 'in vivo blocking',
                'vendor': 'Bio X Cell',
                'reference': 'PMID: 9175789'
            },
            {
                'name': 'Ramucirumab (Cyramza)',
                'species': 'human',
                'type': 'human IgG1',
                'application': 'clinical therapeutic',
                'FDA_approved': True,
                'indications': ['gastric cancer', 'NSCLC', 'CRC', 'HCC']
            }
        ],
        'detection': [
            {
                'name': 'Anti-VEGFR2 (55B11)',
                'application': ['WB', 'IP', 'IHC'],
                'vendor': 'Cell Signaling',
                'catalog': '#2479'
            },
            {
                'name': 'Anti-pVEGFR2 (Y1175)',
                'application': ['WB', 'IHC'],
                'vendor': 'Cell Signaling',
                'catalog': '#2478'
            }
        ]
    },

    'DLL4': {
        'blocking': [
            {
                'name': 'anti-DLL4 (REGN421)',
                'species': 'human',
                'type': 'human IgG1',
                'effect': 'Blocks Notch activation',
                'phenotype': 'Hypervascular, non-productive angiogenesis',
                'clinical_status': 'Phase 1 completed'
            }
        ]
    }
}
```

#### Small Molecule Inhibitors
```python
inhibitor_db = {
    'VEGFR2': {
        'approved': [
            {
                'name': 'Sunitinib',
                'targets': ['VEGFR1', 'VEGFR2', 'VEGFR3', 'PDGFR', 'KIT', 'FLT3'],
                'IC50_VEGFR2': '80 nM',
                'indications': ['RCC', 'GIST', 'pNET']
            },
            {
                'name': 'Sorafenib',
                'targets': ['VEGFR2', 'VEGFR3', 'PDGFR', 'RAF'],
                'IC50_VEGFR2': '90 nM',
                'indications': ['HCC', 'RCC', 'DTC']
            },
            {
                'name': 'Axitinib',
                'targets': ['VEGFR1', 'VEGFR2', 'VEGFR3'],
                'IC50_VEGFR2': '0.2 nM',
                'indications': ['RCC'],
                'note': 'Most selective VEGFR inhibitor'
            }
        ],
        'investigational': [
            {
                'name': 'Tivozanib',
                'IC50_VEGFR2': '0.16 nM',
                'status': 'Phase 3'
            }
        ]
    },

    'gamma_secretase': {
        'tool_compounds': [
            {
                'name': 'DAPT',
                'IC50': '20 nM',
                'use': 'In vitro/in vivo research',
                'vendor': 'Sigma, Tocris'
            },
            {
                'name': 'DBZ',
                'IC50': '10 nM',
                'use': 'In vivo research'
            }
        ]
    }
}
```

#### siRNA/shRNA
```python
rnai_db = {
    'VEGFR2': {
        'validated_sirna': [
            {
                'vendor': 'Dharmacon ON-TARGETplus',
                'catalog': 'L-003148-00-0005',
                'species': 'human',
                'knockdown_efficiency': '>80%'
            },
            {
                'vendor': 'Ambion Silencer Select',
                'catalog': 's7821, s7822, s7823',
                'species': 'human'
            }
        ],
        'shrna': [
            {
                'source': 'MISSION shRNA',
                'vector': 'pLKO.1',
                'clones': 5,
                'species': 'human, mouse'
            }
        ]
    }
}
```

#### CRISPR Resources
```python
crispr_db = {
    'VEGFR2': {
        'human_guides': [
            {'sequence': 'GCCTACCTCACCTGTTTCCT', 'efficiency': 0.82},
            {'sequence': 'AGATCACTACCATCCACCGG', 'efficiency': 0.78}
        ],
        'mouse_guides': [
            {'sequence': 'GCCTACCTCACGTGTTTCCT', 'efficiency': 0.85}
        ],
        'resources': {
            'Addgene': 'Various KO plasmids',
            'Synthego': 'Synthetic sgRNA',
            'IDT': 'Alt-R CRISPR system'
        }
    }
}
```

## Database Query Templates

### PubMed Query for Experimental Evidence
```
("{gene}" AND ("knockout" OR "knockdown" OR "KO" OR "siRNA" OR "shRNA" OR "CRISPR"))
AND ("{pathway}" OR "{process}")
AND ("{cell_type}" OR "{tissue}")
AND ("mouse" OR "mice" OR "human" OR "in vivo" OR "in vitro")
```

### STRING Query for Protein Interactions
```python
def query_string_interactions(gene_list, species=9606, score_threshold=700):
    """Query STRING for protein interactions with high confidence."""
    import requests

    string_api = "https://string-db.org/api"

    # Get interaction network
    params = {
        'identifiers': '%0d'.join(gene_list),
        'species': species,
        'required_score': score_threshold,
        'network_type': 'functional'
    }

    response = requests.get(f"{string_api}/tsv/network", params=params)
    return response.text
```

### Reactome Query for Pathway Details
```python
def query_reactome_pathway(pathway_id):
    """Get detailed pathway information from Reactome."""
    import requests

    base_url = "https://reactome.org/ContentService"

    # Get pathway details
    pathway = requests.get(f"{base_url}/data/query/{pathway_id}").json()

    # Get participating molecules
    molecules = requests.get(
        f"{base_url}/data/pathway/{pathway_id}/participatingPhysicalEntities"
    ).json()

    # Get reactions
    reactions = requests.get(
        f"{base_url}/data/pathway/{pathway_id}/containedEvents"
    ).json()

    return {
        'pathway': pathway,
        'molecules': molecules,
        'reactions': reactions
    }
```
