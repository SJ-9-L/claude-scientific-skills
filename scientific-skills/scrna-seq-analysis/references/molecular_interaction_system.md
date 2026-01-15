# Comprehensive Molecular Interaction System

## Overview

어떤 분자(gene, ligand, receptor, TF, chemical, protein, antibody, cytokine, chemokine)를 입력하면
모든 interaction 정보를 반환하고, 데이터 축적을 통해 KO phenotype을 예측하는 시스템

---

## 1. Molecule Type Definitions

```python
from enum import Enum
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Set, Tuple
import json

class MoleculeType(Enum):
    """All supported molecule types"""
    GENE = "gene"
    PROTEIN = "protein"
    LIGAND = "ligand"
    RECEPTOR = "receptor"
    TRANSCRIPTION_FACTOR = "transcription_factor"
    KINASE = "kinase"
    PHOSPHATASE = "phosphatase"
    CHEMICAL = "chemical"
    ANTIBODY = "antibody"
    CYTOKINE = "cytokine"
    CHEMOKINE = "chemokine"
    GROWTH_FACTOR = "growth_factor"
    ENZYME = "enzyme"
    ION_CHANNEL = "ion_channel"
    GPCR = "gpcr"
    NUCLEAR_RECEPTOR = "nuclear_receptor"
    ADAPTOR = "adaptor_protein"
    SCAFFOLD = "scaffold_protein"
    MIRNA = "mirna"
    LNCRNA = "lncrna"


class InteractionType(Enum):
    """All interaction types"""
    # Physical interactions
    BINDING = "binding"
    COMPLEX_FORMATION = "complex_formation"
    DIMERIZATION = "dimerization"
    OLIGOMERIZATION = "oligomerization"

    # Functional regulation
    ACTIVATION = "activation"
    INHIBITION = "inhibition"
    POSITIVE_REGULATION = "positive_regulation"
    NEGATIVE_REGULATION = "negative_regulation"
    MODULATION = "modulation"

    # Enzymatic modifications
    PHOSPHORYLATION = "phosphorylation"
    DEPHOSPHORYLATION = "dephosphorylation"
    UBIQUITINATION = "ubiquitination"
    SUMOYLATION = "sumoylation"
    ACETYLATION = "acetylation"
    METHYLATION = "methylation"
    GLYCOSYLATION = "glycosylation"
    PROTEOLYTIC_CLEAVAGE = "proteolytic_cleavage"

    # Transcriptional
    TRANSCRIPTION_ACTIVATION = "transcription_activation"
    TRANSCRIPTION_REPRESSION = "transcription_repression"
    CHROMATIN_REMODELING = "chromatin_remodeling"

    # Signaling
    SIGNAL_TRANSDUCTION = "signal_transduction"
    CROSSTALK = "crosstalk"
    FEEDBACK_POSITIVE = "feedback_positive"
    FEEDBACK_NEGATIVE = "feedback_negative"
    FEEDFORWARD = "feedforward"

    # Localization
    NUCLEAR_TRANSLOCATION = "nuclear_translocation"
    MEMBRANE_RECRUITMENT = "membrane_recruitment"
    SECRETION = "secretion"
    ENDOCYTOSIS = "endocytosis"

    # Expression
    UPREGULATION = "upregulation"
    DOWNREGULATION = "downregulation"
    STABILIZATION = "stabilization"
    DEGRADATION = "degradation"


class EvidenceType(Enum):
    """Evidence quality levels"""
    EXPERIMENTAL_DIRECT = "experimental_direct"      # Direct biochemical assay
    EXPERIMENTAL_INDIRECT = "experimental_indirect"  # Genetic/functional evidence
    HIGH_THROUGHPUT = "high_throughput"              # Proteomics, Y2H, etc.
    LITERATURE_CURATED = "literature_curated"        # Manual curation
    DATABASE_INFERRED = "database_inferred"          # Computational prediction
    TEXT_MINING = "text_mining"                      # NLP extracted
    PREDICTED = "predicted"                          # In silico prediction


@dataclass
class Molecule:
    """Universal molecule representation"""
    symbol: str                           # Primary symbol (e.g., "VEGFA", "IL6")
    name: str                             # Full name
    molecule_type: MoleculeType           # Type classification
    aliases: List[str] = field(default_factory=list)  # Alternative names

    # Identifiers
    uniprot_id: Optional[str] = None
    gene_id: Optional[str] = None         # Entrez/NCBI Gene ID
    ensembl_id: Optional[str] = None
    hgnc_id: Optional[str] = None
    chembl_id: Optional[str] = None       # For chemicals
    drugbank_id: Optional[str] = None

    # Classification
    family: Optional[str] = None          # e.g., "RTK", "GPCR", "Interleukin"
    subfamily: Optional[str] = None

    # Properties
    molecular_weight: Optional[float] = None
    subcellular_location: List[str] = field(default_factory=list)
    tissue_expression: Dict[str, float] = field(default_factory=dict)

    # Functional annotation
    function: Optional[str] = None
    biological_process: List[str] = field(default_factory=list)
    pathways: List[str] = field(default_factory=list)

    # Disease associations
    disease_associations: List[Dict] = field(default_factory=list)


@dataclass
class Interaction:
    """Single molecular interaction"""
    source: Molecule                      # Source molecule
    target: Molecule                      # Target molecule
    interaction_type: InteractionType     # Type of interaction

    # Direction and effect
    direction: str = "directed"           # "directed", "undirected", "bidirectional"
    effect: str = "unknown"               # "positive", "negative", "neutral", "context_dependent"

    # Quantitative data
    affinity_kd: Optional[float] = None   # Binding affinity (nM)
    ec50: Optional[float] = None          # Activation EC50
    ic50: Optional[float] = None          # Inhibition IC50
    fold_change: Optional[float] = None   # Expression change

    # Context
    cell_type_specific: List[str] = field(default_factory=list)
    tissue_specific: List[str] = field(default_factory=list)
    disease_context: List[str] = field(default_factory=list)
    species: List[str] = field(default_factory=lambda: ["human", "mouse"])

    # Evidence
    evidence_type: EvidenceType = EvidenceType.DATABASE_INFERRED
    evidence_score: float = 0.5           # 0-1 confidence score
    pubmed_ids: List[str] = field(default_factory=list)
    source_database: List[str] = field(default_factory=list)

    # Mechanism details
    mechanism: Optional[str] = None       # Detailed mechanism description
    residues_involved: List[str] = field(default_factory=list)  # e.g., ["Y1175", "S473"]

    def to_dict(self) -> dict:
        return {
            "source": self.source.symbol,
            "target": self.target.symbol,
            "interaction_type": self.interaction_type.value,
            "effect": self.effect,
            "evidence_score": self.evidence_score,
            "mechanism": self.mechanism,
            "pubmed_ids": self.pubmed_ids[:3],  # Top 3 references
        }
```

---

## 2. Comprehensive Interaction Database

```python
class MolecularInteractionDB:
    """
    Central database for all molecular interactions.
    Aggregates data from multiple sources and enables KO prediction.
    """

    # Major signaling pathway templates
    PATHWAY_TEMPLATES = {
        "RTK_signaling": {
            "canonical_members": {
                "ligands": ["EGF", "VEGF", "FGF", "PDGF", "IGF1", "HGF", "NGF"],
                "receptors": ["EGFR", "VEGFR2", "FGFR1", "PDGFRA", "IGF1R", "MET", "NTRK1"],
                "adaptors": ["GRB2", "SHC1", "GAB1", "IRS1"],
                "kinases": ["SRC", "PLCγ", "PI3K", "RAF", "MEK", "ERK", "AKT", "JNK", "p38"],
                "phosphatases": ["PTPN11", "PTEN", "PP2A", "DUSP1"],
                "transcription_factors": ["ETS1", "AP1", "STAT3", "FOXO1", "MYC"],
            },
            "core_cascade": [
                ("ligand", "receptor", "binding"),
                ("receptor", "receptor", "dimerization"),
                ("receptor", "adaptor", "phosphorylation"),
                ("adaptor", "RAS", "activation"),
                ("RAS", "RAF", "activation"),
                ("RAF", "MEK", "phosphorylation"),
                ("MEK", "ERK", "phosphorylation"),
                ("ERK", "TF", "phosphorylation"),
            ]
        },

        "Notch_signaling": {
            "canonical_members": {
                "ligands": ["DLL1", "DLL3", "DLL4", "JAG1", "JAG2"],
                "receptors": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"],
                "proteases": ["ADAM10", "ADAM17", "γ-secretase", "PSEN1"],
                "transcription_factors": ["RBPJ", "MAML1"],
                "targets": ["HES1", "HES5", "HEY1", "HEY2", "HEYL", "MYC", "CCND1"],
            },
            "core_cascade": [
                ("DLL4", "NOTCH1", "binding"),
                ("ADAM10", "NOTCH1", "proteolytic_cleavage"),
                ("γ-secretase", "NOTCH1", "proteolytic_cleavage"),
                ("NICD", "RBPJ", "binding"),
                ("NICD-RBPJ", "MAML1", "complex_formation"),
                ("complex", "HES1", "transcription_activation"),
            ]
        },

        "Wnt_signaling": {
            "canonical_members": {
                "ligands": ["WNT1", "WNT3A", "WNT5A", "WNT7A", "WNT11"],
                "receptors": ["FZD1", "FZD2", "FZD7", "LRP5", "LRP6", "ROR2", "RYK"],
                "cytoplasmic": ["DVL1", "AXIN1", "APC", "GSK3B", "CK1", "βTrCP"],
                "transcription_factors": ["CTNNB1", "TCF7", "LEF1", "TCF7L2"],
                "targets": ["MYC", "CCND1", "AXIN2", "CD44", "MMP7"],
            }
        },

        "TGFβ_signaling": {
            "canonical_members": {
                "ligands": ["TGFB1", "TGFB2", "TGFB3", "BMP2", "BMP4", "BMP7", "ACTIVIN", "NODAL", "GDF"],
                "receptors": ["TGFBR1", "TGFBR2", "BMPR1A", "BMPR2", "ACVR1", "ACVR2A"],
                "smads": ["SMAD1", "SMAD2", "SMAD3", "SMAD4", "SMAD5", "SMAD6", "SMAD7"],
                "transcription_factors": ["SMAD4", "RUNX", "ATF3", "JUNB"],
                "targets": ["SERPINE1", "COL1A1", "CDKN1A", "CDKN2B", "ID1", "ID2"],
            }
        },

        "NFκB_signaling": {
            "canonical_members": {
                "ligands": ["TNF", "IL1B", "LPS", "CD40L"],
                "receptors": ["TNFR1", "IL1R1", "TLR4", "CD40"],
                "kinases": ["RIPK1", "TRAF2", "TRAF6", "TAK1", "IKKα", "IKKβ", "IKKγ"],
                "inhibitors": ["NFKBIA", "NFKBIB", "A20", "CYLD"],
                "transcription_factors": ["NFKB1", "RELA", "RELB", "REL", "NFKB2"],
                "targets": ["IL6", "IL8", "TNF", "ICAM1", "VCAM1", "COX2", "iNOS", "BCL2", "BIRC3"],
            }
        },

        "JAK_STAT_signaling": {
            "canonical_members": {
                "ligands": ["IFNα", "IFNβ", "IFNγ", "IL2", "IL4", "IL6", "IL7", "IL12", "IL15", "IL21", "EPO", "GH", "PRL"],
                "receptors": ["IFNAR", "IFNGR", "IL2R", "IL4R", "IL6R", "gp130"],
                "kinases": ["JAK1", "JAK2", "JAK3", "TYK2"],
                "stats": ["STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6"],
                "inhibitors": ["SOCS1", "SOCS2", "SOCS3", "CIS", "PIAS1", "PIAS3"],
                "targets": ["IRF1", "IRF9", "MYC", "BCL2", "CCND1", "MCL1", "VEGFA"],
            }
        },

        "PI3K_AKT_mTOR": {
            "canonical_members": {
                "activators": ["RTKs", "GPCRs", "integrins", "RAS"],
                "kinases": ["PI3K", "PDK1", "AKT1", "AKT2", "AKT3", "mTORC1", "mTORC2", "S6K1", "4EBP1"],
                "phosphatases": ["PTEN", "INPP4B", "PHLPP"],
                "substrates": ["GSK3B", "FOXO1", "FOXO3", "BAD", "MDM2", "TSC2"],
                "targets": ["MYC", "CCND1", "HIF1A", "SREBP", "ribosome_biogenesis"],
            }
        },

        "Hippo_signaling": {
            "canonical_members": {
                "upstream": ["NF2", "KIBRA", "TAOK", "cell_polarity", "mechanotransduction"],
                "core_kinases": ["MST1", "MST2", "SAV1", "LATS1", "LATS2", "MOB1A", "MOB1B"],
                "effectors": ["YAP1", "TAZ", "TEAD1", "TEAD2", "TEAD3", "TEAD4"],
                "targets": ["CTGF", "CYR61", "ANKRD1", "BIRC5", "MYC"],
            }
        },
    }

    # Cytokine/Chemokine interaction networks
    CYTOKINE_NETWORKS = {
        "pro_inflammatory": {
            "IL1_family": {
                "members": ["IL1A", "IL1B", "IL1RA", "IL18", "IL33", "IL36", "IL37", "IL38"],
                "receptors": ["IL1R1", "IL1R2", "IL1RAP", "IL18R1", "ST2"],
                "signaling": ["MyD88", "IRAK1", "IRAK4", "TRAF6", "NFκB", "MAPK"],
                "induced_genes": ["IL6", "IL8", "TNF", "COX2", "iNOS", "MMP1", "MMP3"],
            },
            "IL6_family": {
                "members": ["IL6", "IL11", "IL27", "IL31", "LIF", "OSM", "CNTF", "CT1"],
                "receptors": ["IL6R", "IL11R", "gp130", "LIFR", "OSMR"],
                "signaling": ["JAK1", "JAK2", "TYK2", "STAT3", "STAT1", "MAPK", "PI3K"],
                "induced_genes": ["SOCS3", "CRP", "SAA", "hepcidin", "fibrinogen"],
            },
            "TNF_family": {
                "members": ["TNF", "LTA", "LTB", "FASL", "TRAIL", "RANKL", "CD40L", "BAFF", "APRIL"],
                "receptors": ["TNFR1", "TNFR2", "FAS", "TRAILR1", "TRAILR2", "RANK", "CD40", "BAFFR"],
                "signaling": ["TRADD", "TRAF2", "RIPK1", "cIAP", "NFκB", "MAPK", "caspase-8"],
            },
            "IFN_family": {
                "type_I": ["IFNA", "IFNB", "IFNE", "IFNK", "IFNW"],
                "type_II": ["IFNG"],
                "type_III": ["IFNL1", "IFNL2", "IFNL3"],
                "receptors": ["IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNLR1"],
                "signaling": ["JAK1", "TYK2", "JAK2", "STAT1", "STAT2", "IRF9"],
                "ISGs": ["MX1", "OAS1", "ISG15", "IFIT1", "IFI44", "CXCL10"],
            },
        },

        "anti_inflammatory": {
            "IL10_family": {
                "members": ["IL10", "IL19", "IL20", "IL22", "IL24", "IL26", "IL28", "IL29"],
                "receptors": ["IL10R1", "IL10R2", "IL20R1", "IL20R2", "IL22R1"],
                "signaling": ["JAK1", "TYK2", "STAT3", "STAT1"],
                "effects": ["suppress_Th1", "suppress_macrophage", "promote_Treg"],
            },
            "TGFβ_family": {
                "members": ["TGFB1", "TGFB2", "TGFB3"],
                "receptors": ["TGFBR1", "TGFBR2", "TGFBR3"],
                "signaling": ["SMAD2", "SMAD3", "SMAD4", "SMAD7"],
                "effects": ["suppress_proliferation", "induce_Treg", "promote_fibrosis"],
            },
        },

        "chemokines": {
            "CC_family": {
                "members": ["CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL11", "CCL17", "CCL19", "CCL20", "CCL21", "CCL22", "CCL25"],
                "receptors": ["CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR9", "CCR10"],
                "functions": {
                    "CCL2-CCR2": "monocyte_recruitment",
                    "CCL5-CCR5": "T_cell_recruitment",
                    "CCL19/21-CCR7": "lymph_node_homing",
                    "CCL17/22-CCR4": "Th2_recruitment",
                }
            },
            "CXC_family": {
                "members": ["CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8", "CXCL9", "CXCL10", "CXCL11", "CXCL12", "CXCL13"],
                "receptors": ["CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6"],
                "functions": {
                    "CXCL8-CXCR1/2": "neutrophil_recruitment",
                    "CXCL9/10/11-CXCR3": "Th1_recruitment",
                    "CXCL12-CXCR4": "stem_cell_homing",
                    "CXCL13-CXCR5": "B_cell_recruitment",
                }
            },
        }
    }

    # Ligand-Receptor comprehensive database
    LIGAND_RECEPTOR_PAIRS = {
        # Growth factors
        "VEGF_family": {
            "VEGFA": ["KDR", "FLT1", "NRP1", "NRP2"],
            "VEGFB": ["FLT1", "NRP1"],
            "VEGFC": ["FLT4", "KDR", "NRP2"],
            "VEGFD": ["FLT4", "KDR"],
            "PlGF": ["FLT1", "NRP1", "NRP2"],
        },
        "FGF_family": {
            "FGF1": ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
            "FGF2": ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
            "FGF7": ["FGFR2b"],
            "FGF10": ["FGFR2b"],
            "FGF21": ["FGFR1c", "KLB"],
            "FGF23": ["FGFR1c", "KL"],
        },
        "EGF_family": {
            "EGF": ["EGFR"],
            "TGFA": ["EGFR"],
            "HB-EGF": ["EGFR", "ERBB4"],
            "AREG": ["EGFR"],
            "EREG": ["EGFR", "ERBB4"],
            "BTC": ["EGFR", "ERBB4"],
            "NRG1": ["ERBB3", "ERBB4"],
            "NRG2": ["ERBB3", "ERBB4"],
        },
        "PDGF_family": {
            "PDGFA": ["PDGFRA"],
            "PDGFB": ["PDGFRA", "PDGFRB"],
            "PDGFC": ["PDGFRA"],
            "PDGFD": ["PDGFRB"],
        },
        "Angiopoietin_family": {
            "ANGPT1": ["TEK"],
            "ANGPT2": ["TEK", "ITGAVB5"],
            "ANGPT4": ["TEK"],
        },

        # Notch ligands
        "Notch_ligands": {
            "DLL1": ["NOTCH1", "NOTCH2"],
            "DLL3": ["NOTCH1"],
            "DLL4": ["NOTCH1", "NOTCH4"],
            "JAG1": ["NOTCH1", "NOTCH2", "NOTCH3"],
            "JAG2": ["NOTCH1", "NOTCH2"],
        },

        # Wnt ligands
        "Wnt_ligands": {
            "WNT1": ["FZD1", "FZD4", "LRP5", "LRP6"],
            "WNT3A": ["FZD1", "FZD2", "FZD7", "LRP5", "LRP6"],
            "WNT5A": ["FZD2", "FZD5", "ROR2", "RYK"],
            "WNT7A": ["FZD4", "FZD5", "LRP5", "LRP6"],
            "WNT11": ["FZD7", "ROR2"],
        },

        # Hedgehog ligands
        "Hedgehog_ligands": {
            "SHH": ["PTCH1", "PTCH2"],
            "IHH": ["PTCH1"],
            "DHH": ["PTCH1"],
        },

        # Interleukins (selected)
        "Interleukins": {
            "IL1B": ["IL1R1", "IL1RAP"],
            "IL2": ["IL2RA", "IL2RB", "IL2RG"],
            "IL4": ["IL4R", "IL13RA1", "IL2RG"],
            "IL6": ["IL6R", "IL6ST"],
            "IL7": ["IL7R", "IL2RG"],
            "IL10": ["IL10RA", "IL10RB"],
            "IL12": ["IL12RB1", "IL12RB2"],
            "IL15": ["IL15RA", "IL2RB", "IL2RG"],
            "IL17A": ["IL17RA", "IL17RC"],
            "IL21": ["IL21R", "IL2RG"],
            "IL23": ["IL23R", "IL12RB1"],
        },
    }

    def __init__(self):
        self.molecules: Dict[str, Molecule] = {}
        self.interactions: List[Interaction] = []
        self.ko_phenotypes: Dict[str, List[Dict]] = {}

    def query_molecule(self, symbol: str) -> Dict:
        """
        Universal query: returns ALL known information about any molecule
        """
        result = {
            "molecule": None,
            "upstream_regulators": [],      # What regulates this molecule
            "downstream_targets": [],        # What this molecule regulates
            "binding_partners": [],          # Physical interactions
            "pathway_membership": [],        # Which pathways it belongs to
            "crosstalk": [],                 # Cross-pathway interactions
            "cell_type_expression": {},      # Where it's expressed
            "disease_associations": [],      # Disease relevance
            "ko_phenotypes": [],             # Known knockout phenotypes
            "predicted_ko_effects": [],      # Predicted effects
            "molecular_tools": {},           # Available reagents
        }

        # Implementation in next sections
        return result
```

---

## 3. Interaction Query Engine

```python
class InteractionQueryEngine:
    """
    Query engine for comprehensive interaction lookup.
    Returns all interactions for any molecule input.
    """

    def __init__(self, db: MolecularInteractionDB):
        self.db = db

    def get_all_interactions(self, molecule: str,
                             include_indirect: bool = True,
                             max_depth: int = 2) -> Dict:
        """
        Get ALL interactions for a molecule.

        Returns:
        - direct_interactions: Immediate partners
        - indirect_interactions: 2nd/3rd degree connections
        - regulatory_network: Full upstream/downstream network
        """

        result = {
            "query": molecule,
            "molecule_info": self._get_molecule_info(molecule),

            # UPSTREAM: What regulates this molecule
            "upstream": {
                "activators": self._get_activators(molecule),
                "inhibitors": self._get_inhibitors(molecule),
                "transcriptional_regulators": self._get_transcriptional_regulators(molecule),
                "post_translational_modifiers": self._get_ptm_regulators(molecule),
            },

            # DOWNSTREAM: What this molecule regulates
            "downstream": {
                "activates": self._get_activation_targets(molecule),
                "inhibits": self._get_inhibition_targets(molecule),
                "transcriptional_targets": self._get_transcriptional_targets(molecule),
                "substrates": self._get_substrates(molecule),
            },

            # PHYSICAL INTERACTIONS
            "physical": {
                "binding_partners": self._get_binding_partners(molecule),
                "complex_members": self._get_complex_members(molecule),
                "receptor_ligand": self._get_receptor_ligand_pairs(molecule),
            },

            # PATHWAY CONTEXT
            "pathway_context": {
                "pathways": self._get_pathway_membership(molecule),
                "crosstalk": self._get_crosstalk_points(molecule),
                "feedback_loops": self._get_feedback_loops(molecule),
            },

            # CELL TYPE SPECIFICITY
            "cell_type_specific": self._get_cell_type_specificity(molecule),

            # DISEASE CONTEXT
            "disease_context": self._get_disease_context(molecule),
        }

        # Add indirect interactions if requested
        if include_indirect:
            result["indirect_network"] = self._expand_network(molecule, max_depth)

        return result

    def _get_activators(self, molecule: str) -> List[Dict]:
        """Get all molecules that activate/upregulate the target"""
        # Query multiple databases
        activators = []

        # SIGNOR database
        signor_data = self._query_signor(molecule, "activation")
        activators.extend(signor_data)

        # KEGG pathways
        kegg_data = self._query_kegg(molecule, "activation")
        activators.extend(kegg_data)

        # Reactome
        reactome_data = self._query_reactome(molecule, "positive_regulation")
        activators.extend(reactome_data)

        # STRING (high confidence only)
        string_data = self._query_string(molecule, score_threshold=0.7)
        activators.extend([s for s in string_data if s.get("mode") == "activation"])

        return self._deduplicate_and_rank(activators)

    def _get_inhibitors(self, molecule: str) -> List[Dict]:
        """Get all molecules that inhibit/downregulate the target"""
        inhibitors = []

        # Same multi-database approach
        # Include: phosphatases for kinases, negative regulators, antagonists

        return inhibitors

    def _get_crosstalk_points(self, molecule: str) -> List[Dict]:
        """
        Identify pathway crosstalk points.
        Critical for understanding complex signaling.
        """
        crosstalk = []

        # Get pathways this molecule belongs to
        pathways = self._get_pathway_membership(molecule)

        # For each pathway, find connections to other pathways
        for pathway in pathways:
            other_pathways = self._find_pathway_connections(pathway, exclude=molecule)
            for other in other_pathways:
                crosstalk.append({
                    "from_pathway": pathway["name"],
                    "to_pathway": other["name"],
                    "connection_type": other["connection_type"],
                    "mediators": other["mediators"],
                    "effect": other["effect"],  # synergistic, antagonistic, context-dependent
                    "evidence": other["evidence"],
                })

        return crosstalk

    def _expand_network(self, molecule: str, max_depth: int) -> Dict:
        """
        Expand interaction network to specified depth.
        Useful for understanding broader signaling context.
        """
        network = {
            "nodes": [],
            "edges": [],
            "depth_layers": {}
        }

        visited = set()
        current_layer = {molecule}

        for depth in range(max_depth):
            network["depth_layers"][depth] = list(current_layer)
            next_layer = set()

            for mol in current_layer:
                if mol not in visited:
                    visited.add(mol)

                    # Get direct interactions
                    interactions = self._get_direct_interactions(mol)

                    for interaction in interactions:
                        network["edges"].append({
                            "source": interaction["source"],
                            "target": interaction["target"],
                            "type": interaction["type"],
                            "depth": depth
                        })

                        # Add to next layer
                        partner = interaction["target"] if interaction["source"] == mol else interaction["source"]
                        if partner not in visited:
                            next_layer.add(partner)

            current_layer = next_layer

        network["nodes"] = list(visited)
        return network
```

---

## 4. KO Phenotype Prediction Engine

```python
class KOPhenotypePredictor:
    """
    Predicts knockout phenotypes based on accumulated interaction data.
    Uses network topology and known phenotypes to make predictions.
    """

    def __init__(self, interaction_db: MolecularInteractionDB):
        self.db = interaction_db
        self.known_phenotypes = self._load_known_phenotypes()

    def predict_ko_phenotype(self, gene: str,
                             species: str = "mouse",
                             context: str = None) -> Dict:
        """
        Predict phenotypes for a gene knockout.

        Uses:
        1. Known phenotypes for this gene
        2. Known phenotypes for pathway members
        3. Network position analysis
        4. Functional annotation
        """

        prediction = {
            "gene": gene,
            "species": species,
            "context": context,

            # Known data
            "known_phenotypes": self._get_known_phenotypes(gene, species),

            # Predicted based on network
            "predicted_phenotypes": [],
            "confidence_scores": {},

            # Pathway impact
            "affected_pathways": [],
            "downstream_cascade": [],

            # Compensatory mechanisms
            "potential_compensation": [],

            # Experimental recommendations
            "validation_experiments": [],
        }

        # 1. Get network position
        network_position = self._analyze_network_position(gene)

        # 2. Predict based on pathway membership
        for pathway in network_position["pathways"]:
            pathway_phenotypes = self._predict_from_pathway(gene, pathway)
            prediction["predicted_phenotypes"].extend(pathway_phenotypes)

        # 3. Predict based on similar genes
        similar_genes = self._find_similar_genes(gene)
        for similar in similar_genes:
            if similar["known_phenotypes"]:
                transferred = self._transfer_phenotypes(gene, similar)
                prediction["predicted_phenotypes"].extend(transferred)

        # 4. Predict downstream cascade
        cascade = self._predict_downstream_cascade(gene)
        prediction["downstream_cascade"] = cascade

        # 5. Identify compensation mechanisms
        compensation = self._find_compensatory_mechanisms(gene)
        prediction["potential_compensation"] = compensation

        # 6. Calculate confidence
        prediction["confidence_scores"] = self._calculate_confidence(prediction)

        # 7. Suggest validation experiments
        prediction["validation_experiments"] = self._suggest_experiments(gene, prediction)

        return prediction

    def _analyze_network_position(self, gene: str) -> Dict:
        """
        Analyze gene's position in signaling network.
        Critical for phenotype prediction.
        """
        position = {
            "pathways": [],
            "is_hub": False,                    # Many connections
            "is_bottleneck": False,             # Critical path node
            "upstream_count": 0,
            "downstream_count": 0,
            "betweenness_centrality": 0,
            "essential_for_pathways": [],       # Pathways that require this gene
        }

        # Count connections
        interactions = self.db.query_molecule(gene)
        position["upstream_count"] = len(interactions.get("upstream_regulators", []))
        position["downstream_count"] = len(interactions.get("downstream_targets", []))

        # Hub detection (many connections)
        total_connections = position["upstream_count"] + position["downstream_count"]
        if total_connections > 20:
            position["is_hub"] = True

        # Bottleneck detection (critical path)
        # Gene is bottleneck if removing it disconnects major pathway components

        return position

    def _predict_from_pathway(self, gene: str, pathway: Dict) -> List[Dict]:
        """
        Predict phenotypes based on pathway membership.
        """
        predictions = []

        # Get known phenotypes for other genes in same pathway
        pathway_genes = pathway.get("members", [])

        for pg in pathway_genes:
            if pg != gene:
                pg_phenotypes = self._get_known_phenotypes(pg, "mouse")
                for pheno in pg_phenotypes:
                    # Transfer phenotype with adjusted confidence
                    predictions.append({
                        "phenotype": pheno["phenotype"],
                        "source": f"inferred from {pg} in {pathway['name']}",
                        "confidence": pheno.get("severity_score", 0.5) * 0.7,  # Reduce confidence for inference
                        "mechanism": f"{gene} and {pg} are in same pathway ({pathway['name']})",
                    })

        return predictions

    def _predict_downstream_cascade(self, gene: str) -> List[Dict]:
        """
        Predict cascade of effects from KO.

        If gene A activates B, which activates C:
        KO of A -> reduced B activity -> reduced C activity
        """
        cascade = []

        # Get downstream targets
        targets = self.db.query_molecule(gene)["downstream_targets"]

        for target in targets:
            effect = {
                "target": target["symbol"],
                "direct_effect": target["interaction_type"],
                "predicted_change": self._predict_change(gene, target),
            }

            # Recursively get downstream of downstream
            secondary = self._get_secondary_effects(target["symbol"])
            effect["secondary_effects"] = secondary

            cascade.append(effect)

        return cascade

    def _find_compensatory_mechanisms(self, gene: str) -> List[Dict]:
        """
        Find potential compensatory mechanisms that may mask KO phenotype.
        """
        compensation = []

        # 1. Paralog compensation
        paralogs = self._get_paralogs(gene)
        for paralog in paralogs:
            if self._has_functional_overlap(gene, paralog):
                compensation.append({
                    "type": "paralog",
                    "gene": paralog,
                    "overlap_score": self._calculate_overlap(gene, paralog),
                    "note": f"{paralog} may compensate for {gene} loss",
                })

        # 2. Alternative pathway compensation
        pathways = self.db.query_molecule(gene)["pathway_context"]["pathways"]
        for pathway in pathways:
            alternatives = self._find_alternative_routes(gene, pathway)
            compensation.extend(alternatives)

        # 3. Feedback loop compensation
        feedback = self.db.query_molecule(gene)["pathway_context"]["feedback_loops"]
        for fb in feedback:
            if fb["type"] == "negative":
                compensation.append({
                    "type": "feedback_release",
                    "mechanism": fb,
                    "note": f"Loss of {gene} may release negative feedback on {fb['target']}",
                })

        return compensation

    def compare_ko_phenotypes(self, genes: List[str],
                              species: str = "mouse") -> Dict:
        """
        Compare predicted phenotypes across multiple KOs.

        Example: Compare VEGFR2 KO vs VEGFR3 KO vs DLL4 KO
        """
        comparison = {
            "genes": genes,
            "individual_predictions": {},
            "shared_phenotypes": [],
            "unique_phenotypes": {},
            "severity_ranking": [],
            "recommended_experiments": [],
        }

        # Get predictions for each gene
        all_phenotypes = {}
        for gene in genes:
            pred = self.predict_ko_phenotype(gene, species)
            comparison["individual_predictions"][gene] = pred
            all_phenotypes[gene] = set(p["phenotype"] for p in pred["predicted_phenotypes"])

        # Find shared phenotypes
        if len(genes) > 1:
            shared = all_phenotypes[genes[0]]
            for gene in genes[1:]:
                shared = shared.intersection(all_phenotypes[gene])
            comparison["shared_phenotypes"] = list(shared)

        # Find unique phenotypes
        for gene in genes:
            unique = all_phenotypes[gene] - set(comparison["shared_phenotypes"])
            comparison["unique_phenotypes"][gene] = list(unique)

        # Rank by severity
        severity_scores = []
        for gene in genes:
            pred = comparison["individual_predictions"][gene]
            avg_severity = np.mean([p.get("confidence", 0.5) for p in pred["predicted_phenotypes"]])
            severity_scores.append((gene, avg_severity))

        comparison["severity_ranking"] = sorted(severity_scores, key=lambda x: -x[1])

        return comparison
```

---

## 5. Nature Reviews Pathway Data Structure

```python
class NatureReviewsPathway:
    """
    Structured pathway data in Nature Reviews quality format.
    Ready for visualization and publication.
    """

    def __init__(self, pathway_name: str):
        self.name = pathway_name
        self.components = {
            "extracellular": [],
            "membrane": [],
            "cytoplasm": [],
            "nucleus": [],
        }
        self.interactions = []
        self.annotations = []

    def build_pathway_from_query(self, central_molecule: str,
                                 interaction_engine: InteractionQueryEngine) -> Dict:
        """
        Build complete pathway structure from a single molecule query.
        """

        # Get all interactions
        data = interaction_engine.get_all_interactions(central_molecule)

        pathway = {
            "name": f"{central_molecule} Signaling Network",
            "central_molecule": central_molecule,

            # Structured components
            "components": {
                "ligands": self._extract_ligands(data),
                "receptors": self._extract_receptors(data),
                "adaptors": self._extract_adaptors(data),
                "kinases": self._extract_kinases(data),
                "phosphatases": self._extract_phosphatases(data),
                "transcription_factors": self._extract_tfs(data),
                "target_genes": self._extract_targets(data),
            },

            # Interaction map
            "interactions": self._build_interaction_map(data),

            # Regulatory annotations
            "regulation": {
                "positive_regulators": data["upstream"]["activators"],
                "negative_regulators": data["upstream"]["inhibitors"],
                "feedback_loops": data["pathway_context"]["feedback_loops"],
            },

            # Crosstalk with other pathways
            "crosstalk": data["pathway_context"]["crosstalk"],

            # Disease relevance
            "disease_context": data["disease_context"],

            # Therapeutic opportunities
            "therapeutic_targets": self._identify_therapeutic_targets(data),

            # Publication-ready summary
            "summary": self._generate_summary(data),

            # Figure-ready data structure
            "figure_data": self._prepare_figure_data(data),
        }

        return pathway

    def _build_interaction_map(self, data: Dict) -> List[Dict]:
        """
        Build comprehensive interaction map for visualization.
        """
        interactions = []

        # Upstream -> target
        for activator in data["upstream"]["activators"]:
            interactions.append({
                "source": activator["symbol"],
                "target": data["query"],
                "type": "activation",
                "mechanism": activator.get("mechanism"),
                "evidence": activator.get("evidence_score", 0.5),
            })

        for inhibitor in data["upstream"]["inhibitors"]:
            interactions.append({
                "source": inhibitor["symbol"],
                "target": data["query"],
                "type": "inhibition",
                "mechanism": inhibitor.get("mechanism"),
                "evidence": inhibitor.get("evidence_score", 0.5),
            })

        # Target -> downstream
        for target in data["downstream"]["activates"]:
            interactions.append({
                "source": data["query"],
                "target": target["symbol"],
                "type": "activation",
                "mechanism": target.get("mechanism"),
                "evidence": target.get("evidence_score", 0.5),
            })

        for target in data["downstream"]["inhibits"]:
            interactions.append({
                "source": data["query"],
                "target": target["symbol"],
                "type": "inhibition",
                "mechanism": target.get("mechanism"),
                "evidence": target.get("evidence_score", 0.5),
            })

        return interactions

    def _prepare_figure_data(self, data: Dict) -> Dict:
        """
        Prepare data structure for Nature Reviews-style figure generation.
        """
        figure = {
            "layout": "hierarchical",  # or "circular", "force-directed"

            "layers": {
                "extracellular": {
                    "molecules": [],
                    "y_position": 4,
                },
                "membrane": {
                    "molecules": [],
                    "y_position": 3,
                },
                "cytoplasm": {
                    "molecules": [],
                    "y_position": 2,
                },
                "nucleus": {
                    "molecules": [],
                    "y_position": 1,
                },
            },

            "nodes": [],
            "edges": [],

            "legend": {
                "node_types": {
                    "ligand": {"shape": "circle", "color": "#22c55e"},
                    "receptor": {"shape": "Y", "color": "#6366f1"},
                    "kinase": {"shape": "diamond", "color": "#f59e0b"},
                    "phosphatase": {"shape": "diamond", "color": "#ef4444"},
                    "transcription_factor": {"shape": "rectangle", "color": "#8b5cf6"},
                    "target_gene": {"shape": "rectangle", "color": "#06b6d4"},
                },
                "edge_types": {
                    "activation": {"style": "solid", "arrow": "normal"},
                    "inhibition": {"style": "solid", "arrow": "tee"},
                    "binding": {"style": "dashed", "arrow": "none"},
                    "phosphorylation": {"style": "solid", "arrow": "normal", "label": "P"},
                },
            },

            "annotations": [],
        }

        # Populate nodes and edges from interaction data
        # ... implementation ...

        return figure

    def _generate_summary(self, data: Dict) -> str:
        """
        Generate publication-ready pathway summary.
        """
        query = data["query"]
        n_upstream = len(data["upstream"]["activators"]) + len(data["upstream"]["inhibitors"])
        n_downstream = len(data["downstream"]["activates"]) + len(data["downstream"]["inhibits"])
        n_pathways = len(data["pathway_context"]["pathways"])

        summary = f"""
        {query} Signaling Network Summary
        ================================

        {query} is regulated by {n_upstream} upstream molecules and controls
        {n_downstream} downstream targets. It participates in {n_pathways}
        signaling pathways with multiple crosstalk points.

        Key Activators: {', '.join([a['symbol'] for a in data['upstream']['activators'][:5]])}
        Key Inhibitors: {', '.join([i['symbol'] for i in data['upstream']['inhibitors'][:5]])}
        Key Targets: {', '.join([t['symbol'] for t in data['downstream']['activates'][:5]])}

        Pathway Crosstalk: {', '.join([c['to_pathway'] for c in data['pathway_context']['crosstalk'][:3]])}

        Disease Associations: {', '.join([d['disease'] for d in data['disease_context'][:3]])}
        """

        return summary
```

이 구조로 계속 확장할까요? 다음 단계:
1. 실제 데이터베이스 연동 (STRING, KEGG, Reactome API)
2. KO prediction engine 상세 구현
3. Pathway figure 생성 로직
