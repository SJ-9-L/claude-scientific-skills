# Disease-Specific Analysis Templates

## Overview

질환별 특화 분석 템플릿 - Cancer, Immune disorders, Cardiovascular diseases에 대한 pre-configured 분석 워크플로우

---

## 1. Cancer Analysis Templates

### 1.1 Tumor Microenvironment (TME) Analysis

```python
class TumorMicroenvironmentAnalyzer:
    """
    Comprehensive TME analysis from scRNA-seq data.
    Identifies immune infiltration, CAFs, and tumor-stroma interactions.
    """

    # TME Cell Type Markers
    TME_MARKERS = {
        # Immune cells
        "T_cells": {
            "CD8_cytotoxic": ["CD8A", "CD8B", "GZMB", "PRF1", "IFNG"],
            "CD8_exhausted": ["CD8A", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX"],
            "CD4_helper": ["CD4", "IL7R", "CCR7"],
            "Treg": ["CD4", "FOXP3", "IL2RA", "CTLA4", "IKZF2"],
            "Th1": ["CD4", "TBX21", "IFNG", "CXCR3"],
            "Th2": ["CD4", "GATA3", "IL4", "IL13"],
            "Th17": ["CD4", "RORC", "IL17A", "CCR6"],
        },
        "NK_cells": {
            "cytotoxic": ["NCAM1", "NKG7", "GNLY", "KLRD1", "KLRF1"],
            "regulatory": ["NCAM1", "KLRC1", "CD94"],
        },
        "Myeloid": {
            "M1_macrophage": ["CD68", "CD80", "CD86", "NOS2", "IL1B", "TNF"],
            "M2_macrophage": ["CD68", "CD163", "MRC1", "ARG1", "CD206"],
            "TAM": ["CD68", "MARCO", "MSR1", "TREM2", "APOE"],
            "MDSC": ["CD33", "S100A8", "S100A9", "ARG1"],
            "DC_conventional": ["CD1C", "CLEC9A", "XCR1", "BATF3"],
            "DC_plasmacytoid": ["IL3RA", "CLEC4C", "TCF4"],
        },
        "B_cells": {
            "naive": ["CD19", "MS4A1", "IGHD", "TCL1A"],
            "memory": ["CD19", "CD27", "IGHG1"],
            "plasma": ["SDC1", "XBP1", "PRDM1", "IGHG1"],
        },

        # Stromal cells
        "CAF": {
            "myofibroblastic": ["ACTA2", "FAP", "PDGFRA", "COL1A1", "TAGLN"],
            "inflammatory": ["FAP", "IL6", "CXCL12", "PDGFRA", "CFD"],
            "antigen_presenting": ["FAP", "CD74", "HLA-DRA", "SLPI"],
        },
        "Endothelial": {
            "tumor_EC": ["PECAM1", "VWF", "PLVAP", "ACKR1"],
            "tip_cell": ["PECAM1", "DLL4", "CXCR4", "APLN", "ESM1"],
            "stalk_cell": ["PECAM1", "JAG1", "NOTCH1", "HEY1"],
        },
        "Pericyte": ["PDGFRB", "RGS5", "ACTA2", "DES", "MCAM"],
    }

    # Immune checkpoint molecules
    IMMUNE_CHECKPOINTS = {
        "inhibitory": {
            "PD-1/PD-L1": {"receptor": "PDCD1", "ligands": ["CD274", "PDCD1LG2"]},
            "CTLA-4/CD80/86": {"receptor": "CTLA4", "ligands": ["CD80", "CD86"]},
            "LAG-3/MHC-II": {"receptor": "LAG3", "ligands": ["HLA-DRA", "HLA-DRB1"]},
            "TIM-3/Galectin-9": {"receptor": "HAVCR2", "ligands": ["LGALS9"]},
            "TIGIT/CD155": {"receptor": "TIGIT", "ligands": ["PVR", "PVRL2"]},
            "VISTA": {"receptor": "VSIR", "ligands": ["VSIR"]},
            "B7-H3": {"receptor": None, "ligands": ["CD276"]},
            "B7-H4": {"receptor": None, "ligands": ["VTCN1"]},
        },
        "stimulatory": {
            "CD28/CD80/86": {"receptor": "CD28", "ligands": ["CD80", "CD86"]},
            "4-1BB/4-1BBL": {"receptor": "TNFRSF9", "ligands": ["TNFSF9"]},
            "OX40/OX40L": {"receptor": "TNFRSF4", "ligands": ["TNFSF4"]},
            "GITR/GITRL": {"receptor": "TNFRSF18", "ligands": ["TNFSF18"]},
            "ICOS/ICOSL": {"receptor": "ICOS", "ligands": ["ICOSLG"]},
            "CD40/CD40L": {"receptor": "CD40", "ligands": ["CD40LG"]},
        }
    }

    def __init__(self, cancer_type: str):
        self.cancer_type = cancer_type
        self.cancer_specific_markers = self._get_cancer_markers()

    def _get_cancer_markers(self) -> dict:
        """Cancer type-specific tumor markers"""
        CANCER_MARKERS = {
            "breast": {
                "luminal_A": ["ESR1", "PGR", "GATA3", "FOXA1"],
                "luminal_B": ["ESR1", "MKI67", "CCNB1", "MYBL2"],
                "HER2+": ["ERBB2", "GRB7", "STARD3"],
                "basal": ["KRT5", "KRT14", "KRT17", "EGFR"],
                "claudin_low": ["VIM", "CDH2", "SNAI1", "TWIST1"],
            },
            "lung_NSCLC": {
                "adenocarcinoma": ["NKX2-1", "NAPSA", "KRT7", "MUC1"],
                "squamous": ["TP63", "KRT5", "KRT6A", "SOX2"],
                "EGFR_mutant": ["EGFR", "ERBB2", "ERBB3"],
                "KRAS_mutant": ["KRAS", "STK11", "KEAP1"],
                "ALK_fusion": ["ALK", "EML4"],
            },
            "colorectal": {
                "CMS1_immune": ["CD8A", "GZMB", "MLH1", "MSH2"],
                "CMS2_canonical": ["WNT", "MYC", "EGFR", "ERBB2"],
                "CMS3_metabolic": ["KRAS", "PIK3CA"],
                "CMS4_mesenchymal": ["VIM", "ZEB1", "SNAI1", "TGFβ"],
            },
            "melanoma": {
                "proliferative": ["MITF", "SOX10", "TYR", "MLANA"],
                "invasive": ["AXL", "WNT5A", "ZEB1"],
                "BRAF_mutant": ["BRAF", "ERK", "MEK"],
                "NRAS_mutant": ["NRAS", "PI3K"],
            },
            "pancreatic": {
                "classical": ["GATA6", "HNF1A", "HNF4A"],
                "basal_like": ["TP63", "KRT5", "KRT6A"],
                "KRAS_driven": ["KRAS", "CDKN2A", "TP53", "SMAD4"],
            },
            "ovarian": {
                "high_grade_serous": ["TP53", "BRCA1", "BRCA2", "PAX8"],
                "clear_cell": ["HNF1B", "ARID1A"],
                "endometrioid": ["CTNNB1", "PTEN", "PIK3CA"],
            },
            "glioblastoma": {
                "proneural": ["OLIG2", "PDGFRA", "IDH1"],
                "classical": ["EGFR", "PTEN"],
                "mesenchymal": ["CHI3L1", "CD44", "NF1"],
            },
        }
        return CANCER_MARKERS.get(self.cancer_type, {})

    def analyze_immune_infiltration(self, adata) -> dict:
        """
        Analyze immune cell infiltration in TME.

        Returns:
        - Cell type proportions
        - Immune score
        - Hot/Cold tumor classification
        - Spatial distribution (if spatial data)
        """
        results = {
            "cell_proportions": {},
            "immune_score": 0,
            "tumor_classification": "",
            "checkpoint_expression": {},
            "cytokine_profile": {},
        }

        # Calculate immune score
        immune_genes = ["CD3D", "CD8A", "GZMB", "PRF1", "IFNG", "CXCL9", "CXCL10"]
        suppressive_genes = ["FOXP3", "IL10", "TGFB1", "IDO1", "ARG1"]

        # Classify tumor
        # Hot: High CD8+ T cell infiltration, high IFN-gamma signature
        # Cold: Low immune infiltration
        # Excluded: Immune cells at margin but not infiltrating
        # Immunosuppressed: Infiltration but with suppressive signals

        return results

    def analyze_resistance_mechanisms(self, treatment: str) -> dict:
        """
        Analyze drug resistance mechanisms.

        treatment: "anti-PD1", "chemotherapy", "targeted_therapy", etc.
        """
        RESISTANCE_PATHWAYS = {
            "anti-PD1": {
                "primary": {
                    "antigen_presentation": ["B2M", "HLA-A", "HLA-B", "TAP1", "TAP2"],
                    "IFN_signaling": ["JAK1", "JAK2", "STAT1", "IRF1"],
                    "tumor_intrinsic": ["PTEN", "STK11", "KEAP1", "beta-catenin"],
                },
                "acquired": {
                    "immune_escape": ["B2M_loss", "HLA_LOH", "JAK_mutation"],
                    "alternative_checkpoints": ["LAG3", "TIM3", "VISTA", "TIGIT"],
                    "metabolic": ["IDO1", "ARG1", "adenosine_pathway"],
                },
            },
            "chemotherapy": {
                "drug_efflux": ["ABCB1", "ABCC1", "ABCG2"],
                "DNA_repair": ["ERCC1", "BRCA1", "BRCA2", "RAD51"],
                "apoptosis_evasion": ["BCL2", "MCL1", "XIAP", "survivin"],
                "EMT": ["VIM", "CDH2", "SNAI1", "ZEB1", "TWIST1"],
            },
            "EGFR_TKI": {
                "secondary_mutations": ["T790M", "C797S"],
                "bypass_activation": ["MET", "HER2", "AXL", "IGF1R"],
                "phenotype_switch": ["SCLC_transformation", "EMT"],
                "pathway_rewiring": ["FGFR", "SRC", "AKT"],
            },
        }
        return RESISTANCE_PATHWAYS.get(treatment, {})


### 1.2 Cancer Metastasis Analysis

```python
class MetastasisAnalyzer:
    """
    Analyze metastatic potential and metastasis-related pathways.
    """

    METASTASIS_SIGNATURES = {
        "EMT": {
            "epithelial": ["CDH1", "EPCAM", "KRT8", "KRT18", "KRT19", "CLDN1"],
            "mesenchymal": ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1"],
            "hybrid_EMT": ["partial_E", "partial_M"],  # Hybrid phenotype markers
        },
        "invasion": {
            "proteases": ["MMP2", "MMP9", "MMP14", "ADAM17", "CTSL"],
            "motility": ["RAC1", "CDC42", "RHOA", "ROCK1", "LIMK1"],
            "adhesion_switch": ["ITGAV", "ITGB3", "ITGB1", "CD44"],
        },
        "intravasation": {
            "vascular_mimicry": ["VE-cadherin", "VEGFA", "LAMC2"],
            "endothelial_interaction": ["PECAM1", "ICAM1", "VCAM1"],
        },
        "survival_circulation": {
            "anoikis_resistance": ["BCL2", "BCLXL", "TrkB"],
            "platelet_interaction": ["TBXA2R", "P2RY12"],
        },
        "extravasation": {
            "adhesion": ["SELE", "SELP", "ICAM1"],
            "transendothelial": ["MMP2", "MMP9", "PECAM1"],
        },
        "colonization": {
            "pre_metastatic_niche": ["LOX", "VEGFA", "S100A8", "S100A9"],
            "dormancy": ["NR2F1", "DEC2", "p27"],
            "reactivation": ["VCAM1", "POSTN", "TNC"],
        },

        # Organ-specific metastasis signatures
        "organ_tropism": {
            "bone": ["CXCR4", "RANKL", "PTHrP", "CTGF", "IL11"],
            "lung": ["VCAM1", "ANGPTL4", "MMP1", "MMP2"],
            "liver": ["CEACAM1", "KRT19", "claudins"],
            "brain": ["ST6GALNAC5", "COX2", "HBEGF", "L1CAM"],
            "lymph_node": ["CCR7", "VEGFC", "VEGFD"],
        }
    }

    def analyze_metastatic_potential(self, adata, primary_site: str) -> dict:
        """
        Score metastatic potential and predict likely metastatic sites.
        """
        return {
            "EMT_score": 0,
            "invasion_score": 0,
            "predicted_sites": [],
            "dormancy_signature": False,
            "therapeutic_vulnerabilities": [],
        }
```

---

## 2. Immune Disorder Templates

### 2.1 Autoimmune Disease Analysis

```python
class AutoimmuneAnalyzer:
    """
    Analyze autoimmune disease mechanisms from scRNA-seq.
    """

    AUTOIMMUNE_SIGNATURES = {
        # T cell abnormalities
        "T_cell_dysregulation": {
            "Th17_expansion": ["IL17A", "IL17F", "IL22", "RORC", "CCR6", "IL23R"],
            "Treg_dysfunction": ["FOXP3", "IL2RA", "CTLA4", "functional_markers"],
            "Tfh_expansion": ["CXCR5", "PD1", "BCL6", "IL21", "ICOS"],
            "exhaustion": ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"],
        },

        # B cell abnormalities
        "B_cell_dysregulation": {
            "autoantibody_producing": ["plasma_markers", "autoreactive_BCR"],
            "age_associated_B": ["TBET", "CD11c", "CD21low"],
            "memory_expansion": ["CD27", "class_switched"],
        },

        # Innate immune
        "innate_dysregulation": {
            "type_I_IFN": ["ISG15", "MX1", "IFIT1", "OAS1", "IFI44"],
            "inflammasome": ["NLRP3", "CASP1", "IL1B", "IL18", "GSDMD"],
            "neutrophil_NETs": ["MPO", "ELANE", "PR3", "citrullinated_histones"],
        },

        # Disease-specific
        "disease_markers": {
            "rheumatoid_arthritis": {
                "synovial_fibroblast": ["PDPN", "FAP", "THY1", "CDH11"],
                "citrullination": ["PAD4", "citrullinated_proteins"],
                "bone_erosion": ["RANKL", "MMP3", "cathepsin_K"],
            },
            "SLE": {
                "IFN_signature": ["MX1", "OAS1", "ISG15", "IFI44", "SIGLEC1"],
                "B_cell_hyperactivity": ["BAFF", "APRIL", "BLyS"],
                "complement": ["C3", "C4", "C1q"],
            },
            "multiple_sclerosis": {
                "BBB_disruption": ["MMP9", "CXCL10", "CCL2"],
                "demyelination": ["MBP", "MOG", "reactive_astrocyte"],
                "Th17_CNS": ["IL17", "GM-CSF", "CCR6"],
            },
            "IBD": {
                "barrier_dysfunction": ["MUC2", "tight_junctions", "antimicrobial"],
                "Th17_mucosal": ["IL17A", "IL22", "IL23R"],
                "macrophage_M1": ["TNF", "IL1B", "IL6", "NOS2"],
            },
            "type1_diabetes": {
                "beta_cell_stress": ["INS", "PDX1", "NKX6.1", "stress_markers"],
                "islet_infiltration": ["CD8", "CD4", "autoreactive_T"],
            },
            "psoriasis": {
                "IL17_axis": ["IL17A", "IL17F", "IL22", "IL23"],
                "keratinocyte": ["KRT16", "S100A7", "S100A8", "S100A9"],
            },
        }
    }

    THERAPEUTIC_TARGETS = {
        "TNF_blockade": {
            "targets": ["TNF", "TNFR1", "TNFR2"],
            "drugs": ["Infliximab", "Adalimumab", "Etanercept", "Golimumab", "Certolizumab"],
            "indications": ["RA", "IBD", "Psoriasis", "AS"],
        },
        "IL6_blockade": {
            "targets": ["IL6", "IL6R", "gp130"],
            "drugs": ["Tocilizumab", "Sarilumab", "Siltuximab"],
            "indications": ["RA", "JIA", "CRS"],
        },
        "IL17_blockade": {
            "targets": ["IL17A", "IL17F", "IL17RA"],
            "drugs": ["Secukinumab", "Ixekizumab", "Brodalumab", "Bimekizumab"],
            "indications": ["Psoriasis", "PsA", "AS"],
        },
        "IL23_blockade": {
            "targets": ["IL23p19", "IL12/23p40"],
            "drugs": ["Guselkumab", "Risankizumab", "Tildrakizumab", "Ustekinumab"],
            "indications": ["Psoriasis", "PsA", "IBD"],
        },
        "JAK_inhibition": {
            "targets": ["JAK1", "JAK2", "JAK3", "TYK2"],
            "drugs": ["Tofacitinib", "Baricitinib", "Upadacitinib", "Filgotinib", "Deucravacitinib"],
            "indications": ["RA", "UC", "AD", "Psoriasis"],
        },
        "B_cell_depletion": {
            "targets": ["CD20", "CD19", "BAFF"],
            "drugs": ["Rituximab", "Ocrelizumab", "Ofatumumab", "Belimumab"],
            "indications": ["RA", "MS", "SLE"],
        },
        "T_cell_costim_blockade": {
            "targets": ["CD80/86", "CTLA4-Ig"],
            "drugs": ["Abatacept"],
            "indications": ["RA", "JIA"],
        },
    }

    def analyze_disease_activity(self, adata, disease: str) -> dict:
        """Analyze disease activity score from single-cell data."""
        pass

    def predict_treatment_response(self, adata, drug_class: str) -> dict:
        """Predict response to specific drug class."""
        pass


### 2.2 Inflammation Analysis

```python
class InflammationAnalyzer:
    """
    Analyze inflammatory pathways and resolution mechanisms.
    """

    INFLAMMATION_PATHWAYS = {
        "initiation": {
            "PRR_sensing": {
                "TLRs": ["TLR2", "TLR4", "TLR7", "TLR9"],
                "NLRs": ["NOD1", "NOD2", "NLRP3"],
                "RLRs": ["RIG-I", "MDA5"],
                "cGAS-STING": ["CGAS", "STING1", "TBK1", "IRF3"],
            },
            "NFkB_activation": ["NFKB1", "RELA", "IKBKB", "NFKBIA"],
            "inflammasome": ["NLRP3", "ASC", "CASP1", "GSDMD", "IL1B", "IL18"],
        },

        "amplification": {
            "cytokines": {
                "pro_inflammatory": ["IL1B", "IL6", "TNF", "IL12", "IL23", "IL17", "IFNG"],
                "chemokines": ["CXCL8", "CCL2", "CCL3", "CCL5", "CXCL10"],
            },
            "lipid_mediators": {
                "pro_inflammatory": ["prostaglandins", "leukotrienes", "PAF"],
                "enzymes": ["PTGS2", "ALOX5", "PLA2"],
            },
        },

        "resolution": {
            "SPMs": {
                "resolvins": ["RvD1", "RvD2", "RvE1"],
                "lipoxins": ["LXA4", "LXB4"],
                "protectins": ["PD1", "PDX"],
                "maresins": ["MaR1", "MaR2"],
            },
            "anti_inflammatory_cytokines": ["IL10", "TGFB1", "IL1RA", "IL37"],
            "efferocytosis": ["MERTK", "AXL", "TYRO3", "CD36", "TIM4"],
            "M2_macrophage": ["ARG1", "MRC1", "CD163", "RETNLA"],
        },

        "chronic_inflammation": {
            "markers": ["CRP", "SAA", "fibrinogen"],
            "tissue_remodeling": ["MMPs", "TIMPs", "fibrosis_markers"],
            "metabolic_inflammation": ["obesity", "insulin_resistance", "NAFLD"],
        }
    }

    def calculate_inflammation_score(self, adata) -> dict:
        """Calculate composite inflammation score."""
        return {
            "acute_score": 0,
            "chronic_score": 0,
            "resolution_capacity": 0,
            "dominant_pathway": "",
        }
```

---

## 3. Cardiovascular Disease Templates

### 3.1 Vascular Analysis

```python
class VascularAnalyzer:
    """
    Analyze vascular biology and angiogenesis from scRNA-seq.
    """

    VASCULAR_CELL_TYPES = {
        "endothelial": {
            "arterial": ["EFNB2", "DLL4", "GJA5", "GJA4", "HEY1", "NOTCH1", "CXCR4"],
            "venous": ["NR2F2", "EPHB4", "APLNR", "ACKR1"],
            "capillary": ["RGCC", "CA4", "FCN3"],
            "tip_cell": ["DLL4", "CXCR4", "APLN", "ESM1", "ANGPT2", "PDGFB"],
            "stalk_cell": ["JAG1", "NOTCH1", "HEY1", "HEY2", "VEGFR1"],
            "lymphatic": ["PROX1", "LYVE1", "PDPN", "FLT4", "CCL21"],
            "HEV": ["ACKR1", "GLYCAM1", "CHST4", "MADCAM1"],
        },
        "mural": {
            "pericyte": ["PDGFRB", "RGS5", "KCNJ8", "ABCC9", "DES"],
            "vSMC": ["ACTA2", "MYH11", "CNN1", "TAGLN", "SMTN"],
        },
        "specialized": {
            "fenestrated": ["PLVAP", "liver_sinusoid", "glomerular"],
            "BBB": ["CLDN5", "SLC2A1", "MFSD2A", "tight_junctions"],
        }
    }

    ANGIOGENESIS_PATHWAYS = {
        "VEGF": {
            "ligands": ["VEGFA", "VEGFB", "VEGFC", "VEGFD", "PlGF"],
            "receptors": ["KDR", "FLT1", "FLT4"],
            "coreceptors": ["NRP1", "NRP2"],
            "downstream": ["PLCG1", "AKT", "ERK", "SRC", "FAK"],
        },
        "Notch": {
            "ligands": ["DLL4", "DLL1", "JAG1", "JAG2"],
            "receptors": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"],
            "downstream": ["HES1", "HEY1", "HEY2", "HEYL"],
        },
        "Angiopoietin": {
            "ligands": ["ANGPT1", "ANGPT2"],
            "receptor": ["TEK"],
            "function": {"ANGPT1": "vessel_stabilization", "ANGPT2": "vessel_destabilization"},
        },
        "Ephrin": {
            "arterial": {"ligand": "EFNB2", "receptor": "EPHB4"},
            "function": "arteriovenous_specification",
        },
        "TGFβ_BMP": {
            "ligands": ["TGFB1", "BMP9", "BMP10"],
            "receptors": ["ACVRL1", "ENG", "TGFBR2"],
            "downstream": ["SMAD1", "SMAD5", "SMAD2", "SMAD3"],
        },
        "Wnt": {
            "ligands": ["WNT7A", "WNT7B"],
            "receptors": ["FZD4", "LRP5", "LRP6"],
            "downstream": ["CTNNB1", "LEF1", "TCF7"],
            "function": "BBB_formation",
        },
    }

    VASCULAR_PATHOLOGIES = {
        "atherosclerosis": {
            "stages": {
                "initiation": ["endothelial_dysfunction", "LDL_retention", "monocyte_recruitment"],
                "progression": ["foam_cell_formation", "SMC_migration", "fibrous_cap"],
                "complication": ["plaque_rupture", "thrombosis", "calcification"],
            },
            "cell_markers": {
                "activated_EC": ["VCAM1", "ICAM1", "SELE", "CCL2"],
                "foam_cell": ["CD68", "ABCA1", "CD36", "MSR1"],
                "plaque_macrophage": ["TREM2", "GPNMB", "SPP1", "LGALS3"],
            },
        },
        "aneurysm": {
            "markers": ["MMP2", "MMP9", "elastin_degradation", "SMC_apoptosis"],
            "pathways": ["TGFβ", "NOTCH", "inflammation"],
        },
        "vascular_malformation": {
            "AVM": ["KRAS", "BRAF", "MAP2K1"],
            "HHT": ["ENG", "ACVRL1", "SMAD4"],
            "CCM": ["CCM1", "CCM2", "CCM3"],
        },
    }

    def analyze_angiogenic_state(self, adata) -> dict:
        """Analyze angiogenic activation state."""
        pass

    def compare_ec_phenotypes(self, adata, conditions: list) -> dict:
        """Compare endothelial cell phenotypes across conditions."""
        pass


### 3.2 Cardiac Analysis

```python
class CardiacAnalyzer:
    """
    Analyze cardiac cell types and disease states.
    """

    CARDIAC_CELL_TYPES = {
        "cardiomyocyte": {
            "ventricular": ["MYH7", "MYL2", "TNNT2", "TNNI3", "NPPA", "NPPB"],
            "atrial": ["MYL7", "MYL4", "NPPA", "PITX2"],
            "nodal": ["HCN4", "TBX3", "SHOX2"],
            "purkinje": ["GJA5", "SCN5A", "IRX3"],
        },
        "cardiac_fibroblast": {
            "quiescent": ["TCF21", "PDGFRA", "DDR2"],
            "activated": ["POSTN", "FAP", "ACTA2", "COL1A1"],
            "myofibroblast": ["ACTA2", "TAGLN", "COL1A1", "COL3A1"],
        },
        "endocardial": ["NPR3", "CDH5", "EMCN"],
        "epicardial": ["WT1", "TBX18", "TCF21", "RALDH2"],
        "immune": {
            "resident_macrophage": ["LYVE1", "TIMD4", "MRC1", "F13A1"],
            "recruited_macrophage": ["CCR2", "LY6C", "CD14"],
        },
    }

    CARDIAC_PATHWAYS = {
        "hypertrophy": {
            "pathological": ["NFAT", "calcineurin", "HDAC", "MEF2"],
            "physiological": ["PI3K", "AKT", "mTOR"],
            "markers": ["NPPA", "NPPB", "MYH7", "ACTA1"],
        },
        "fibrosis": {
            "drivers": ["TGFB1", "CTGF", "PDGF", "angiotensin_II"],
            "markers": ["COL1A1", "COL3A1", "FN1", "POSTN", "ACTA2"],
        },
        "ischemia_reperfusion": {
            "cell_death": ["apoptosis", "necroptosis", "ferroptosis", "pyroptosis"],
            "inflammation": ["DAMPs", "NLRP3", "IL1B", "IL6"],
            "repair": ["angiogenesis", "scar_formation", "regeneration"],
        },
        "arrhythmia": {
            "ion_channels": ["SCN5A", "KCNQ1", "KCNH2", "KCNJ2", "CACNA1C"],
            "gap_junctions": ["GJA1", "GJA5"],
            "calcium_handling": ["RYR2", "ATP2A2", "PLN", "CASQ2"],
        },
    }

    def analyze_cardiac_remodeling(self, adata, condition: str) -> dict:
        """Analyze cardiac remodeling patterns."""
        pass
```

---

## 4. Usage Examples

### 4.1 Tumor Microenvironment Analysis

```python
from disease_specific_templates import TumorMicroenvironmentAnalyzer, MetastasisAnalyzer

# Initialize analyzer
tme = TumorMicroenvironmentAnalyzer(cancer_type="breast")

# Load scRNA-seq data
adata = sc.read_h5ad("breast_cancer_scrnaseq.h5ad")

# Analyze immune infiltration
infiltration = tme.analyze_immune_infiltration(adata)
print(f"Tumor classification: {infiltration['tumor_classification']}")
print(f"Immune score: {infiltration['immune_score']}")
print(f"CD8+ T cell fraction: {infiltration['cell_proportions']['CD8_T']}")

# Analyze checkpoint expression
checkpoints = infiltration['checkpoint_expression']
# Identifies patients likely to respond to immunotherapy

# Analyze resistance to anti-PD1
resistance = tme.analyze_resistance_mechanisms("anti-PD1")
# Returns primary and acquired resistance pathways

# Metastasis analysis
meta = MetastasisAnalyzer()
metastatic_potential = meta.analyze_metastatic_potential(adata, "breast")
print(f"EMT score: {metastatic_potential['EMT_score']}")
print(f"Predicted metastatic sites: {metastatic_potential['predicted_sites']}")
```

### 4.2 Autoimmune Disease Analysis

```python
from disease_specific_templates import AutoimmuneAnalyzer

# Initialize for rheumatoid arthritis
analyzer = AutoimmuneAnalyzer(disease="rheumatoid_arthritis")

# Load synovial tissue scRNA-seq
adata = sc.read_h5ad("synovial_scrnaseq.h5ad")

# Analyze disease activity
activity = analyzer.analyze_disease_activity(adata)
print(f"Inflammation score: {activity['inflammation_score']}")
print(f"Fibroblast activation: {activity['fibroblast_activation']}")

# Predict treatment response
anti_tnf_response = analyzer.predict_treatment_response(adata, "TNF_blockade")
jak_inhibitor_response = analyzer.predict_treatment_response(adata, "JAK_inhibition")

# Compare and recommend
if anti_tnf_response['predicted_response'] > jak_inhibitor_response['predicted_response']:
    print("Recommend TNF blockade")
else:
    print("Recommend JAK inhibition")
```

### 4.3 Vascular Analysis

```python
from disease_specific_templates import VascularAnalyzer

# Initialize analyzer
vascular = VascularAnalyzer()

# Analyze tumor angiogenesis
tumor_adata = sc.read_h5ad("tumor_ec_scrnaseq.h5ad")

# Compare to normal vasculature
normal_adata = sc.read_h5ad("normal_ec_scrnaseq.h5ad")

comparison = vascular.compare_ec_phenotypes(
    [tumor_adata, normal_adata],
    conditions=["tumor", "normal"]
)

# Identify tumor-specific EC markers
print("Tumor EC markers:", comparison['differential_markers'])
print("Tip cell fraction:", comparison['tip_cell_fraction'])
print("Angiogenic score:", comparison['angiogenic_score'])

# Analyze response to anti-VEGF
anti_vegf_response = vascular.predict_anti_vegf_response(tumor_adata)
print(f"Predicted response: {anti_vegf_response['response_probability']}")
print(f"Resistance mechanisms: {anti_vegf_response['resistance_pathways']}")
```

---

## 5. Integration with Main Pipeline

```python
# Complete analysis workflow
from scripts.scrna_pipeline import ScRNAAnalyzer
from scripts.pathway_analyzer import PathwayDeepDive
from disease_specific_templates import TumorMicroenvironmentAnalyzer

# Step 1: Standard scRNA-seq analysis
analyzer = ScRNAAnalyzer()
adata = analyzer.full_pipeline("tumor_sample.h5ad")

# Step 2: Disease-specific annotation
tme = TumorMicroenvironmentAnalyzer(cancer_type="lung_NSCLC")
immune_analysis = tme.analyze_immune_infiltration(adata)

# Step 3: Pathway deep dive
pathway = PathwayDeepDive(
    target="EGFR signaling",
    context="lung_adenocarcinoma"
)
resistance = pathway.analyze_resistance(treatment="EGFR_TKI")

# Step 4: Generate comprehensive report
report = {
    "immune_contexture": immune_analysis,
    "molecular_subtype": tme.classify_subtype(adata),
    "therapeutic_targets": pathway.get_druggable_targets(),
    "resistance_mechanisms": resistance,
    "clinical_recommendations": tme.generate_recommendations(adata)
}
```
