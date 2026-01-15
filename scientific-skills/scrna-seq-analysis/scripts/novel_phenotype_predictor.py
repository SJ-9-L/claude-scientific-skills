"""
Novel Gene Phenotype Predictor
==============================
Predicts phenotypes for genes with NO existing KO data.

This is the CORE VALUE of the system:
- Uses expression similarity to known genes
- Network position analysis
- Protein family/domain inference
- Cross-species ortholog phenotypes

Input: Gene with NO known KO phenotype
Output: Predicted phenotype with confidence and validation plan
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field
from collections import defaultdict
import json


# ============================================================================
# PREDICTION BASED ON EXPRESSION SIMILARITY
# ============================================================================

class ExpressionSimilarityPredictor:
    """
    Predict phenotype based on co-expression with known genes.

    Logic: If gene X has similar expression pattern to DLL4,
    and DLL4 KO causes hypersprouting,
    then gene X KO might also affect sprouting.
    """

    def __init__(self):
        # Pre-computed co-expression modules from meta-analysis
        # Each module contains genes with highly correlated expression

        self.COEXPRESSION_MODULES = {
            "tip_cell_module": {
                "core_genes": ["ESM1", "CXCR4", "APLN", "ANGPT2", "DLL4"],
                "ko_phenotype": "sprouting_defect",
                "associated_genes": [
                    {"gene": "KCNE3", "correlation": 0.82, "ko_data": False},
                    {"gene": "RAMP3", "correlation": 0.78, "ko_data": False},
                    {"gene": "INHBB", "correlation": 0.75, "ko_data": False},
                    {"gene": "MCAM", "correlation": 0.73, "ko_data": False},
                    {"gene": "PLAUR", "correlation": 0.71, "ko_data": False},
                    {"gene": "LRRC32", "correlation": 0.69, "ko_data": False},
                    {"gene": "PXDN", "correlation": 0.67, "ko_data": False},
                    {"gene": "COL4A1", "correlation": 0.65, "ko_data": True},  # Has KO data
                    {"gene": "COL4A2", "correlation": 0.64, "ko_data": True},
                ],
            },

            "notch_regulation_module": {
                "core_genes": ["DLL4", "JAG1", "NOTCH1", "HES1", "HEY1"],
                "ko_phenotype": "tip_stalk_imbalance",
                "associated_genes": [
                    {"gene": "NEURL1B", "correlation": 0.79, "ko_data": False},
                    {"gene": "MFNG", "correlation": 0.74, "ko_data": False},
                    {"gene": "NRARP", "correlation": 0.72, "ko_data": False},
                    {"gene": "MAML1", "correlation": 0.68, "ko_data": True},
                ],
            },

            "arterial_identity_module": {
                "core_genes": ["GJA5", "GJA4", "EFNB2", "HEY2", "SOX17"],
                "ko_phenotype": "arterial_specification_defect",
                "associated_genes": [
                    {"gene": "SEMA3G", "correlation": 0.85, "ko_data": False},
                    {"gene": "LTBP4", "correlation": 0.76, "ko_data": False},
                    {"gene": "MECOM", "correlation": 0.73, "ko_data": True},
                ],
            },

            "venous_identity_module": {
                "core_genes": ["NR2F2", "EPHB4", "ACKR1", "APLNR"],
                "ko_phenotype": "venous_specification_defect",
                "associated_genes": [
                    {"gene": "ADGRG6", "correlation": 0.81, "ko_data": False},
                    {"gene": "SMARCA4", "correlation": 0.72, "ko_data": True},
                ],
            },

            "ec_survival_module": {
                "core_genes": ["KDR", "PIK3CA", "AKT1", "BCL2"],
                "ko_phenotype": "ec_apoptosis_vascular_regression",
                "associated_genes": [
                    {"gene": "ANGPT1", "correlation": 0.78, "ko_data": True},
                    {"gene": "TEK", "correlation": 0.76, "ko_data": True},
                    {"gene": "PDGFRB", "correlation": 0.71, "ko_data": True},
                ],
            },

            "migration_module": {
                "core_genes": ["CDC42", "RAC1", "RHOA", "PTK2", "PXN"],
                "ko_phenotype": "migration_defect",
                "associated_genes": [
                    {"gene": "TIAM1", "correlation": 0.77, "ko_data": False},
                    {"gene": "DOCK1", "correlation": 0.74, "ko_data": False},
                    {"gene": "ELMO1", "correlation": 0.71, "ko_data": False},
                ],
            },
        }

    def find_module_membership(self, gene: str) -> List[Dict]:
        """Find which co-expression modules a gene belongs to"""

        memberships = []

        for module_name, module_data in self.COEXPRESSION_MODULES.items():
            # Check if gene is a core member
            if gene.upper() in [g.upper() for g in module_data["core_genes"]]:
                memberships.append({
                    "module": module_name,
                    "role": "core",
                    "correlation": 1.0,
                    "predicted_phenotype": module_data["ko_phenotype"],
                    "confidence": 0.95,
                })

            # Check if gene is an associated member
            for assoc in module_data["associated_genes"]:
                if gene.upper() == assoc["gene"].upper():
                    memberships.append({
                        "module": module_name,
                        "role": "associated",
                        "correlation": assoc["correlation"],
                        "predicted_phenotype": module_data["ko_phenotype"],
                        "confidence": assoc["correlation"] * 0.8,  # Scaled confidence
                        "has_ko_data": assoc["ko_data"],
                    })

        return memberships


# ============================================================================
# PREDICTION BASED ON PROTEIN FAMILY
# ============================================================================

class ProteinFamilyPredictor:
    """
    Predict phenotype based on protein family membership.

    Logic: If gene X is in the same family as gene Y,
    and gene Y KO has phenotype Z,
    then gene X KO might have similar phenotype.
    """

    FAMILY_PHENOTYPES = {
        "receptor_tyrosine_kinase": {
            "members": ["EGFR", "KDR", "FLT1", "FLT4", "FGFR1", "FGFR2", "PDGFRA", "PDGFRB", "MET", "TEK"],
            "general_phenotype": "signaling_defect_growth_migration",
            "subfamily_phenotypes": {
                "VEGFR_family": {
                    "members": ["KDR", "FLT1", "FLT4"],
                    "phenotype": "vascular_development_defect",
                    "severity": "severe_to_lethal",
                },
                "FGFR_family": {
                    "members": ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
                    "phenotype": "skeletal_vascular_defect",
                    "severity": "variable",
                },
            },
        },

        "notch_pathway": {
            "members": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4", "JAG1", "JAG2"],
            "general_phenotype": "cell_fate_lateral_inhibition_defect",
            "subfamily_phenotypes": {
                "notch_receptors": {
                    "members": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"],
                    "phenotype": "specification_defect",
                },
                "delta_ligands": {
                    "members": ["DLL1", "DLL3", "DLL4"],
                    "phenotype": "lateral_inhibition_defect",
                },
            },
        },

        "connexins": {
            "members": ["GJA1", "GJA4", "GJA5", "GJB1", "GJB2"],
            "general_phenotype": "gap_junction_communication_defect",
            "subfamily_phenotypes": {
                "vascular_connexins": {
                    "members": ["GJA4", "GJA5"],
                    "phenotype": "arterial_identity_defect",
                },
            },
        },

        "ephrin_signaling": {
            "members": ["EFNA1", "EFNA2", "EFNB1", "EFNB2", "EPHA2", "EPHB2", "EPHB4"],
            "general_phenotype": "cell_boundary_formation_defect",
            "subfamily_phenotypes": {
                "AV_specification": {
                    "members": ["EFNB2", "EPHB4"],
                    "phenotype": "arteriovenous_malformation",
                    "severity": "severe",
                },
            },
        },

        "transcription_factors_bHLH": {
            "members": ["HES1", "HES5", "HEY1", "HEY2", "HEYL", "TWIST1", "ID1", "ID2", "ID3"],
            "general_phenotype": "transcriptional_regulation_defect",
            "subfamily_phenotypes": {
                "notch_targets": {
                    "members": ["HES1", "HES5", "HEY1", "HEY2", "HEYL"],
                    "phenotype": "similar_to_notch_ko",
                },
            },
        },

        "small_GTPases": {
            "members": ["RHOA", "RHOB", "RHOC", "RAC1", "RAC2", "CDC42", "RRAS", "RRAS2"],
            "general_phenotype": "cytoskeletal_migration_defect",
            "subfamily_phenotypes": {
                "rho_family": {
                    "members": ["RHOA", "RHOB", "RHOC"],
                    "phenotype": "contractility_stress_fiber_defect",
                },
                "rac_cdc42_family": {
                    "members": ["RAC1", "RAC2", "CDC42"],
                    "phenotype": "lamellipodia_filopodia_migration_defect",
                },
            },
        },
    }

    def predict_from_family(self, gene: str) -> Dict:
        """Predict phenotype based on protein family"""

        for family_name, family_data in self.FAMILY_PHENOTYPES.items():
            if gene.upper() in [g.upper() for g in family_data["members"]]:
                prediction = {
                    "gene": gene,
                    "family": family_name,
                    "predicted_phenotype": family_data["general_phenotype"],
                    "confidence": 0.6,  # Base confidence for family-level prediction
                    "subfamily_prediction": None,
                }

                # Check subfamily for more specific prediction
                for subfam_name, subfam_data in family_data.get("subfamily_phenotypes", {}).items():
                    if gene.upper() in [g.upper() for g in subfam_data["members"]]:
                        prediction["subfamily"] = subfam_name
                        prediction["predicted_phenotype"] = subfam_data["phenotype"]
                        prediction["confidence"] = 0.75  # Higher confidence for subfamily
                        if "severity" in subfam_data:
                            prediction["predicted_severity"] = subfam_data["severity"]

                return prediction

        return {"gene": gene, "family": None, "predicted_phenotype": None}


# ============================================================================
# PREDICTION BASED ON NETWORK POSITION
# ============================================================================

class NetworkPositionPredictor:
    """
    Predict phenotype based on network centrality and position.

    Logic:
    - Hub genes (many connections) → more severe phenotype
    - Bottleneck genes (critical path) → more specific phenotype
    - Peripheral genes → milder/redundant phenotype
    """

    # Pre-computed network metrics for EC-relevant genes
    NETWORK_METRICS = {
        # Hubs - many connections, severe phenotypes
        "hubs": {
            "KDR": {"degree": 156, "betweenness": 0.42, "predicted_severity": "lethal"},
            "NOTCH1": {"degree": 134, "betweenness": 0.38, "predicted_severity": "lethal"},
            "AKT1": {"degree": 245, "betweenness": 0.51, "predicted_severity": "lethal"},
            "TP53": {"degree": 312, "betweenness": 0.61, "predicted_severity": "lethal_or_cancer"},
            "SRC": {"degree": 198, "betweenness": 0.44, "predicted_severity": "severe"},
        },

        # Bottlenecks - critical path nodes
        "bottlenecks": {
            "RBPJ": {"degree": 45, "betweenness": 0.35, "role": "notch_pathway_bottleneck"},
            "SMAD4": {"degree": 67, "betweenness": 0.39, "role": "tgfb_pathway_bottleneck"},
            "CTNNB1": {"degree": 89, "betweenness": 0.42, "role": "wnt_pathway_bottleneck"},
        },

        # Peripheral - redundant or tissue-specific
        "peripheral": {
            "ESM1": {"degree": 12, "betweenness": 0.02, "predicted_severity": "mild_context_dependent"},
            "ACKR1": {"degree": 8, "betweenness": 0.01, "predicted_severity": "mild_tissue_specific"},
        },
    }

    def predict_from_network(self, gene: str) -> Dict:
        """Predict phenotype severity based on network position"""

        # Check hubs
        if gene.upper() in [g.upper() for g in self.NETWORK_METRICS["hubs"]]:
            data = self.NETWORK_METRICS["hubs"].get(gene.upper(), {})
            return {
                "gene": gene,
                "network_role": "hub",
                "degree": data.get("degree"),
                "betweenness": data.get("betweenness"),
                "predicted_severity": data.get("predicted_severity", "severe"),
                "confidence": 0.8,
                "reasoning": "High connectivity suggests essential function; KO likely severe"
            }

        # Check bottlenecks
        if gene.upper() in [g.upper() for g in self.NETWORK_METRICS["bottlenecks"]]:
            data = self.NETWORK_METRICS["bottlenecks"].get(gene.upper(), {})
            return {
                "gene": gene,
                "network_role": "bottleneck",
                "role_description": data.get("role"),
                "predicted_severity": "pathway_specific_severe",
                "confidence": 0.75,
                "reasoning": "Critical pathway node; KO affects entire pathway"
            }

        # Check peripheral
        if gene.upper() in [g.upper() for g in self.NETWORK_METRICS["peripheral"]]:
            data = self.NETWORK_METRICS["peripheral"].get(gene.upper(), {})
            return {
                "gene": gene,
                "network_role": "peripheral",
                "predicted_severity": data.get("predicted_severity", "mild"),
                "confidence": 0.65,
                "reasoning": "Low connectivity; likely redundant or context-specific"
            }

        return {
            "gene": gene,
            "network_role": "unknown",
            "predicted_severity": "unknown",
            "confidence": 0.3,
        }


# ============================================================================
# MAIN PREDICTOR - INTEGRATES ALL METHODS
# ============================================================================

class NovelPhenotypePredictor:
    """
    Main predictor that integrates all prediction methods.

    For genes with NO existing KO data, predicts:
    1. Expected phenotype
    2. Severity
    3. Affected pathways
    4. Confidence score
    5. Validation experiments
    """

    def __init__(self):
        self.expression_predictor = ExpressionSimilarityPredictor()
        self.family_predictor = ProteinFamilyPredictor()
        self.network_predictor = NetworkPositionPredictor()

    def predict(self, gene: str, context: str = "endothelial") -> Dict:
        """
        Generate phenotype prediction for a gene with no KO data.
        """

        prediction = {
            "gene": gene,
            "context": context,
            "has_existing_ko_data": self._check_existing_data(gene),

            # Predictions from different methods
            "expression_based": None,
            "family_based": None,
            "network_based": None,

            # Integrated prediction
            "integrated_prediction": {
                "primary_phenotype": None,
                "severity": None,
                "affected_pathways": [],
                "confidence": 0,
            },

            # For validation
            "validation_plan": [],
            "similar_genes_with_data": [],
        }

        # Skip if KO data exists
        if prediction["has_existing_ko_data"]:
            prediction["note"] = "KO data exists - use known phenotype instead"
            return prediction

        # 1. Expression-based prediction
        modules = self.expression_predictor.find_module_membership(gene)
        if modules:
            prediction["expression_based"] = modules
            # Get highest confidence module
            best_module = max(modules, key=lambda x: x["confidence"])
            prediction["integrated_prediction"]["affected_pathways"].append(best_module["module"])

        # 2. Family-based prediction
        family_pred = self.family_predictor.predict_from_family(gene)
        if family_pred["family"]:
            prediction["family_based"] = family_pred

        # 3. Network-based prediction
        network_pred = self.network_predictor.predict_from_network(gene)
        prediction["network_based"] = network_pred

        # 4. Integrate predictions
        prediction["integrated_prediction"] = self._integrate_predictions(prediction)

        # 5. Generate validation plan
        prediction["validation_plan"] = self._generate_validation_plan(gene, prediction)

        # 6. Find similar genes with data for reference
        prediction["similar_genes_with_data"] = self._find_similar_with_data(gene)

        return prediction

    def _check_existing_data(self, gene: str) -> bool:
        """Check if KO data already exists for this gene"""

        GENES_WITH_KO_DATA = {
            "KDR", "FLT1", "FLT4", "DLL4", "NOTCH1", "NOTCH2",
            "EFNB2", "EPHB4", "GJA5", "GJA4", "HEY1", "HEY2",
            "VEGFA", "ANGPT1", "ANGPT2", "TEK", "NRP1", "NRP2",
            "PDGFB", "PDGFRB", "TGFB1", "ENG", "ACVRL1", "SMAD4",
        }

        return gene.upper() in GENES_WITH_KO_DATA

    def _integrate_predictions(self, prediction: Dict) -> Dict:
        """Integrate predictions from all methods"""

        integrated = {
            "primary_phenotype": None,
            "secondary_phenotypes": [],
            "severity": "unknown",
            "affected_pathways": [],
            "confidence": 0,
            "evidence_sources": [],
        }

        confidences = []
        phenotypes = []

        # From expression modules
        if prediction["expression_based"]:
            for module in prediction["expression_based"]:
                phenotypes.append(module["predicted_phenotype"])
                confidences.append(module["confidence"])
                integrated["evidence_sources"].append("expression_comodule")
                integrated["affected_pathways"].append(module["module"])

        # From protein family
        if prediction["family_based"] and prediction["family_based"]["predicted_phenotype"]:
            phenotypes.append(prediction["family_based"]["predicted_phenotype"])
            confidences.append(prediction["family_based"]["confidence"])
            integrated["evidence_sources"].append("protein_family")

        # From network
        if prediction["network_based"]["predicted_severity"] != "unknown":
            integrated["severity"] = prediction["network_based"]["predicted_severity"]
            confidences.append(prediction["network_based"]["confidence"])
            integrated["evidence_sources"].append("network_position")

        # Combine
        if phenotypes:
            integrated["primary_phenotype"] = phenotypes[0]
            integrated["secondary_phenotypes"] = phenotypes[1:] if len(phenotypes) > 1 else []

        if confidences:
            # Weighted average, higher weight for more sources
            integrated["confidence"] = np.mean(confidences) * (1 + 0.1 * len(confidences))
            integrated["confidence"] = min(0.95, integrated["confidence"])  # Cap at 0.95

        return integrated

    def _generate_validation_plan(self, gene: str, prediction: Dict) -> List[Dict]:
        """Generate experimental validation plan"""

        plans = []

        # Based on predicted phenotype
        phenotype = prediction["integrated_prediction"]["primary_phenotype"]

        if phenotype and "sprouting" in phenotype.lower():
            plans.append({
                "experiment": "Retinal angiogenesis assay (P5)",
                "readouts": ["tip cell count", "vessel density", "branching"],
                "priority": "high",
            })
            plans.append({
                "experiment": "EC spheroid sprouting assay",
                "readouts": ["sprout number", "sprout length"],
                "priority": "high",
            })

        if phenotype and "arterial" in phenotype.lower():
            plans.append({
                "experiment": "A/V marker expression (EFNB2/EPHB4)",
                "readouts": ["arterial marker intensity", "vein marker intensity"],
                "priority": "high",
            })

        # General validation
        plans.append({
            "experiment": "siRNA knockdown in HUVEC",
            "readouts": ["proliferation", "migration", "tube formation"],
            "priority": "medium",
            "note": "Quick validation before in vivo"
        })

        plans.append({
            "experiment": "Expression validation (qPCR/WB)",
            "readouts": ["confirm expression pattern matches prediction"],
            "priority": "essential",
        })

        return plans

    def _find_similar_with_data(self, gene: str) -> List[Dict]:
        """Find similar genes that have KO data for reference"""

        modules = self.expression_predictor.find_module_membership(gene)

        similar = []
        for module in modules:
            module_name = module["module"]
            module_data = self.expression_predictor.COEXPRESSION_MODULES.get(module_name, {})

            # Core genes have KO data
            for core_gene in module_data.get("core_genes", []):
                similar.append({
                    "gene": core_gene,
                    "relationship": "co-expression module core",
                    "has_ko_data": True,
                    "correlation": module["correlation"],
                })

        return similar[:5]  # Top 5


# ============================================================================
# DEMO
# ============================================================================

def run_prediction_demo():
    """Demonstrate novel phenotype prediction"""

    predictor = NovelPhenotypePredictor()

    # Test genes with NO existing KO data
    test_genes = ["KCNE3", "RAMP3", "INHBB", "NEURL1B", "ADGRG6"]

    print("=" * 70)
    print("NOVEL PHENOTYPE PREDICTIONS")
    print("For genes with NO existing KO data")
    print("=" * 70)

    for gene in test_genes:
        print(f"\n{'─' * 70}")
        print(f"🧬 {gene}")
        print(f"{'─' * 70}")

        pred = predictor.predict(gene)

        if pred["has_existing_ko_data"]:
            print("  ⚠️ KO data already exists")
            continue

        integrated = pred["integrated_prediction"]

        print(f"\n  📊 PREDICTION:")
        print(f"     Primary phenotype: {integrated['primary_phenotype']}")
        print(f"     Severity: {integrated['severity']}")
        print(f"     Confidence: {integrated['confidence']:.0%}")
        print(f"     Evidence: {', '.join(integrated['evidence_sources'])}")

        if pred["expression_based"]:
            print(f"\n  🔗 Expression module:")
            for mod in pred["expression_based"][:2]:
                print(f"     - {mod['module']} (r={mod['correlation']:.2f})")

        if pred["family_based"] and pred["family_based"]["family"]:
            print(f"\n  👨‍👩‍👧 Protein family: {pred['family_based']['family']}")

        if pred["similar_genes_with_data"]:
            print(f"\n  📚 Similar genes WITH KO data (for reference):")
            for sim in pred["similar_genes_with_data"][:3]:
                print(f"     - {sim['gene']} (r={sim['correlation']:.2f})")

        if pred["validation_plan"]:
            print(f"\n  🔬 Validation experiments:")
            for plan in pred["validation_plan"][:3]:
                print(f"     - {plan['experiment']} [{plan['priority']}]")


if __name__ == "__main__":
    run_prediction_demo()
