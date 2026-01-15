"""
KO Phenotype Predictor
======================
Predicts knockout phenotypes based on molecular interaction network analysis.
Uses accumulated interaction data to infer phenotypic consequences.
"""

import requests
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from functools import lru_cache
import json


# ============================================================================
# PHENOTYPE DATABASE CLIENTS
# ============================================================================

class MGIClient:
    """Mouse Genome Informatics database for phenotype data"""

    BASE_URL = "http://www.informatics.jax.org/api"

    # MP (Mammalian Phenotype) ontology categories
    PHENOTYPE_CATEGORIES = {
        "MP:0005385": "cardiovascular",
        "MP:0005386": "behavior/neurological",
        "MP:0005387": "immune",
        "MP:0005388": "respiratory",
        "MP:0005389": "reproductive",
        "MP:0005390": "skeleton",
        "MP:0005391": "vision/eye",
        "MP:0005394": "taste/olfaction",
        "MP:0005395": "hearing/vestibular",
        "MP:0005397": "hematopoietic",
        "MP:0005398": "mortality/aging",
        "MP:0005376": "homeostasis/metabolism",
        "MP:0005377": "growth/size",
        "MP:0005378": "liver/biliary",
        "MP:0005379": "endocrine/exocrine",
        "MP:0005380": "embryo",
        "MP:0005381": "digestive/alimentary",
        "MP:0005382": "craniofacial",
        "MP:0005384": "cellular",
    }

    @lru_cache(maxsize=500)
    def get_gene_phenotypes(self, gene_symbol: str) -> Dict:
        """Get all phenotypes for a gene from MGI"""

        phenotypes = {
            "gene": gene_symbol,
            "alleles": [],
            "phenotypes": [],
            "lethal": False,
            "viable": True,
        }

        # MGI uses a different API structure
        # This is a simplified representation
        try:
            # Search for gene
            search_url = f"http://www.informatics.jax.org/searchtool/Search.do"
            params = {
                "query": gene_symbol,
                "submit": "Quick Search",
            }

            # Note: MGI doesn't have a clean REST API
            # In production, you'd use MouseMine API or download files
            # Here's MouseMine approach:
            mousemine_url = "https://www.mousemine.org/mousemine/service/query/results"

            query = f'''
            <query model="genomic" view="Gene.symbol Gene.primaryIdentifier
                   Gene.alleles.symbol Gene.alleles.alleleType
                   Gene.alleles.phenotypeSummary">
                <constraint path="Gene.symbol" op="=" value="{gene_symbol}"/>
                <constraint path="Gene.organism.shortName" op="=" value="M. musculus"/>
            </query>
            '''

            response = requests.post(
                mousemine_url,
                data={'query': query, 'format': 'json'},
                timeout=30
            )

            if response.ok:
                data = response.json()
                results = data.get("results", [])

                for result in results:
                    allele_info = {
                        "symbol": result[2] if len(result) > 2 else None,
                        "type": result[3] if len(result) > 3 else None,
                        "phenotype_summary": result[4] if len(result) > 4 else None,
                    }
                    if allele_info["symbol"]:
                        phenotypes["alleles"].append(allele_info)

                        # Check for lethality
                        if allele_info.get("phenotype_summary"):
                            summary = allele_info["phenotype_summary"].lower()
                            if "lethal" in summary or "lethality" in summary:
                                phenotypes["lethal"] = True
                                phenotypes["viable"] = False

        except Exception as e:
            print(f"MGI query error: {e}")

        return phenotypes


class IMPCClient:
    """International Mouse Phenotyping Consortium database"""

    BASE_URL = "https://www.ebi.ac.uk/mi/impc/solr"

    @lru_cache(maxsize=500)
    def get_gene_phenotypes(self, gene_symbol: str) -> Dict:
        """Get standardized phenotyping data from IMPC"""

        result = {
            "gene": gene_symbol,
            "phenotypes": [],
            "significant_parameters": [],
            "viability": None,
            "fertility": None,
        }

        try:
            # Get gene info
            gene_url = f"{self.BASE_URL}/gene/select"
            params = {
                'q': f'marker_symbol:{gene_symbol}',
                'wt': 'json',
                'rows': 1
            }
            gene_response = requests.get(gene_url, params=params, timeout=30)

            if not gene_response.ok:
                return result

            gene_data = gene_response.json()
            if gene_data['response']['numFound'] == 0:
                return result

            gene_info = gene_data['response']['docs'][0]
            mgi_id = gene_info.get('mgi_accession_id')

            # Get phenotype data
            phenotype_url = f"{self.BASE_URL}/genotype-phenotype/select"
            params = {
                'q': f'marker_accession_id:"{mgi_id}"',
                'wt': 'json',
                'rows': 1000
            }
            phenotype_response = requests.get(phenotype_url, params=params, timeout=30)

            if phenotype_response.ok:
                phenotype_data = phenotype_response.json()

                for doc in phenotype_data['response']['docs']:
                    phenotype = {
                        "mp_term": doc.get("mp_term_name"),
                        "mp_id": doc.get("mp_term_id"),
                        "parameter": doc.get("parameter_name"),
                        "procedure": doc.get("procedure_name"),
                        "p_value": doc.get("p_value"),
                        "effect_size": doc.get("effect_size"),
                        "sex": doc.get("sex"),
                        "zygosity": doc.get("zygosity"),
                    }
                    result["phenotypes"].append(phenotype)

                    # Track significant parameters
                    if doc.get("p_value") and float(doc.get("p_value", 1)) < 0.0001:
                        result["significant_parameters"].append(doc.get("parameter_name"))

            # Get viability/fertility
            result["viability"] = gene_info.get("viability")
            result["fertility"] = gene_info.get("fertility")

        except Exception as e:
            print(f"IMPC query error: {e}")

        return result


class HPOClient:
    """Human Phenotype Ontology for human disease phenotypes"""

    BASE_URL = "https://hpo.jax.org/api/hpo"

    @lru_cache(maxsize=500)
    def get_gene_phenotypes(self, gene_symbol: str) -> Dict:
        """Get human phenotypes associated with a gene"""

        result = {
            "gene": gene_symbol,
            "phenotypes": [],
            "diseases": [],
        }

        try:
            # Search for gene
            search_url = f"{self.BASE_URL}/search"
            params = {"q": gene_symbol, "max": 10, "category": "genes"}

            response = requests.get(search_url, params=params, timeout=30)

            if response.ok:
                data = response.json()

                # Get associated phenotypes
                for term in data.get("terms", []):
                    if term.get("id", "").startswith("HP:"):
                        result["phenotypes"].append({
                            "hpo_id": term.get("id"),
                            "name": term.get("name"),
                        })

                # Get associated diseases
                for disease in data.get("diseases", []):
                    result["diseases"].append({
                        "disease_id": disease.get("id"),
                        "name": disease.get("name"),
                    })

        except Exception as e:
            print(f"HPO query error: {e}")

        return result


# ============================================================================
# PHENOTYPE PREDICTOR
# ============================================================================

class KOPhenotypePredictor:
    """
    Predicts knockout phenotypes using:
    1. Known phenotypes from MGI/IMPC
    2. Network analysis (downstream effects)
    3. Pathway membership
    4. Homolog phenotypes
    """

    def __init__(self):
        self.mgi = MGIClient()
        self.impc = IMPCClient()
        self.hpo = HPOClient()

        # Known essential genes and their phenotypes
        self.ESSENTIAL_GENE_CATEGORIES = self._load_essential_genes()

    def _load_essential_genes(self) -> Dict:
        """Load curated essential gene data"""
        return {
            # Core signaling - embryonic lethal
            "signaling_core": {
                "genes": ["KRAS", "HRAS", "NRAS", "RAF1", "BRAF", "MAP2K1", "MAPK1", "MAPK3",
                         "PIK3CA", "PIK3CB", "AKT1", "MTOR", "PTEN", "TP53"],
                "ko_outcome": "embryonic_lethal_or_severe",
            },

            # RTK signaling
            "rtk_critical": {
                "genes": ["EGFR", "ERBB2", "KDR", "FLT1", "FLT4", "PDGFRA", "PDGFRB", "FGFR1", "FGFR2"],
                "ko_outcome": "developmental_defects",
            },

            # Notch pathway
            "notch_core": {
                "genes": ["NOTCH1", "NOTCH2", "DLL1", "DLL4", "JAG1", "RBPJ", "HES1"],
                "ko_outcome": "embryonic_lethal_vascular",
            },

            # Wnt pathway
            "wnt_core": {
                "genes": ["CTNNB1", "APC", "GSK3B", "AXIN1", "LRP5", "LRP6"],
                "ko_outcome": "embryonic_lethal_axis",
            },

            # TGF-β pathway
            "tgfb_core": {
                "genes": ["TGFB1", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4"],
                "ko_outcome": "embryonic_lethal_or_immune",
            },

            # Cell cycle
            "cell_cycle": {
                "genes": ["CDK1", "CDK2", "CCND1", "CCNE1", "RB1", "E2F1"],
                "ko_outcome": "embryonic_lethal_proliferation",
            },

            # Apoptosis
            "apoptosis": {
                "genes": ["BCL2", "BAX", "BAK1", "CASP3", "CASP8", "CASP9"],
                "ko_outcome": "variable",
            },

            # Cytokines - immune phenotypes
            "cytokines_major": {
                "genes": ["TNF", "IL1B", "IL6", "IL10", "IFNG", "IL2", "IL4", "IL17A"],
                "ko_outcome": "immune_defects",
            },
        }

    def predict(self, gene: str, species: str = "mouse",
                context: str = None,
                include_network_effects: bool = True) -> Dict:
        """
        Predict phenotypes for gene knockout.

        Args:
            gene: Gene symbol
            species: mouse or human
            context: Optional context (e.g., "angiogenesis", "immune", "cancer")
            include_network_effects: Whether to predict downstream cascade effects

        Returns:
            Comprehensive phenotype prediction
        """

        prediction = {
            "gene": gene,
            "species": species,
            "context": context,

            # Known phenotypes from databases
            "known_phenotypes": {
                "mgi": {},
                "impc": {},
                "hpo": {},
            },

            # Predicted phenotypes
            "predicted_phenotypes": [],

            # Viability prediction
            "viability": {
                "predicted": "viable",
                "confidence": 0.5,
                "evidence": [],
            },

            # Severity assessment
            "severity": {
                "score": 0,  # 0-10
                "category": "unknown",  # mild, moderate, severe, lethal
            },

            # Affected systems
            "affected_systems": [],

            # Compensatory mechanisms
            "potential_compensation": [],

            # Network effects (if requested)
            "downstream_cascade": [],

            # Experimental recommendations
            "recommended_experiments": [],

            # Confidence metrics
            "confidence": {
                "overall": 0.5,
                "data_sources": 0,
                "network_support": 0,
            },
        }

        # 1. Get known phenotypes
        print(f"Fetching known phenotypes for {gene}...")

        if species == "mouse":
            prediction["known_phenotypes"]["mgi"] = self.mgi.get_gene_phenotypes(gene)
            prediction["known_phenotypes"]["impc"] = self.impc.get_gene_phenotypes(gene)

        prediction["known_phenotypes"]["hpo"] = self.hpo.get_gene_phenotypes(gene)

        # 2. Check essential gene categories
        essential_category = self._check_essential_category(gene)
        if essential_category:
            prediction["viability"]["predicted"] = essential_category["ko_outcome"]
            prediction["viability"]["evidence"].append(f"Member of {essential_category['category']} genes")
            prediction["viability"]["confidence"] = 0.8

        # 3. Analyze IMPC data for phenotype prediction
        impc_data = prediction["known_phenotypes"]["impc"]
        if impc_data.get("phenotypes"):
            prediction = self._analyze_impc_phenotypes(prediction, impc_data)

        # 4. Predict affected systems
        prediction["affected_systems"] = self._predict_affected_systems(gene, prediction)

        # 5. Network-based predictions
        if include_network_effects:
            prediction["downstream_cascade"] = self._predict_network_effects(gene)

        # 6. Find compensatory mechanisms
        prediction["potential_compensation"] = self._find_compensation(gene)

        # 7. Calculate severity
        prediction["severity"] = self._calculate_severity(prediction)

        # 8. Generate experiment recommendations
        prediction["recommended_experiments"] = self._recommend_experiments(gene, prediction)

        # 9. Calculate overall confidence
        prediction["confidence"] = self._calculate_confidence(prediction)

        return prediction

    def _check_essential_category(self, gene: str) -> Optional[Dict]:
        """Check if gene is in a known essential category"""
        gene_upper = gene.upper()
        for category, data in self.ESSENTIAL_GENE_CATEGORIES.items():
            if gene_upper in [g.upper() for g in data["genes"]]:
                return {"category": category, "ko_outcome": data["ko_outcome"]}
        return None

    def _analyze_impc_phenotypes(self, prediction: Dict, impc_data: Dict) -> Dict:
        """Analyze IMPC phenotype data"""

        # Count phenotypes by system
        system_counts = {}
        severe_phenotypes = []

        for pheno in impc_data.get("phenotypes", []):
            mp_term = pheno.get("mp_term", "")
            p_value = pheno.get("p_value", 1)

            # Categorize by system
            for mp_id, system in MGIClient.PHENOTYPE_CATEGORIES.items():
                if system.lower() in mp_term.lower():
                    system_counts[system] = system_counts.get(system, 0) + 1

            # Track severe phenotypes (highly significant)
            if p_value and float(p_value) < 1e-6:
                severe_phenotypes.append(mp_term)

        # Add to prediction
        if system_counts:
            top_systems = sorted(system_counts.items(), key=lambda x: -x[1])[:5]
            prediction["affected_systems"] = [s[0] for s in top_systems]

        if severe_phenotypes:
            for pheno in severe_phenotypes[:10]:
                prediction["predicted_phenotypes"].append({
                    "phenotype": pheno,
                    "source": "IMPC",
                    "confidence": 0.95,
                    "evidence": "highly_significant",
                })

        # Check viability
        if impc_data.get("viability"):
            prediction["viability"]["predicted"] = impc_data["viability"]
            prediction["viability"]["confidence"] = 0.95
            prediction["viability"]["evidence"].append("IMPC viability screen")

        return prediction

    def _predict_affected_systems(self, gene: str, prediction: Dict) -> List[str]:
        """Predict which biological systems will be affected"""

        # Known gene-system associations
        GENE_SYSTEM_MAP = {
            # Cardiovascular
            "cardiovascular": ["VEGFA", "KDR", "FLT1", "FLT4", "DLL4", "NOTCH1", "NRP1", "NRP2",
                             "ANGPT1", "ANGPT2", "TEK", "PECAM1", "CDH5", "GATA2", "ETV2"],

            # Immune
            "immune": ["TNF", "IL1B", "IL6", "IL10", "IFNG", "CD4", "CD8A", "FOXP3",
                      "NFKB1", "RELA", "TLR4", "MYD88", "TRAF6", "JAK1", "JAK2", "STAT3"],

            # Nervous system
            "nervous": ["BDNF", "NGF", "NTRK1", "NTRK2", "GDNF", "GFRA1", "RET",
                       "SHH", "GLI1", "NOTCH1", "DLL1", "NUMB", "NEUROD1"],

            # Metabolic
            "metabolic": ["INS", "INSR", "IRS1", "AKT1", "MTOR", "PTEN", "FOXO1",
                         "PPARG", "PPARA", "HNF4A", "PDX1", "GLUT4"],

            # Skeletal
            "skeletal": ["RUNX2", "SP7", "SOX9", "COL1A1", "COL2A1", "BMP2", "BMP4",
                        "TGFB1", "SMAD3", "WNT3A", "CTNNB1"],

            # Reproductive
            "reproductive": ["ESR1", "ESR2", "AR", "FSHR", "LHR", "GNRHR", "AMH",
                           "INHBA", "BMP15", "GDF9"],
        }

        affected = []
        gene_upper = gene.upper()

        for system, genes in GENE_SYSTEM_MAP.items():
            if gene_upper in [g.upper() for g in genes]:
                affected.append(system)

        # Add systems from IMPC/MGI data
        for system in prediction.get("affected_systems", []):
            if system not in affected:
                affected.append(system)

        return affected

    def _predict_network_effects(self, gene: str) -> List[Dict]:
        """Predict downstream cascade effects from KO"""

        # This would integrate with the MolecularInteractionEngine
        # For now, return known pathway effects

        PATHWAY_CASCADES = {
            "VEGFR2": [
                {"target": "PLCγ", "effect": "reduced_activity", "consequence": "decreased_calcium_signaling"},
                {"target": "ERK1/2", "effect": "reduced_phosphorylation", "consequence": "decreased_proliferation"},
                {"target": "AKT", "effect": "reduced_phosphorylation", "consequence": "increased_apoptosis"},
                {"target": "eNOS", "effect": "reduced_activity", "consequence": "decreased_NO_production"},
                {"target": "SRC", "effect": "reduced_activity", "consequence": "decreased_migration"},
            ],
            "NOTCH1": [
                {"target": "HES1", "effect": "absent", "consequence": "loss_of_lateral_inhibition"},
                {"target": "HEY1", "effect": "absent", "consequence": "arterial_specification_defect"},
                {"target": "MYC", "effect": "reduced", "consequence": "decreased_proliferation"},
                {"target": "CCND1", "effect": "reduced", "consequence": "cell_cycle_defect"},
            ],
            "DLL4": [
                {"target": "NOTCH1", "effect": "no_activation", "consequence": "excess_tip_cells"},
                {"target": "HES1", "effect": "reduced", "consequence": "hypersprouting"},
                {"target": "VEGFR2", "effect": "increased", "consequence": "enhanced_VEGF_response"},
            ],
        }

        return PATHWAY_CASCADES.get(gene.upper(), [])

    def _find_compensation(self, gene: str) -> List[Dict]:
        """Find potential compensatory mechanisms"""

        # Known paralog families
        PARALOG_FAMILIES = {
            "VEGFR": ["KDR", "FLT1", "FLT4"],
            "NOTCH": ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"],
            "JAK": ["JAK1", "JAK2", "JAK3", "TYK2"],
            "STAT": ["STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6"],
            "AKT": ["AKT1", "AKT2", "AKT3"],
            "ERK": ["MAPK1", "MAPK3"],
            "SMAD": ["SMAD1", "SMAD2", "SMAD3", "SMAD5", "SMAD8"],
            "FGFR": ["FGFR1", "FGFR2", "FGFR3", "FGFR4"],
            "PDGFR": ["PDGFRA", "PDGFRB"],
            "ERBB": ["EGFR", "ERBB2", "ERBB3", "ERBB4"],
        }

        compensation = []
        gene_upper = gene.upper()

        for family, members in PARALOG_FAMILIES.items():
            if gene_upper in [m.upper() for m in members]:
                paralogs = [m for m in members if m.upper() != gene_upper]
                for paralog in paralogs:
                    compensation.append({
                        "type": "paralog",
                        "gene": paralog,
                        "family": family,
                        "likelihood": "moderate",
                        "note": f"{paralog} may partially compensate for {gene} loss",
                    })

        return compensation

    def _calculate_severity(self, prediction: Dict) -> Dict:
        """Calculate overall severity score"""

        severity = {
            "score": 5,  # Default moderate
            "category": "moderate",
            "factors": [],
        }

        # Factor 1: Viability
        viability = prediction["viability"]["predicted"]
        if "lethal" in viability.lower():
            severity["score"] = 10
            severity["category"] = "lethal"
            severity["factors"].append("embryonic_lethal")
        elif "severe" in viability.lower():
            severity["score"] = 8
            severity["category"] = "severe"
            severity["factors"].append("severe_developmental_defects")

        # Factor 2: Number of affected systems
        n_systems = len(prediction["affected_systems"])
        if n_systems >= 3:
            severity["score"] = min(10, severity["score"] + 2)
            severity["factors"].append(f"affects_{n_systems}_systems")

        # Factor 3: Number of significant phenotypes
        n_phenotypes = len(prediction["predicted_phenotypes"])
        if n_phenotypes >= 10:
            severity["score"] = min(10, severity["score"] + 1)
            severity["factors"].append(f"{n_phenotypes}_significant_phenotypes")

        # Factor 4: Compensation available
        if prediction["potential_compensation"]:
            severity["score"] = max(1, severity["score"] - 1)
            severity["factors"].append("compensation_possible")

        # Update category based on final score
        if severity["score"] <= 3:
            severity["category"] = "mild"
        elif severity["score"] <= 6:
            severity["category"] = "moderate"
        elif severity["score"] <= 8:
            severity["category"] = "severe"
        else:
            severity["category"] = "lethal"

        return severity

    def _recommend_experiments(self, gene: str, prediction: Dict) -> List[Dict]:
        """Generate experimental recommendations"""

        experiments = []

        # 1. Based on affected systems
        for system in prediction["affected_systems"][:3]:
            if system == "cardiovascular":
                experiments.append({
                    "assay": "Retinal angiogenesis (P5)",
                    "readout": ["vessel_density", "tip_cell_number", "branching"],
                    "priority": "high",
                })
                experiments.append({
                    "assay": "Aortic ring sprouting",
                    "readout": ["sprout_length", "sprout_number"],
                    "priority": "high",
                })
            elif system == "immune":
                experiments.append({
                    "assay": "Flow cytometry immune panel",
                    "readout": ["T_cell_subsets", "B_cells", "myeloid_cells"],
                    "priority": "high",
                })
                experiments.append({
                    "assay": "Cytokine profiling",
                    "readout": ["serum_cytokines", "tissue_cytokines"],
                    "priority": "medium",
                })

        # 2. Based on viability
        if "lethal" in prediction["viability"]["predicted"]:
            experiments.append({
                "assay": "Timed embryo analysis",
                "readout": ["developmental_stage", "morphology"],
                "priority": "critical",
                "note": "Determine lethal stage",
            })
            experiments.append({
                "assay": "Conditional KO (Cre-lox)",
                "readout": ["tissue_specific_phenotype"],
                "priority": "high",
                "note": "Bypass embryonic lethality",
            })

        # 3. Standard validation
        experiments.append({
            "assay": "Western blot confirmation",
            "readout": ["protein_absence", "pathway_activation"],
            "priority": "essential",
        })
        experiments.append({
            "assay": "qPCR pathway genes",
            "readout": ["downstream_target_expression"],
            "priority": "essential",
        })

        return experiments

    def _calculate_confidence(self, prediction: Dict) -> Dict:
        """Calculate confidence metrics"""

        confidence = {
            "overall": 0.5,
            "data_sources": 0,
            "factors": [],
        }

        # Factor 1: Data sources available
        sources = 0
        if prediction["known_phenotypes"]["mgi"].get("alleles"):
            sources += 1
            confidence["factors"].append("MGI_data_available")
        if prediction["known_phenotypes"]["impc"].get("phenotypes"):
            sources += 1
            confidence["factors"].append("IMPC_data_available")
        if prediction["known_phenotypes"]["hpo"].get("phenotypes"):
            sources += 1
            confidence["factors"].append("HPO_data_available")

        confidence["data_sources"] = sources

        # Calculate overall confidence
        base_confidence = 0.3 + (sources * 0.2)  # 0.3 to 0.9 based on data

        # Boost if viability data is solid
        if prediction["viability"]["confidence"] > 0.8:
            base_confidence += 0.1

        confidence["overall"] = min(0.95, base_confidence)

        return confidence

    def compare_knockouts(self, genes: List[str], species: str = "mouse") -> Dict:
        """
        Compare phenotypes across multiple gene knockouts.

        Useful for: VEGFR2 KO vs VEGFR3 KO vs DLL4 KO comparison
        """

        comparison = {
            "genes": genes,
            "individual_predictions": {},
            "phenotype_comparison": {
                "shared": [],
                "unique": {},
            },
            "severity_ranking": [],
            "viability_comparison": {},
            "recommended_experiment_groups": [],
        }

        # Get predictions for each gene
        all_phenotypes = {}
        for gene in genes:
            pred = self.predict(gene, species)
            comparison["individual_predictions"][gene] = pred

            # Collect phenotypes
            phenotype_set = set()
            for p in pred["predicted_phenotypes"]:
                phenotype_set.add(p["phenotype"])
            for p in pred["known_phenotypes"]["impc"].get("phenotypes", []):
                if p.get("mp_term"):
                    phenotype_set.add(p["mp_term"])
            all_phenotypes[gene] = phenotype_set

            # Viability
            comparison["viability_comparison"][gene] = pred["viability"]["predicted"]

        # Find shared vs unique phenotypes
        if len(genes) > 1:
            shared = all_phenotypes[genes[0]]
            for gene in genes[1:]:
                shared = shared.intersection(all_phenotypes[gene])
            comparison["phenotype_comparison"]["shared"] = list(shared)

            for gene in genes:
                unique = all_phenotypes[gene] - shared
                comparison["phenotype_comparison"]["unique"][gene] = list(unique)

        # Severity ranking
        severity_list = []
        for gene in genes:
            score = comparison["individual_predictions"][gene]["severity"]["score"]
            severity_list.append((gene, score))
        comparison["severity_ranking"] = sorted(severity_list, key=lambda x: -x[1])

        # Experiment recommendations for comparison
        comparison["recommended_experiment_groups"] = [
            {
                "design": f"Compare {' vs '.join(genes)} vs WT",
                "n_per_group": 8,
                "assays": ["phenotype_screen", "pathway_analysis", "transcriptomics"],
            }
        ]

        return comparison


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if __name__ == "__main__":
    predictor = KOPhenotypePredictor()

    # Single gene prediction
    print("=== VEGFR2 KO Prediction ===")
    vegfr2_pred = predictor.predict("KDR", species="mouse", context="angiogenesis")

    print(f"\nViability: {vegfr2_pred['viability']['predicted']}")
    print(f"Severity: {vegfr2_pred['severity']['category']} (score: {vegfr2_pred['severity']['score']})")
    print(f"Affected systems: {vegfr2_pred['affected_systems']}")

    print("\nPredicted phenotypes:")
    for p in vegfr2_pred['predicted_phenotypes'][:5]:
        print(f"  - {p['phenotype']}")

    print("\nDownstream cascade:")
    for c in vegfr2_pred['downstream_cascade'][:3]:
        print(f"  - {c['target']}: {c['effect']} -> {c['consequence']}")

    # Comparison
    print("\n\n=== KO Comparison: VEGFR2 vs VEGFR3 vs DLL4 ===")
    comparison = predictor.compare_knockouts(["KDR", "FLT4", "DLL4"])

    print(f"\nSeverity ranking:")
    for gene, score in comparison['severity_ranking']:
        print(f"  {gene}: {score}")

    print(f"\nShared phenotypes: {len(comparison['phenotype_comparison']['shared'])}")
    print(f"Unique to each:")
    for gene, unique in comparison['phenotype_comparison']['unique'].items():
        print(f"  {gene}: {len(unique)} unique phenotypes")
