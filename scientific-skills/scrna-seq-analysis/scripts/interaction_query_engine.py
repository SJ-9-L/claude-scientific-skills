"""
Interaction Query Engine
========================
Real-time database integration for comprehensive molecular interaction queries.
Supports any molecule type: gene, protein, ligand, receptor, TF, chemical, etc.
"""

import requests
import json
import time
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
import re


# ============================================================================
# DATA CLASSES
# ============================================================================

class MoleculeType(Enum):
    GENE = "gene"
    PROTEIN = "protein"
    LIGAND = "ligand"
    RECEPTOR = "receptor"
    TRANSCRIPTION_FACTOR = "tf"
    KINASE = "kinase"
    PHOSPHATASE = "phosphatase"
    CHEMICAL = "chemical"
    ANTIBODY = "antibody"
    CYTOKINE = "cytokine"
    CHEMOKINE = "chemokine"
    GROWTH_FACTOR = "growth_factor"
    ENZYME = "enzyme"


class InteractionType(Enum):
    ACTIVATION = "activation"
    INHIBITION = "inhibition"
    BINDING = "binding"
    PHOSPHORYLATION = "phosphorylation"
    DEPHOSPHORYLATION = "dephosphorylation"
    TRANSCRIPTION_ACTIVATION = "transcription_activation"
    TRANSCRIPTION_REPRESSION = "transcription_repression"
    UBIQUITINATION = "ubiquitination"
    COMPLEX_FORMATION = "complex_formation"
    CROSSTALK = "crosstalk"
    REGULATION = "regulation"


@dataclass
class Interaction:
    """Single molecular interaction"""
    source: str
    target: str
    interaction_type: InteractionType
    effect: str  # positive, negative, unknown
    score: float  # confidence 0-1
    mechanism: Optional[str] = None
    pubmed_ids: List[str] = field(default_factory=list)
    source_db: str = ""
    cell_type: Optional[str] = None
    tissue: Optional[str] = None
    species: str = "human"


# ============================================================================
# DATABASE API CLIENTS
# ============================================================================

class STRINGClient:
    """STRING Database API for protein-protein interactions"""

    BASE_URL = "https://string-db.org/api"
    VERSION = "11.5"

    # STRING action types mapping
    ACTION_MAPPING = {
        "activation": InteractionType.ACTIVATION,
        "inhibition": InteractionType.INHIBITION,
        "binding": InteractionType.BINDING,
        "ptmod": InteractionType.PHOSPHORYLATION,  # post-translational modification
        "expression": InteractionType.TRANSCRIPTION_ACTIVATION,
        "catalysis": InteractionType.ACTIVATION,
        "reaction": InteractionType.ACTIVATION,
    }

    def __init__(self, species: int = 9606):  # 9606 = human, 10090 = mouse
        self.species = species

    @lru_cache(maxsize=1000)
    def get_interactions(self, protein: str, score_threshold: float = 0.7) -> List[Interaction]:
        """Get protein interactions from STRING"""

        interactions = []

        # Get network
        url = f"{self.BASE_URL}/json/network"
        params = {
            "identifiers": protein,
            "species": self.species,
            "required_score": int(score_threshold * 1000),
            "network_type": "functional",
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            for item in data:
                interaction = Interaction(
                    source=item.get("preferredName_A", item.get("stringId_A", "")),
                    target=item.get("preferredName_B", item.get("stringId_B", "")),
                    interaction_type=InteractionType.BINDING,
                    effect="unknown",
                    score=item.get("score", 0) / 1000,
                    source_db="STRING",
                    species="human" if self.species == 9606 else "mouse",
                )
                interactions.append(interaction)

        except Exception as e:
            print(f"STRING query error: {e}")

        # Get functional associations with directionality
        url_actions = f"{self.BASE_URL}/json/actions"
        params_actions = {
            "identifiers": protein,
            "species": self.species,
        }

        try:
            response = requests.get(url_actions, params=params_actions, timeout=30)
            response.raise_for_status()
            actions = response.json()

            for action in actions:
                action_type = action.get("mode", "binding")
                int_type = self.ACTION_MAPPING.get(action_type, InteractionType.BINDING)

                # Determine effect direction
                effect = "unknown"
                if action.get("is_directional"):
                    if action.get("a_is_acting"):
                        effect = "positive" if action_type in ["activation", "expression"] else "negative"

                interaction = Interaction(
                    source=action.get("preferredName_A", ""),
                    target=action.get("preferredName_B", ""),
                    interaction_type=int_type,
                    effect=effect,
                    score=action.get("score", 0) / 1000,
                    mechanism=action.get("mode"),
                    source_db="STRING",
                )
                interactions.append(interaction)

        except Exception as e:
            print(f"STRING actions query error: {e}")

        return interactions

    def get_enrichment(self, genes: List[str]) -> Dict:
        """Get pathway enrichment for gene list"""

        url = f"{self.BASE_URL}/json/enrichment"
        params = {
            "identifiers": "\r".join(genes),
            "species": self.species,
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            print(f"STRING enrichment error: {e}")
            return {}


class ReactomeClient:
    """Reactome Database API for pathway information"""

    BASE_URL = "https://reactome.org/ContentService"

    @lru_cache(maxsize=1000)
    def get_pathways(self, gene: str, species: str = "Homo sapiens") -> List[Dict]:
        """Get pathways containing a gene"""

        pathways = []

        # Search for entity
        url = f"{self.BASE_URL}/search/query"
        params = {
            "query": gene,
            "species": species,
            "types": "Protein,Gene",
            "cluster": "true",
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            # Get pathway membership for each result
            for result in data.get("results", [])[:5]:  # Limit to top 5
                entity_id = result.get("stId")
                if entity_id:
                    pathway_url = f"{self.BASE_URL}/data/pathways/low/entity/{entity_id}"
                    pathway_response = requests.get(pathway_url, timeout=30)
                    if pathway_response.ok:
                        for pathway in pathway_response.json():
                            pathways.append({
                                "id": pathway.get("stId"),
                                "name": pathway.get("displayName"),
                                "species": pathway.get("speciesName"),
                            })

        except Exception as e:
            print(f"Reactome query error: {e}")

        return pathways

    def get_pathway_participants(self, pathway_id: str) -> List[Dict]:
        """Get all participants in a pathway"""

        url = f"{self.BASE_URL}/data/participants/{pathway_id}"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            print(f"Reactome participants error: {e}")
            return []


class KEGGClient:
    """KEGG Database API for pathway and interaction data"""

    BASE_URL = "https://rest.kegg.jp"

    # KEGG relation types
    RELATION_MAPPING = {
        "activation": InteractionType.ACTIVATION,
        "inhibition": InteractionType.INHIBITION,
        "expression": InteractionType.TRANSCRIPTION_ACTIVATION,
        "repression": InteractionType.TRANSCRIPTION_REPRESSION,
        "phosphorylation": InteractionType.PHOSPHORYLATION,
        "dephosphorylation": InteractionType.DEPHOSPHORYLATION,
        "ubiquitination": InteractionType.UBIQUITINATION,
        "binding/association": InteractionType.BINDING,
        "compound": InteractionType.BINDING,
    }

    @lru_cache(maxsize=1000)
    def get_gene_info(self, gene: str, organism: str = "hsa") -> Dict:
        """Get KEGG gene information"""

        # First, find KEGG gene ID
        url = f"{self.BASE_URL}/find/genes/{gene}"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            # Parse response
            lines = response.text.strip().split("\n")
            for line in lines:
                if organism in line:
                    kegg_id = line.split("\t")[0]
                    break
            else:
                return {}

            # Get gene details
            detail_url = f"{self.BASE_URL}/get/{kegg_id}"
            detail_response = requests.get(detail_url, timeout=30)

            if detail_response.ok:
                return self._parse_kegg_entry(detail_response.text)

        except Exception as e:
            print(f"KEGG gene info error: {e}")

        return {}

    def get_pathway_info(self, pathway_id: str) -> Dict:
        """Get KEGG pathway information"""

        url = f"{self.BASE_URL}/get/{pathway_id}/kgml"

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return self._parse_kgml(response.text)
        except Exception as e:
            print(f"KEGG pathway error: {e}")
            return {}

    def _parse_kegg_entry(self, text: str) -> Dict:
        """Parse KEGG flat file format"""
        result = {}
        current_field = None

        for line in text.split("\n"):
            if line.startswith(" "):
                if current_field:
                    result[current_field] += " " + line.strip()
            else:
                parts = line.split(None, 1)
                if len(parts) >= 2:
                    current_field = parts[0]
                    result[current_field] = parts[1]

        return result

    def _parse_kgml(self, xml_text: str) -> Dict:
        """Parse KGML XML format (simplified)"""
        # Would use xml.etree for full implementation
        return {"raw": xml_text}


class UniProtClient:
    """UniProt API for protein information"""

    BASE_URL = "https://rest.uniprot.org/uniprotkb"

    @lru_cache(maxsize=1000)
    def get_protein_info(self, gene_or_id: str, organism: str = "human") -> Dict:
        """Get comprehensive protein information"""

        # Search UniProt
        organism_id = "9606" if organism == "human" else "10090"
        query = f"(gene:{gene_or_id}) AND (organism_id:{organism_id})"

        url = f"{self.BASE_URL}/search"
        params = {
            "query": query,
            "format": "json",
            "size": 1,
        }

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            if data.get("results"):
                protein = data["results"][0]
                return {
                    "accession": protein.get("primaryAccession"),
                    "name": protein.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                    "gene": protein.get("genes", [{}])[0].get("geneName", {}).get("value"),
                    "function": self._extract_function(protein),
                    "subcellular_location": self._extract_location(protein),
                    "domains": self._extract_domains(protein),
                    "go_terms": self._extract_go(protein),
                    "interactions": self._extract_interactions(protein),
                    "pathways": self._extract_pathways(protein),
                }

        except Exception as e:
            print(f"UniProt query error: {e}")

        return {}

    def _extract_function(self, protein: Dict) -> str:
        comments = protein.get("comments", [])
        for comment in comments:
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    return texts[0].get("value", "")
        return ""

    def _extract_location(self, protein: Dict) -> List[str]:
        locations = []
        comments = protein.get("comments", [])
        for comment in comments:
            if comment.get("commentType") == "SUBCELLULAR LOCATION":
                for loc in comment.get("subcellularLocations", []):
                    location = loc.get("location", {}).get("value")
                    if location:
                        locations.append(location)
        return locations

    def _extract_domains(self, protein: Dict) -> List[Dict]:
        domains = []
        features = protein.get("features", [])
        for feature in features:
            if feature.get("type") == "Domain":
                domains.append({
                    "name": feature.get("description"),
                    "start": feature.get("location", {}).get("start", {}).get("value"),
                    "end": feature.get("location", {}).get("end", {}).get("value"),
                })
        return domains

    def _extract_go(self, protein: Dict) -> List[Dict]:
        go_terms = []
        references = protein.get("uniProtKBCrossReferences", [])
        for ref in references:
            if ref.get("database") == "GO":
                go_terms.append({
                    "id": ref.get("id"),
                    "term": ref.get("properties", [{}])[0].get("value", "") if ref.get("properties") else "",
                })
        return go_terms

    def _extract_interactions(self, protein: Dict) -> List[str]:
        interactors = []
        comments = protein.get("comments", [])
        for comment in comments:
            if comment.get("commentType") == "INTERACTION":
                for interaction in comment.get("interactions", []):
                    partner = interaction.get("interactantTwo", {}).get("geneName")
                    if partner:
                        interactors.append(partner)
        return interactors

    def _extract_pathways(self, protein: Dict) -> List[str]:
        pathways = []
        references = protein.get("uniProtKBCrossReferences", [])
        for ref in references:
            if ref.get("database") == "Reactome":
                props = ref.get("properties", [])
                if props:
                    pathways.append(props[0].get("value", ""))
        return pathways


class SIGNORClient:
    """SIGNOR Database for signaling interactions"""

    BASE_URL = "https://signor.uniroma2.it/api"

    # SIGNOR effect mapping
    EFFECT_MAPPING = {
        "up-regulates": ("positive", InteractionType.ACTIVATION),
        "up-regulates activity": ("positive", InteractionType.ACTIVATION),
        "up-regulates quantity": ("positive", InteractionType.TRANSCRIPTION_ACTIVATION),
        "down-regulates": ("negative", InteractionType.INHIBITION),
        "down-regulates activity": ("negative", InteractionType.INHIBITION),
        "down-regulates quantity": ("negative", InteractionType.TRANSCRIPTION_REPRESSION),
        "form complex": ("neutral", InteractionType.COMPLEX_FORMATION),
    }

    @lru_cache(maxsize=1000)
    def get_interactions(self, gene: str) -> List[Interaction]:
        """Get signaling interactions from SIGNOR"""

        interactions = []

        # Search for entity
        url = f"{self.BASE_URL}/getInteractions"
        params = {"entity": gene}

        try:
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            for item in data:
                effect_info = self.EFFECT_MAPPING.get(
                    item.get("effect", ""),
                    ("unknown", InteractionType.REGULATION)
                )

                interaction = Interaction(
                    source=item.get("entityA", ""),
                    target=item.get("entityB", ""),
                    interaction_type=effect_info[1],
                    effect=effect_info[0],
                    score=0.9 if item.get("directInteraction") else 0.7,
                    mechanism=item.get("mechanism"),
                    pubmed_ids=[item.get("pmid")] if item.get("pmid") else [],
                    source_db="SIGNOR",
                    cell_type=item.get("cellType"),
                    tissue=item.get("tissue"),
                )
                interactions.append(interaction)

        except Exception as e:
            print(f"SIGNOR query error: {e}")

        return interactions


class ChEMBLClient:
    """ChEMBL Database for chemical-target interactions"""

    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"

    @lru_cache(maxsize=1000)
    def get_target_compounds(self, gene: str) -> List[Dict]:
        """Get compounds targeting a gene/protein"""

        compounds = []

        # First get target ID
        target_url = f"{self.BASE_URL}/target/search.json"
        params = {"q": gene, "limit": 5}

        try:
            response = requests.get(target_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()

            targets = data.get("targets", [])
            for target in targets:
                target_chembl_id = target.get("target_chembl_id")

                if target_chembl_id:
                    # Get activities for this target
                    activity_url = f"{self.BASE_URL}/activity.json"
                    activity_params = {
                        "target_chembl_id": target_chembl_id,
                        "limit": 50,
                        "pchembl_value__gte": 6,  # Active compounds only
                    }

                    activity_response = requests.get(activity_url, params=activity_params, timeout=30)
                    if activity_response.ok:
                        activities = activity_response.json().get("activities", [])

                        for act in activities:
                            compounds.append({
                                "compound_id": act.get("molecule_chembl_id"),
                                "compound_name": act.get("molecule_pref_name"),
                                "activity_type": act.get("standard_type"),
                                "activity_value": act.get("standard_value"),
                                "activity_unit": act.get("standard_units"),
                                "pchembl": act.get("pchembl_value"),
                                "assay_type": act.get("assay_type"),
                            })

        except Exception as e:
            print(f"ChEMBL query error: {e}")

        return compounds


# ============================================================================
# MAIN QUERY ENGINE
# ============================================================================

class MolecularInteractionEngine:
    """
    Main engine for comprehensive molecular interaction queries.
    Aggregates data from multiple databases.
    """

    def __init__(self, species: str = "human"):
        self.species = species
        species_id = 9606 if species == "human" else 10090

        # Initialize clients
        self.string = STRINGClient(species=species_id)
        self.reactome = ReactomeClient()
        self.kegg = KEGGClient()
        self.uniprot = UniProtClient()
        self.signor = SIGNORClient()
        self.chembl = ChEMBLClient()

    def query(self, molecule: str, molecule_type: MoleculeType = None) -> Dict:
        """
        Comprehensive query for any molecule.

        Returns ALL known interactions:
        - Upstream regulators (what activates/inhibits this)
        - Downstream targets (what this activates/inhibits)
        - Physical interactions (binding partners)
        - Pathway context
        - Crosstalk
        - Chemical modulators
        """

        print(f"Querying all databases for: {molecule}")

        result = {
            "query": molecule,
            "molecule_type": molecule_type.value if molecule_type else "unknown",
            "species": self.species,

            # Basic info
            "basic_info": {},

            # UPSTREAM (what regulates this molecule)
            "upstream": {
                "activators": [],
                "inhibitors": [],
                "transcriptional_regulators": [],
            },

            # DOWNSTREAM (what this molecule regulates)
            "downstream": {
                "activation_targets": [],
                "inhibition_targets": [],
                "transcriptional_targets": [],
            },

            # PHYSICAL INTERACTIONS
            "physical_interactions": {
                "binding_partners": [],
                "complex_members": [],
            },

            # PATHWAY CONTEXT
            "pathways": [],
            "crosstalk": [],

            # CHEMICAL MODULATORS
            "chemical_modulators": {
                "inhibitors": [],
                "activators": [],
            },

            # RAW DATA from each source
            "_sources": {},
        }

        # 1. UniProt - Basic protein info
        try:
            uniprot_data = self.uniprot.get_protein_info(molecule, self.species)
            result["basic_info"] = uniprot_data
            result["_sources"]["uniprot"] = uniprot_data
        except Exception as e:
            print(f"UniProt error: {e}")

        # 2. STRING - Protein interactions
        try:
            string_interactions = self.string.get_interactions(molecule)
            result["_sources"]["string"] = [i.__dict__ for i in string_interactions]

            for interaction in string_interactions:
                self._categorize_interaction(result, molecule, interaction)
        except Exception as e:
            print(f"STRING error: {e}")

        # 3. SIGNOR - Signaling interactions (directional)
        try:
            signor_interactions = self.signor.get_interactions(molecule)
            result["_sources"]["signor"] = [i.__dict__ for i in signor_interactions]

            for interaction in signor_interactions:
                self._categorize_interaction(result, molecule, interaction)
        except Exception as e:
            print(f"SIGNOR error: {e}")

        # 4. Reactome - Pathways
        try:
            pathways = self.reactome.get_pathways(molecule, "Homo sapiens" if self.species == "human" else "Mus musculus")
            result["pathways"] = pathways
            result["_sources"]["reactome"] = pathways
        except Exception as e:
            print(f"Reactome error: {e}")

        # 5. ChEMBL - Chemical modulators
        try:
            compounds = self.chembl.get_target_compounds(molecule)
            result["_sources"]["chembl"] = compounds

            for compound in compounds:
                if compound.get("activity_type") in ["IC50", "Ki", "Kd"]:
                    result["chemical_modulators"]["inhibitors"].append(compound)
                elif compound.get("activity_type") in ["EC50", "AC50"]:
                    result["chemical_modulators"]["activators"].append(compound)
        except Exception as e:
            print(f"ChEMBL error: {e}")

        # 6. Identify crosstalk points
        result["crosstalk"] = self._identify_crosstalk(result)

        # 7. Deduplicate and score
        result = self._deduplicate_results(result)

        return result

    def _categorize_interaction(self, result: Dict, query: str, interaction: Interaction):
        """Categorize interaction as upstream, downstream, or physical"""

        if interaction.source.upper() == query.upper():
            # Query molecule is the source -> downstream interaction
            if interaction.interaction_type == InteractionType.ACTIVATION:
                result["downstream"]["activation_targets"].append({
                    "target": interaction.target,
                    "mechanism": interaction.mechanism,
                    "score": interaction.score,
                    "evidence": interaction.pubmed_ids,
                    "source_db": interaction.source_db,
                })
            elif interaction.interaction_type == InteractionType.INHIBITION:
                result["downstream"]["inhibition_targets"].append({
                    "target": interaction.target,
                    "mechanism": interaction.mechanism,
                    "score": interaction.score,
                    "evidence": interaction.pubmed_ids,
                    "source_db": interaction.source_db,
                })
            elif interaction.interaction_type in [InteractionType.TRANSCRIPTION_ACTIVATION, InteractionType.TRANSCRIPTION_REPRESSION]:
                result["downstream"]["transcriptional_targets"].append({
                    "target": interaction.target,
                    "effect": "activation" if interaction.interaction_type == InteractionType.TRANSCRIPTION_ACTIVATION else "repression",
                    "score": interaction.score,
                    "evidence": interaction.pubmed_ids,
                    "source_db": interaction.source_db,
                })
            else:
                result["physical_interactions"]["binding_partners"].append({
                    "partner": interaction.target,
                    "type": interaction.interaction_type.value,
                    "score": interaction.score,
                    "source_db": interaction.source_db,
                })

        elif interaction.target.upper() == query.upper():
            # Query molecule is the target -> upstream interaction
            if interaction.interaction_type == InteractionType.ACTIVATION:
                result["upstream"]["activators"].append({
                    "regulator": interaction.source,
                    "mechanism": interaction.mechanism,
                    "score": interaction.score,
                    "evidence": interaction.pubmed_ids,
                    "source_db": interaction.source_db,
                })
            elif interaction.interaction_type == InteractionType.INHIBITION:
                result["upstream"]["inhibitors"].append({
                    "regulator": interaction.source,
                    "mechanism": interaction.mechanism,
                    "score": interaction.score,
                    "evidence": interaction.pubmed_ids,
                    "source_db": interaction.source_db,
                })
            else:
                result["physical_interactions"]["binding_partners"].append({
                    "partner": interaction.source,
                    "type": interaction.interaction_type.value,
                    "score": interaction.score,
                    "source_db": interaction.source_db,
                })

    def _identify_crosstalk(self, result: Dict) -> List[Dict]:
        """Identify pathway crosstalk points"""

        crosstalk = []
        pathways = result.get("pathways", [])

        # Known crosstalk patterns
        CROSSTALK_PATTERNS = {
            "RAS": ["MAPK", "PI3K", "RAL"],
            "AKT": ["mTOR", "FOXO", "GSK3", "NFκB"],
            "STAT3": ["NFκB", "HIF1A", "SRC"],
            "ERK": ["RSK", "ELK1", "MYC"],
            "GSK3B": ["Wnt", "mTOR", "NFκB"],
            "NOTCH": ["Wnt", "TGFβ", "Hippo"],
            "YAP": ["Wnt", "TGFβ", "Hippo"],
        }

        # Check if query molecule is a known crosstalk point
        query = result["query"].upper()
        if query in CROSSTALK_PATTERNS:
            for target_pathway in CROSSTALK_PATTERNS[query]:
                crosstalk.append({
                    "from": query,
                    "to_pathway": target_pathway,
                    "type": "direct_crosstalk",
                })

        return crosstalk

    def _deduplicate_results(self, result: Dict) -> Dict:
        """Remove duplicates and rank by confidence score"""

        for category in ["upstream", "downstream"]:
            for subcategory in result[category]:
                items = result[category][subcategory]
                # Deduplicate by target/regulator name
                seen = set()
                unique = []
                for item in items:
                    key = item.get("target") or item.get("regulator") or item.get("partner")
                    if key and key not in seen:
                        seen.add(key)
                        unique.append(item)
                # Sort by score
                result[category][subcategory] = sorted(unique, key=lambda x: -x.get("score", 0))

        return result

    def get_pathway_interactions(self, pathway_name: str) -> Dict:
        """Get all interactions within a specific pathway"""

        # This would query Reactome/KEGG for pathway members
        # Then get all interactions between pathway members
        pass

    def compare_molecules(self, molecules: List[str]) -> Dict:
        """
        Compare interactions across multiple molecules.
        Useful for KO comparison experiments.
        """

        comparison = {
            "molecules": molecules,
            "shared_targets": [],
            "unique_targets": {},
            "shared_regulators": [],
            "unique_regulators": {},
            "pathway_overlap": [],
        }

        all_results = {}
        for mol in molecules:
            all_results[mol] = self.query(mol)

        # Find shared vs unique targets
        all_targets = {mol: set() for mol in molecules}
        for mol in molecules:
            for target in all_results[mol]["downstream"]["activation_targets"]:
                all_targets[mol].add(target["target"])
            for target in all_results[mol]["downstream"]["inhibition_targets"]:
                all_targets[mol].add(target["target"])

        # Intersection
        if len(molecules) > 1:
            shared = all_targets[molecules[0]]
            for mol in molecules[1:]:
                shared = shared.intersection(all_targets[mol])
            comparison["shared_targets"] = list(shared)

            # Unique to each
            for mol in molecules:
                unique = all_targets[mol] - shared
                comparison["unique_targets"][mol] = list(unique)

        return comparison


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if __name__ == "__main__":
    # Initialize engine
    engine = MolecularInteractionEngine(species="human")

    # Query any molecule
    result = engine.query("VEGFR2", MoleculeType.RECEPTOR)

    print("\n=== VEGFR2 Interaction Summary ===")
    print(f"\nUpstream Activators: {len(result['upstream']['activators'])}")
    for act in result['upstream']['activators'][:5]:
        print(f"  - {act['regulator']} (score: {act['score']:.2f})")

    print(f"\nUpstream Inhibitors: {len(result['upstream']['inhibitors'])}")
    for inh in result['upstream']['inhibitors'][:5]:
        print(f"  - {inh['regulator']} (score: {inh['score']:.2f})")

    print(f"\nDownstream Targets: {len(result['downstream']['activation_targets'])}")
    for tgt in result['downstream']['activation_targets'][:5]:
        print(f"  - {tgt['target']} (score: {tgt['score']:.2f})")

    print(f"\nPathways: {len(result['pathways'])}")
    for pathway in result['pathways'][:5]:
        print(f"  - {pathway['name']}")

    print(f"\nChemical Inhibitors: {len(result['chemical_modulators']['inhibitors'])}")
    for chem in result['chemical_modulators']['inhibitors'][:5]:
        print(f"  - {chem.get('compound_name', 'Unknown')} ({chem.get('activity_type')}: {chem.get('activity_value')} {chem.get('activity_unit', '')})")
