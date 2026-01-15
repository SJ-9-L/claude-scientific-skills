# End-to-End Workflow Examples

## Overview

실제 연구 시나리오를 기반으로 한 완전한 워크플로우 예제. 데이터 로딩부터 publication-ready output까지 전체 과정을 단계별로 설명합니다.

---

## Workflow 1: Novel Marker Discovery in Tumor Angiogenesis

### Scenario
연구 목표: 종양 혈관신생에서 새로운 tip cell marker를 발굴하고, 이를 항암 치료 표적으로 검증

### Step 1: Data Loading & Quality Control

```python
import scanpy as sc
import pandas as pd
import numpy as np
from scripts.scrna_pipeline import ScRNAAnalyzer
from scripts.database_integrator import DatabaseIntegrator

# Initialize analyzer
analyzer = ScRNAAnalyzer()

# Load tumor endothelial cell data (예: GEO GSE189357)
adata = analyzer.load_data(
    path="tumor_endothelial.h5ad",
    # 또는 10X format
    # path="filtered_feature_bc_matrix/",
    # format="10x_mtx"
)

print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
# Expected output: Loaded 15,234 cells, 33,538 genes

# Quality control
adata = analyzer.quality_control(
    adata,
    min_genes=200,      # 최소 발현 유전자 수
    max_genes=6000,     # 최대 발현 유전자 수 (doublet 제거)
    max_mt_pct=15,      # 미토콘드리아 유전자 비율 (죽은 세포 제거)
    min_cells=3         # 최소 세포 수 (rare gene 제거)
)

# QC visualization
analyzer.plot_qc(adata, save="qc_plots.png")
```

**Expected QC Output:**
```
QC Summary:
- Initial cells: 15,234
- After filtering: 12,847 (84.3%)
- Removed: 2,387 cells (15.7%)
  - Low gene count: 1,234
  - High gene count (doublets): 456
  - High MT%: 697
- Final gene count: 18,234 genes
```

### Step 2: Preprocessing & Clustering

```python
# Normalization and feature selection
adata = analyzer.preprocess(
    adata,
    normalize_total=True,     # Library size normalization
    target_sum=1e4,           # Target counts per cell
    log_transform=True,       # Log1p transformation
    n_top_genes=3000,         # Highly variable genes
    regress_out=['total_counts', 'pct_counts_mt']  # Batch effect 제거
)

# Dimensionality reduction
adata = analyzer.reduce_dimensions(
    adata,
    n_pcs=50,                 # PCA components
    n_neighbors=15,           # UMAP neighbors
    min_dist=0.3             # UMAP minimum distance
)

# Clustering
adata = analyzer.cluster(
    adata,
    resolution=0.8,           # Leiden resolution
    # Multiple resolutions for comparison
    resolutions=[0.4, 0.6, 0.8, 1.0, 1.2]
)

# Visualize clusters
sc.pl.umap(adata, color=['leiden'], save='_clusters.png')
```

**Expected Clustering Output:**
```
Clustering Summary (resolution=0.8):
- Total clusters: 12
- Cluster sizes:
  - Cluster 0: 2,341 cells (18.2%)
  - Cluster 1: 2,156 cells (16.8%)
  - Cluster 2: 1,834 cells (14.3%)
  ...
```

### Step 3: Cell Type Annotation

```python
from references.cell_type_markers import CellTypeMarkerDB

marker_db = CellTypeMarkerDB()

# Automatic annotation based on markers
annotations = analyzer.annotate_clusters(
    adata,
    marker_db=marker_db,
    tissue_context="tumor_vasculature"
)

# Manual verification with known markers
ec_markers = {
    "pan_EC": ["PECAM1", "CDH5", "VWF", "ERG"],
    "tip_cell": ["DLL4", "CXCR4", "APLN", "ESM1", "ANGPT2"],
    "stalk_cell": ["JAG1", "NOTCH1", "HEY1", "VEGFR1"],
    "arterial": ["EFNB2", "GJA5", "DLL4"],
    "venous": ["NR2F2", "EPHB4", "ACKR1"],
    "lymphatic": ["PROX1", "LYVE1", "PDPN", "FLT4"],
}

sc.pl.dotplot(adata, ec_markers, groupby='leiden', save='_ec_markers.png')

# Assign final annotations
adata.obs['cell_type'] = annotations['predicted_type']
```

**Expected Annotation:**
```
Cell Type Annotations:
- Tip cells: Clusters 0, 7 (n=3,215)
- Stalk cells: Clusters 1, 3 (n=2,890)
- Arterial EC: Cluster 2 (n=1,834)
- Venous EC: Clusters 4, 5 (n=2,456)
- Capillary EC: Cluster 6 (n=1,523)
- Lymphatic EC: Cluster 8 (n=456)
- Pericytes: Cluster 9 (n=234)
- Other: Clusters 10, 11 (n=239)
```

### Step 4: Differential Expression & Novel Gene Discovery

```python
# Find markers for tip cells
tip_markers = analyzer.find_markers(
    adata,
    groupby='cell_type',
    groups='Tip cells',
    reference='rest',       # vs all other cells
    method='wilcoxon',
    n_genes=200,
    min_fold_change=1.5,
    max_pval=0.01
)

# Novel gene discovery
novel_candidates = analyzer.discover_novel_genes(
    adata,
    cluster='Tip cells',
    known_markers=["DLL4", "CXCR4", "APLN", "ESM1", "ANGPT2", "PDGFB"],
    reference_databases=['cellxgene', 'human_cell_atlas', 'pubmed'],
    novelty_criteria={
        'not_in_literature': True,      # PubMed에 없음
        'not_in_reference_atlas': True, # HCA에서 tip cell marker로 없음
        'high_specificity': True,       # 다른 세포 타입에서 낮은 발현
        'min_expression': 0.3,          # 최소 30% 세포에서 발현
    }
)

print(f"Found {len(novel_candidates)} novel tip cell marker candidates")
```

**Expected Novel Gene Output:**
```
Novel Tip Cell Marker Candidates (Top 10):

Gene      | Log2FC | Pct_tip | Pct_other | Specificity | PubMed_hits
----------|--------|---------|-----------|-------------|------------
GENE_A    | 3.24   | 78.5%   | 12.3%     | 0.92        | 0
GENE_B    | 2.89   | 65.2%   | 8.7%      | 0.88        | 2*
GENE_C    | 2.67   | 82.1%   | 15.4%     | 0.85        | 0
GENE_D    | 2.45   | 71.3%   | 11.2%     | 0.84        | 1*
...

* PubMed hits not related to tip cells or angiogenesis
```

### Step 5: Database Integration & Validation

```python
from scripts.database_integrator import DatabaseIntegrator

db = DatabaseIntegrator()

# Comprehensive validation for top candidate
top_candidate = "GENE_A"

validation = db.validate_target(
    gene=top_candidate,
    queries={
        # Protein information
        "uniprot": {
            "info": ["function", "subcellular_location", "domains"],
            "expected": "membrane_protein_or_secreted"  # For therapeutic target
        },

        # Expression specificity
        "gtex": {
            "query": "tissue_specificity",
            "expected": "vascular_enriched"
        },
        "hpa": {
            "query": ["protein_expression", "cell_type_specificity"],
        },

        # Pathway context
        "string": {
            "query": "interaction_partners",
            "filter": "VEGF_pathway_genes"
        },
        "reactome": {
            "query": "pathway_membership",
        },

        # Disease relevance
        "opentargets": {
            "query": ["disease_associations", "tractability"],
        },
        "cosmic": {
            "query": "cancer_associations",
        },

        # Druggability
        "chembl": {
            "query": "existing_compounds",
        },
        "drugbank": {
            "query": "approved_drugs",
        }
    }
)

# Generate validation report
report = db.generate_validation_report(validation)
print(report)
```

**Expected Validation Output:**
```
Target Validation Report: GENE_A
====================================

1. Protein Properties:
   - Type: Single-pass type I membrane protein
   - Domains: Extracellular domain (aa 1-450), TM (451-471), Cytoplasmic (472-623)
   - Location: Cell surface
   - Druggability: HIGH (extracellular target)

2. Expression Profile:
   - GTEx: Highest in blood vessels, heart, placenta
   - HPA: Endothelial-specific expression confirmed
   - Single-cell: Tip cell enriched (this study)

3. Pathway Context:
   - STRING partners: VEGFR2, DLL4, NOTCH1, NRP1 (score > 0.7)
   - Reactome: "Signaling by VEGF" (R-HSA-194138)
   - KEGG: Not annotated (NOVEL)

4. Disease Relevance:
   - Open Targets: Associated with tumor angiogenesis (score: 0.72)
   - COSMIC: Overexpressed in 45% of solid tumors

5. Therapeutic Potential:
   - Existing compounds: NONE (novel target)
   - Antibody feasibility: HIGH (extracellular domain)
   - Small molecule feasibility: MEDIUM (no known binding pocket)

CONCLUSION: Strong novel therapeutic target candidate
```

### Step 6: Pathway Deep Dive

```python
from scripts.pathway_analyzer import PathwayDeepDive

pathway = PathwayDeepDive(
    target="GENE_A",
    species="human",
    context="tumor_angiogenesis"
)

# Analyze position in VEGF signaling
position = pathway.analyze_pathway_position(
    pathway_name="VEGF signaling",
    determine_upstream_downstream=True
)

# Get experimental evidence
experimental = pathway.get_experimental_evidence(
    evidence_types=[
        "knockout_studies",
        "knockdown_studies",
        "overexpression_studies",
        "inhibitor_studies"
    ],
    species=["human", "mouse"],
    context="angiogenesis"
)

# Design validation experiments
from references.experiment_design_templates import KOComparisonExperiment

validation_exp = KOComparisonExperiment(
    target_genes=["GENE_A"],
    comparison_genes=["DLL4", "VEGFR2"],  # Known tip cell regulators
    species="mouse",
    context="tumor_angiogenesis"
)

exp_design = validation_exp.generate_design()
print(exp_design)
```

**Expected Experiment Design:**
```
Validation Experiment Design for GENE_A
========================================

Experimental Groups:
1. Wild-type (WT)
2. GENE_A KO (constitutive or EC-specific Cdh5-Cre)
3. DLL4 het (DLL4+/-, known tip cell regulator)
4. VEGFR2 EC-KO (positive control)

Phenotypes to Assess:
1. Tip cell formation
   - Method: Retina angiogenesis assay (P5)
   - Readout: Tip cell number, filopodia count

2. Tumor angiogenesis
   - Method: LLC or B16 subcutaneous tumor
   - Readout: Vessel density, tip cell markers

3. Sprouting angiogenesis
   - Method: Aortic ring assay, EC spheroid sprouting
   - Readout: Sprout number, length, branching

Expected Outcomes if GENE_A is tip cell regulator:
- GENE_A KO: Reduced tip cells, impaired sprouting
- Similar phenotype to DLL4 het or distinct?

Timeline: 8-12 weeks
Statistics: n=6-8 per group, ANOVA with post-hoc
```

### Step 7: Competitive Analysis

```python
from scripts.competitive_analysis import CompetitiveLandscape

landscape = CompetitiveLandscape(
    topic=f"{top_candidate} AND (angiogenesis OR tumor)",
    timeframe="2015-2024"
)

# Check if anyone else is working on this
competitors = landscape.search_literature(
    databases=["pubmed", "biorxiv", "patents"],
    query_expansion=True  # Include synonyms
)

# Research group analysis
groups = landscape.identify_key_groups(
    min_publications=3,
    include_industry=True
)

# Patent landscape
patents = landscape.search_patents(
    databases=["USPTO", "EPO", "WIPO"],
    query=f"{top_candidate} therapeutic antibody"
)

# Gap analysis
gaps = landscape.find_research_gaps(
    current_knowledge=validation,
    identify_opportunities=True
)

print(f"Competitors working on {top_candidate}: {len(groups)}")
print(f"Related patents: {len(patents)}")
print(f"Identified research gaps: {len(gaps)}")
```

**Expected Competitive Output:**
```
Competitive Landscape Analysis: GENE_A
======================================

Literature Search:
- PubMed: 12 papers (none on angiogenesis function)
- bioRxiv: 2 preprints
- Total mentions: 14

Key Research Groups: 3
1. [University Lab A] - Structural biology focus
2. [Pharma Company B] - General target screening
3. [University Lab C] - Expression profiling

Patent Landscape:
- Therapeutic patents: 0
- Diagnostic patents: 1 (broad, not specific)
- FREEDOM TO OPERATE: HIGH

Research Gaps Identified:
1. No functional studies in angiogenesis
2. No knockout mouse phenotype published
3. No therapeutic antibodies developed
4. Mechanism of action unknown

OPPORTUNITY SCORE: 9.2/10 (Highly novel, low competition)
```

### Step 8: Figure Generation

```python
from scripts.figure_generator import NatureReviewsFigure

# Create main figure
fig = NatureReviewsFigure(
    style="nature_reviews",
    dimensions=(180, 150),  # mm
    title="Novel Tip Cell Marker GENE_A in Tumor Angiogenesis"
)

# Panel A: UMAP showing expression
fig.add_panel(
    panel_id="A",
    type="umap",
    data=adata,
    color_by=top_candidate,
    title="GENE_A expression in tumor endothelium"
)

# Panel B: Dot plot of tip cell markers
fig.add_panel(
    panel_id="B",
    type="dotplot",
    genes=[top_candidate, "DLL4", "CXCR4", "APLN", "ESM1"],
    groupby="cell_type"
)

# Panel C: Pathway schematic
fig.add_panel(
    panel_id="C",
    type="pathway",
    pathway_data={
        "VEGF-A": {"type": "ligand", "position": "extracellular"},
        "VEGFR2": {"type": "receptor", "position": "membrane"},
        "GENE_A": {"type": "novel", "position": "membrane", "highlight": True},
        "DLL4": {"type": "ligand", "position": "membrane"},
        "NOTCH1": {"type": "receptor", "position": "membrane"},
    },
    interactions=[
        ("VEGF-A", "VEGFR2", "activation"),
        ("VEGFR2", "GENE_A", "regulation", "?"),  # Proposed
        ("GENE_A", "DLL4", "regulation", "?"),
        ("DLL4", "NOTCH1", "lateral_inhibition"),
    ]
)

# Panel D: Expression comparison
fig.add_panel(
    panel_id="D",
    type="violin",
    gene=top_candidate,
    groupby="cell_type",
    order=["Tip cells", "Stalk cells", "Other EC"]
)

# Export
fig.export("figure_1_novel_marker.svg")
fig.export("figure_1_novel_marker.pdf")
fig.export("figure_1_novel_marker.png", dpi=600)

print("Figure exported successfully")
```

### Final Output Summary

```
===========================================
WORKFLOW COMPLETE: Novel Marker Discovery
===========================================

Input: Tumor endothelial scRNA-seq (15,234 cells)

Key Findings:
1. Identified 12 endothelial subpopulations
2. Discovered 8 novel tip cell marker candidates
3. Top candidate: GENE_A
   - Tip cell specific (78.5% vs 12.3%)
   - Not previously associated with angiogenesis
   - Membrane protein (druggable)
   - Low competition landscape

Generated Outputs:
├── Data/
│   ├── processed_adata.h5ad
│   ├── cluster_markers.csv
│   └── novel_candidates.csv
├── Figures/
│   ├── qc_plots.png
│   ├── umap_clusters.png
│   ├── ec_markers_dotplot.png
│   └── figure_1_novel_marker.svg/pdf/png
├── Reports/
│   ├── validation_report.md
│   ├── competitive_analysis.md
│   └── experiment_design.md
└── Supplementary/
    ├── all_DE_genes.xlsx
    └── database_queries.json

Next Steps:
1. Experimental validation (mouse KO)
2. Antibody development
3. Mechanism studies
```

---

## Workflow 2: Drug Resistance Mechanism Analysis

### Scenario
연구 목표: Anti-PD1 치료 후 획득 내성 환자의 scRNA-seq 데이터에서 내성 메커니즘 규명

### Step 1: Paired Sample Analysis

```python
from scripts.scrna_pipeline import ScRNAAnalyzer
from references.disease_specific_templates import TumorMicroenvironmentAnalyzer

analyzer = ScRNAAnalyzer()
tme = TumorMicroenvironmentAnalyzer(cancer_type="melanoma")

# Load paired samples (pre- and post-treatment)
adata_pre = analyzer.load_data("patient_01_pre_treatment.h5ad")
adata_post = analyzer.load_data("patient_01_post_resistance.h5ad")

# Add sample labels
adata_pre.obs['timepoint'] = 'pre_treatment'
adata_post.obs['timepoint'] = 'post_resistance'

# Merge datasets
adata = analyzer.merge_samples(
    [adata_pre, adata_post],
    batch_key='timepoint',
    batch_correction='harmony'  # or 'scVI', 'combat'
)

# Joint processing
adata = analyzer.full_pipeline(
    adata,
    batch_key='timepoint',
    resolution=1.0
)
```

### Step 2: Immune Cell Analysis

```python
# Identify immune cell populations
immune_annotation = tme.annotate_immune_cells(adata)

# Compare immune infiltration
comparison = tme.compare_immune_infiltration(
    adata,
    group_by='timepoint',
    groups=['pre_treatment', 'post_resistance']
)

print("Immune Cell Changes:")
for cell_type, change in comparison['fold_changes'].items():
    direction = "↑" if change > 1 else "↓"
    print(f"  {cell_type}: {direction} {abs(change):.2f}x (p={comparison['pvalues'][cell_type]:.3f})")
```

**Expected Output:**
```
Immune Cell Changes (post-resistance vs pre-treatment):
  CD8+ T cells (cytotoxic): ↓ 2.3x (p=0.001)
  CD8+ T cells (exhausted): ↑ 3.5x (p<0.001)
  Tregs: ↑ 2.1x (p=0.003)
  M2 macrophages: ↑ 1.8x (p=0.012)
  MDSC: ↑ 2.7x (p=0.002)
  NK cells: ↓ 1.5x (p=0.045)
```

### Step 3: Resistance Mechanism Identification

```python
# Analyze resistance mechanisms
resistance = tme.analyze_resistance_mechanisms(
    adata,
    treatment="anti-PD1",
    compare_groups=('pre_treatment', 'post_resistance')
)

# Check known resistance pathways
resistance_pathways = {
    "antigen_presentation": ["B2M", "HLA-A", "HLA-B", "TAP1", "TAP2"],
    "IFN_signaling": ["JAK1", "JAK2", "STAT1", "IRF1"],
    "alternative_checkpoints": ["LAG3", "HAVCR2", "TIGIT", "VISTA"],
    "immunosuppression": ["IDO1", "ARG1", "CD73", "adenosine"],
}

for pathway, genes in resistance_pathways.items():
    de_result = analyzer.test_pathway_differential(
        adata,
        genes=genes,
        groupby='timepoint',
        groups=['pre_treatment', 'post_resistance']
    )
    print(f"{pathway}: {de_result['direction']} (p={de_result['pvalue']:.3f})")
```

**Expected Output:**
```
Resistance Pathway Analysis:

1. Antigen Presentation:
   - B2M: ↓ 3.2x in tumor cells (p<0.001)
   - HLA-A: ↓ 2.1x (p=0.002)
   - Likely mechanism: Immune evasion via MHC-I loss

2. IFN Signaling:
   - JAK1: No change
   - JAK2: ↓ 1.8x (p=0.015)
   - STAT1: ↓ 1.5x (p=0.034)
   - Possible JAK/STAT pathway defect

3. Alternative Checkpoints:
   - LAG3: ↑ 2.8x on CD8+ T cells (p<0.001)
   - TIM3: ↑ 2.3x on CD8+ T cells (p=0.002)
   - TIGIT: ↑ 1.9x on CD8+ T cells (p=0.008)
   - Upregulation of compensatory checkpoints

4. Immunosuppressive Mechanisms:
   - IDO1: ↑ 4.2x in tumor/stromal cells (p<0.001)
   - ARG1: ↑ 2.5x in MDSCs (p=0.003)
   - Enhanced immunosuppressive microenvironment

PRIMARY RESISTANCE MECHANISM:
- MHC-I downregulation (B2M loss)
- Alternative checkpoint upregulation (LAG3/TIM3)
```

### Step 4: Therapeutic Recommendations

```python
# Generate therapeutic recommendations
recommendations = tme.generate_therapeutic_recommendations(
    resistance_analysis=resistance,
    current_treatment="anti-PD1"
)

print(recommendations)
```

**Expected Output:**
```
Therapeutic Recommendations Based on Resistance Analysis
=========================================================

Primary Resistance Mechanisms Identified:
1. B2M/MHC-I downregulation
2. LAG3/TIM3 checkpoint upregulation
3. IDO1-mediated immunosuppression

Recommended Combination Strategies:

Option 1: Anti-PD1 + Anti-LAG3
- Rationale: LAG3 upregulation is the most significant alternative checkpoint
- Clinical evidence: RELATIVITY-047 trial showed benefit
- Drugs: Relatlimab (approved), Fianlimab (Phase 3)

Option 2: Anti-PD1 + Anti-TIM3
- Rationale: TIM3 co-upregulated with LAG3
- Drugs: Cobolimab, Sabatolimab (Phase 2)

Option 3: Anti-PD1 + IDO1 inhibitor
- Rationale: High IDO1 expression in TME
- Caveat: Previous failures (epacadostat), but patient selection may help
- Consider: BMS-986205, Linrodostat

Option 4: Cell therapy consideration
- Rationale: MHC-I loss may preclude T cell-based approaches
- Consider: NK cell therapy (MHC-I independent)

NOT Recommended:
- Anti-CTLA4 monotherapy (no CTLA4 upregulation seen)
- Standard chemotherapy (no benefit for immune evasion)

Biomarker Monitoring:
- Track B2M expression in cfDNA/ctDNA
- Monitor LAG3/TIM3 on circulating T cells
```

---

## Workflow 3: Cross-Species Pathway Conservation Analysis

### Scenario
연구 목표: Mouse에서 발견한 새로운 pathway를 Human에서 검증하고 therapeutic target 가능성 평가

### Complete Workflow

```python
from scripts.scrna_pipeline import ScRNAAnalyzer
from scripts.pathway_analyzer import PathwayDeepDive
from scripts.database_integrator import DatabaseIntegrator

# Step 1: Load mouse discovery data
mouse_adata = analyzer.load_data("mouse_discovery_data.h5ad")
mouse_markers = analyzer.find_markers(mouse_adata, groupby='cell_type')

# Step 2: Identify human orthologs
db = DatabaseIntegrator()
orthologs = db.get_orthologs(
    genes=mouse_markers['gene'].tolist(),
    from_species='mouse',
    to_species='human',
    confidence='high'  # 1:1 orthologs only
)

print(f"Found {len(orthologs)} high-confidence human orthologs")

# Step 3: Validate in human data
human_adata = analyzer.load_data("human_reference_data.h5ad")  # e.g., from HCA

# Check conservation of expression patterns
conservation = analyzer.check_expression_conservation(
    mouse_adata=mouse_adata,
    human_adata=human_adata,
    ortholog_mapping=orthologs,
    cell_type_mapping={
        'mouse_tip_cell': 'human_tip_cell',
        'mouse_stalk_cell': 'human_stalk_cell'
    }
)

# Step 4: Identify conserved and divergent pathways
pathway_comparison = PathwayDeepDive.compare_species(
    mouse_pathway=mouse_pathway_analysis,
    human_pathway=human_pathway_analysis,
    identify_divergence=True
)

# Step 5: Prioritize human therapeutic targets
targets = db.prioritize_therapeutic_targets(
    candidates=conserved_genes,
    criteria={
        'human_expression': True,
        'druggable': True,
        'disease_association': 'cancer',
        'safety_profile': 'favorable'
    }
)

# Final output
print("\nTherapeutic Target Prioritization:")
for rank, target in enumerate(targets[:5], 1):
    print(f"{rank}. {target['gene']}")
    print(f"   Conservation score: {target['conservation']:.2f}")
    print(f"   Druggability: {target['druggability']}")
    print(f"   Safety: {target['safety_score']:.2f}")
```

---

## Workflow 4: Multi-Sample Integration & Batch Effect Correction

### Scenario
연구 목표: 여러 환자 샘플을 통합 분석하여 공통 및 환자-특이적 발현 패턴 식별

```python
from scripts.scrna_pipeline import ScRNAAnalyzer
import scanpy as sc

analyzer = ScRNAAnalyzer()

# Load multiple patient samples
samples = {
    'patient_01': 'data/patient_01.h5ad',
    'patient_02': 'data/patient_02.h5ad',
    'patient_03': 'data/patient_03.h5ad',
    'patient_04': 'data/patient_04.h5ad',
    'patient_05': 'data/patient_05.h5ad',
}

# Process each sample
adatas = {}
for name, path in samples.items():
    adata = analyzer.load_data(path)
    adata = analyzer.quality_control(adata)
    adata.obs['patient'] = name
    adatas[name] = adata
    print(f"{name}: {adata.n_obs} cells after QC")

# Concatenate
adata = sc.concat(adatas, label='patient')
print(f"\nTotal: {adata.n_obs} cells from {len(samples)} patients")

# Batch correction options
# Option 1: Harmony (fast, good for mild batch effects)
adata = analyzer.batch_correct(
    adata,
    method='harmony',
    batch_key='patient'
)

# Option 2: scVI (deep learning, better for strong batch effects)
# adata = analyzer.batch_correct(
#     adata,
#     method='scvi',
#     batch_key='patient',
#     n_latent=30
# )

# Option 3: BBKNN (graph-based)
# adata = analyzer.batch_correct(
#     adata,
#     method='bbknn',
#     batch_key='patient'
# )

# Evaluate batch correction
batch_metrics = analyzer.evaluate_batch_correction(
    adata,
    batch_key='patient',
    metrics=['silhouette', 'kBET', 'LISI']
)

print("\nBatch Correction Quality:")
print(f"  Batch mixing (LISI): {batch_metrics['LISI']:.2f}")
print(f"  Biological conservation: {batch_metrics['bio_conservation']:.2f}")

# Identify conserved vs patient-specific clusters
clustering_analysis = analyzer.analyze_cluster_composition(
    adata,
    cluster_key='leiden',
    batch_key='patient'
)

for cluster, composition in clustering_analysis.items():
    if composition['type'] == 'conserved':
        print(f"Cluster {cluster}: Conserved across patients")
    else:
        print(f"Cluster {cluster}: Enriched in {composition['enriched_patients']}")
```

---

## Output Templates

### Publication-Ready Table Format

```python
def generate_publication_table(markers_df, top_n=20):
    """Generate publication-ready marker table."""

    table = markers_df.head(top_n)[[
        'gene', 'log2fc', 'pval_adj', 'pct_in', 'pct_out', 'specificity'
    ]].copy()

    # Format columns
    table['log2fc'] = table['log2fc'].round(2)
    table['pval_adj'] = table['pval_adj'].apply(lambda x: f"{x:.2e}")
    table['pct_in'] = (table['pct_in'] * 100).round(1).astype(str) + '%'
    table['pct_out'] = (table['pct_out'] * 100).round(1).astype(str) + '%'
    table['specificity'] = table['specificity'].round(3)

    # Rename for publication
    table.columns = [
        'Gene', 'Log2 FC', 'Adj. P-value',
        '% in cluster', '% other', 'Specificity'
    ]

    return table.to_markdown(index=False)
```

### Methods Section Template

```
Single-cell RNA sequencing analysis

scRNA-seq data were analyzed using the scrna-seq-analysis skill pipeline
(v2.17.0). Quality control removed cells with fewer than 200 detected genes,
more than 6,000 detected genes (potential doublets), or greater than 15%
mitochondrial gene expression. Data were normalized using library size
normalization with a target sum of 10,000 counts per cell, followed by
log1p transformation. Highly variable genes (n=3,000) were identified
using the Seurat v3 method. Dimensionality reduction was performed using
PCA (50 components) followed by UMAP visualization. Cell clustering was
performed using the Leiden algorithm at resolution 0.8. Differential
expression analysis was performed using the Wilcoxon rank-sum test with
Benjamini-Hochberg correction. Novel marker discovery was performed by
comparing cluster-specific markers against reference atlases (Human Cell
Atlas, CellxGene) and literature (PubMed). Target validation was performed
by integrating data from UniProt, STRING, Reactome, and Open Targets
databases.
```
