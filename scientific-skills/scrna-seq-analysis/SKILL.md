---
name: scrna-seq-analysis
description: Comprehensive scRNA-seq analysis pipeline for novel gene/protein discovery, pathway analysis, competitive landscape research, and Nature Reviews-quality figure generation. Integrates 40+ scientific databases for deep pathway validation with experimental evidence including KO models, cell-type specifics, and molecular interactions.
license: MIT
allowed-tools: [Read, Write, Edit, Bash, WebSearch, WebFetch]
metadata:
    skill-author: K-Dense Inc.
---

# scRNA-seq Analysis & Pathway Discovery

## Overview

이 스킬은 단일세포 RNA 시퀀싱(scRNA-seq) 데이터 분석부터 새로운 유전자/단백질 발견, pathway 규명, 경쟁 동향 분석, 그리고 Nature Reviews급 publication figure 생성까지 전체 연구 워크플로우를 지원합니다.

**핵심 기능:**
1. **scRNA-seq 분석**: 데이터 전처리, 클러스터링, 차별 발현 유전자 발굴
2. **Novel Discovery**: 새로운 유전자/단백질 후보 발굴 및 검증
3. **Pathway Deep Dive**: 실험적 검증 가능한 수준의 pathway 분석
4. **Competitive Intelligence**: 경쟁 그룹 및 최신 연구 동향 분석
5. **Figure Generation**: Nature Reviews 수준의 pathway schematic 생성

## When to Use This Skill

이 스킬을 사용해야 하는 경우:
- scRNA-seq 데이터에서 새로운 마커 유전자나 발현 단백질을 발굴할 때
- 특정 pathway의 crosstalk, interaction, regulation을 deep dive할 때
- KO mouse, specific cell type 등 구체적 실험 조건에서의 pathway 변화를 분석할 때
- 경쟁 연구 그룹의 최신 동향과 unmet needs를 파악할 때
- Nature Reviews급 main figure에 들어갈 pathway schematic이 필요할 때
- 실험 설계를 위한 ligands, receptors, antibodies, siRNAs 정보가 필요할 때

## 통합 데이터베이스 (40개)

### 생물정보학 데이터베이스 (28개)
| 카테고리 | 데이터베이스 | 용도 |
|---------|------------|------|
| **Protein** | UniProt, PDB, AlphaFold, BRENDA | 단백질 서열, 구조, 효소 정보 |
| **Pathway** | KEGG, Reactome, STRING | Pathway, interaction network |
| **Genomics** | Ensembl, NCBI Gene, GEO | 유전자 정보, 발현 데이터 |
| **Single-cell** | CellxGene, Human Cell Atlas | scRNA-seq reference |
| **Chemical** | ChEMBL, PubChem, DrugBank, ZINC | 화합물, 약물 정보 |
| **Literature** | PubMed, bioRxiv, OpenAlex | 논문 검색 |
| **Target** | Open Targets, PharmGKB | 약물 타겟 정보 |

### 임상/의료 데이터베이스 (12개)
| 데이터베이스 | 용도 |
|------------|------|
| ClinicalTrials.gov | 진행 중인 임상 시험 |
| ClinVar | 유전자 변이 임상 의미 |
| COSMIC | 암 체세포 변이 |
| FDA FAERS | 약물 부작용 |
| GWAS Catalog | 유전자-질병 연관성 |
| DisGeNET | 유전자-질병 네트워크 |
| OMIM | 유전 질환 정보 |
| Orphanet | 희귀질환 정보 |
| CTD | 화학물질-유전자-질병 상호작용 |
| PharmGKB | 약물유전체학 |
| GTEx | 조직별 유전자 발현 |
| TCGA | 암 유전체 데이터 |

## Core Workflow

### Phase 1: scRNA-seq Data Analysis

```python
# Step 1: Data Loading & QC
import scanpy as sc
import anndata as ad
from scripts.scrna_pipeline import ScRNAAnalyzer

analyzer = ScRNAAnalyzer()
adata = analyzer.load_data("data/sample.h5ad")
adata = analyzer.quality_control(adata,
    min_genes=200,
    max_genes=5000,
    max_mt_pct=10
)

# Step 2: Normalization & Feature Selection
adata = analyzer.preprocess(adata,
    normalize_total=True,
    log_transform=True,
    n_top_genes=3000
)

# Step 3: Dimensionality Reduction & Clustering
adata = analyzer.reduce_dimensions(adata, n_pcs=50)
adata = analyzer.cluster(adata, resolution=0.8)

# Step 4: Differential Expression Analysis
markers = analyzer.find_markers(adata,
    groupby='leiden',
    method='wilcoxon',
    n_genes=100
)

# Step 5: Novel Gene Discovery
novel_candidates = analyzer.discover_novel_genes(adata,
    reference_db=['cellxgene', 'hca'],
    novelty_threshold=0.8
)
```

### Phase 2: Pathway Deep Dive Analysis

```python
from scripts.pathway_analyzer import PathwayDeepDive

# Initialize with target pathway
pathway = PathwayDeepDive(
    target="VEGF signaling",
    species="mouse",
    context="angiogenesis"
)

# Get comprehensive pathway analysis
analysis = pathway.analyze(
    # Experimental conditions
    knockouts=["VEGFR2", "VEGFR3", "Notch1"],
    treatments=["DAPT", "VEGF-A", "anti-VEGF"],
    cell_types=["tip cell", "stalk cell", "phalanx cell"],

    # Phenotypes of interest
    phenotypes=[
        "tip cell migration",
        "venous endothelial cell migration",
        "vascular formation",
        "sprouting angiogenesis"
    ],

    # Analysis depth
    include_crosstalk=True,
    include_regulation=True,
    experimental_evidence_only=True
)

# Get interaction details
interactions = pathway.get_interactions(
    types=["ligand-receptor", "protein-protein", "signaling"],
    evidence_sources=["experimental", "database"],
    min_confidence=0.7
)

# Get molecular tools for validation
tools = pathway.get_validation_tools(
    include=["antibodies", "siRNAs", "inhibitors", "activators"],
    validated_only=True
)
```

### Phase 3: Competitive Landscape Analysis

```python
from scripts.competitive_analysis import CompetitiveLandscape

landscape = CompetitiveLandscape(
    topic="VEGF signaling in tumor angiogenesis",
    timeframe="2020-2024"
)

# Identify key research groups
groups = landscape.identify_key_groups(
    min_publications=5,
    include_industry=True
)

# Analyze research trends
trends = landscape.analyze_trends(
    include_preprints=True,
    include_patents=True
)

# Find unmet needs and gaps
gaps = landscape.find_research_gaps(
    current_knowledge=analysis,
    identify_opportunities=True
)

# Generate competitive report
report = landscape.generate_report(
    format="detailed",
    include_visualizations=True
)
```

### Phase 4: Nature Reviews Figure Generation

```python
from scripts.figure_generator import NatureReviewsFigure

# Create publication-quality pathway figure
figure = NatureReviewsFigure(
    style="nature_reviews",
    dimensions=(180, 120),  # mm
    resolution=300  # dpi
)

# Add pathway components
figure.add_pathway(
    pathway_data=analysis,
    layout="hierarchical",
    show_interactions=True,
    show_regulation=True,
    highlight_novel=novel_candidates
)

# Add experimental evidence annotations
figure.add_annotations(
    ko_effects=analysis.ko_effects,
    treatment_effects=analysis.treatment_effects,
    cell_type_specificity=True
)

# Add legend and labels
figure.add_legend(
    include_symbols=True,
    include_evidence_types=True
)

# Export in multiple formats
figure.export("vegf_pathway.svg")
figure.export("vegf_pathway.pdf")
figure.export("vegf_pathway.png", dpi=600)
```

## Pathway Deep Dive: Detailed Analysis

### Experimental Evidence Categories

이 스킬은 단순한 "upregulation/downregulation" 정보가 아닌, 실험적으로 검증 가능한 구체적 정보를 제공합니다:

#### 1. Knockout/Knockdown Studies
```python
# Example: VEGFR2 KO analysis
ko_analysis = pathway.analyze_knockout(
    gene="VEGFR2",
    species="mouse",
    context="developmental_angiogenesis"
)

# Returns:
# {
#     "phenotypes": [
#         {"name": "tip cell migration", "effect": "abolished", "evidence": "PMID:12345"},
#         {"name": "vascular density", "effect": "reduced 80%", "evidence": "PMID:12346"}
#     ],
#     "compensatory_pathways": ["VEGFR1", "VEGFR3"],
#     "downstream_targets_affected": ["ERK1/2", "AKT", "PLCγ"],
#     "cell_types_affected": ["tip cells", "stalk cells"]
# }
```

#### 2. Treatment/Perturbation Studies
```python
# Example: DAPT (gamma-secretase inhibitor) analysis
treatment = pathway.analyze_treatment(
    compound="DAPT",
    target_pathway="Notch signaling",
    crosstalk_with="VEGF signaling"
)

# Returns detailed effects on:
# - Direct targets: Notch1 ICD cleavage blocked
# - Downstream effects: HES1, HEY1 downregulation
# - Crosstalk effects: Enhanced tip cell formation
# - Cell-type specific responses
# - Dosage-dependent effects
```

#### 3. Cell-Type Specific Analysis
```python
# Analyze pathway in specific cell types
cell_analysis = pathway.analyze_by_cell_type(
    pathway="VEGF-Notch crosstalk",
    cell_types={
        "tip_cell": ["high DLL4", "high VEGFR2", "low Notch activity"],
        "stalk_cell": ["high Jagged1", "high Notch activity", "low VEGFR2"],
        "phalanx_cell": ["low proliferation", "high VE-cadherin"]
    }
)
```

### Molecular Interaction Categories

```python
# Get comprehensive molecular interactions
interactions = pathway.get_molecular_interactions(
    categories=[
        "ligand_receptor",      # VEGF-A → VEGFR2
        "receptor_coreceptor",  # VEGFR2 + NRP1
        "kinase_substrate",     # VEGFR2 → PLCγ
        "transcription",        # ERK → ETS factors
        "protein_complex",      # DLL4-Notch1
        "signaling_cascade"     # Full signaling chains
    ]
)
```

### Available Molecular Tools

```python
# Get validated reagents for experimental validation
reagents = pathway.get_reagents(
    target="VEGFR2",
    types={
        "antibodies": {
            "blocking": ["DC101 (mouse)", "Ramucirumab (human)"],
            "detection": ["AF357 (WB)", "Clone 55B11 (IHC)"]
        },
        "small_molecules": {
            "inhibitors": ["Sunitinib", "Sorafenib", "Axitinib"],
            "activators": ["VEGF-A165"]
        },
        "genetic_tools": {
            "siRNA": ["Validated sequences from Dharmacon"],
            "shRNA": ["MISSION shRNA library"],
            "CRISPR": ["Guide RNA sequences"]
        },
        "recombinant_proteins": ["rhVEGF-A", "rhVEGFR2-Fc"]
    }
)
```

## Database Integration Examples

### Multi-Database Query for Novel Target

```python
from scripts.database_integrator import DatabaseIntegrator

db = DatabaseIntegrator()

# Comprehensive target analysis
target_analysis = db.analyze_target(
    gene="DLL4",
    queries={
        # Protein information
        "uniprot": ["sequence", "domains", "modifications"],
        "pdb": ["structures", "binding_sites"],
        "alphafold": ["predicted_structure"],

        # Expression data
        "gtex": ["tissue_expression"],
        "hpa": ["protein_expression", "subcellular_location"],
        "cellxgene": ["single_cell_expression"],

        # Pathway context
        "kegg": ["pathways", "interactions"],
        "reactome": ["reactions", "regulators"],
        "string": ["protein_network", "functional_partners"],

        # Disease relevance
        "clinvar": ["variants", "clinical_significance"],
        "cosmic": ["cancer_mutations"],
        "disgenet": ["disease_associations"],

        # Drug/Chemical
        "chembl": ["bioactivities", "assays"],
        "drugbank": ["drug_interactions"],
        "opentargets": ["tractability", "safety"]
    }
)

# Clinical context
clinical = db.get_clinical_context(
    gene="DLL4",
    queries={
        "clinicaltrials": ["ongoing_trials", "completed_trials"],
        "fda": ["approved_drugs", "adverse_events"],
        "pharmgkb": ["pharmacogenomics"]
    }
)
```

### Literature & Competitive Analysis

```python
from scripts.competitive_analysis import LiteratureAnalyzer

lit = LiteratureAnalyzer()

# Comprehensive literature search
papers = lit.search(
    query="DLL4 Notch angiogenesis",
    databases=["pubmed", "biorxiv", "pmc"],
    date_range="2020-2024",
    filters={
        "article_type": ["research", "review"],
        "organism": ["human", "mouse"]
    }
)

# Identify key research groups
groups = lit.identify_research_groups(
    papers=papers,
    min_publications=3,
    include_affiliations=True
)

# Analyze citation networks
network = lit.analyze_citations(
    seed_papers=papers[:20],
    depth=2
)

# Find research gaps
gaps = lit.identify_gaps(
    current_literature=papers,
    known_pathways=["VEGF", "Notch", "Wnt"],
    suggest_experiments=True
)
```

## Nature Reviews Figure Generation

### Figure Components

```python
from scripts.figure_generator import PathwaySchematic

# Create comprehensive pathway figure
fig = PathwaySchematic(
    title="VEGF-Notch Crosstalk in Angiogenesis",
    style="nature_reviews"
)

# Add cellular compartments
fig.add_compartment("tip_cell", position="left")
fig.add_compartment("stalk_cell", position="right")
fig.add_compartment("extracellular", position="top")

# Add molecules with proper symbols
molecules = {
    "ligands": ["VEGF-A", "DLL4", "Jagged1"],
    "receptors": ["VEGFR2", "NRP1", "Notch1"],
    "kinases": ["PLCγ", "ERK1/2", "AKT"],
    "transcription_factors": ["HES1", "HEY1", "ETS"]
}

for category, mols in molecules.items():
    fig.add_molecules(mols, category=category)

# Add interactions with evidence types
fig.add_interaction("VEGF-A", "VEGFR2",
    type="activation",
    evidence="experimental",
    refs=["PMID:12345"]
)

fig.add_interaction("DLL4", "Notch1",
    type="lateral_inhibition",
    evidence="experimental",
    refs=["PMID:23456"]
)

# Add experimental annotations
fig.add_annotation(
    text="DAPT blocks",
    target="Notch1_cleavage",
    style="inhibitor"
)

fig.add_annotation(
    text="VEGFR2 KO: No tip cells",
    position="tip_cell",
    style="ko_effect"
)

# Add legend
fig.add_legend(
    show_interaction_types=True,
    show_evidence_levels=True,
    show_cell_types=True
)

# Export
fig.save("vegf_notch_crosstalk.svg", format="svg")
fig.save("vegf_notch_crosstalk.pdf", format="pdf")
```

### SVG Output Structure

생성되는 SVG는 다음 구조를 가집니다:
- **Molecules**: 표준 기호 사용 (수용체, 키나아제, 전사인자 등)
- **Interactions**: 화살표, 억제선, 복합체 형성 표시
- **Annotations**: 실험 결과, KO 효과, 약물 표적 등
- **Legend**: 모든 기호 설명
- **Scale**: Nature Reviews 가이드라인 준수

## Best Practices

### 1. Data Analysis
- 항상 QC 리포트를 먼저 확인하고 threshold 조정
- Multiple resolution으로 clustering 수행 후 생물학적 의미 확인
- Novel gene 발굴 시 multiple reference database와 비교

### 2. Pathway Analysis
- 단순 pathway enrichment가 아닌 experimental evidence 기반 분석
- Crosstalk 분석 시 cell-type specificity 고려
- KO/treatment effect 해석 시 compensatory mechanism 확인

### 3. Competitive Analysis
- Preprint 포함하여 최신 동향 파악
- Patent 검색으로 산업계 동향 확인
- Citation network으로 핵심 연구자/논문 식별

### 4. Figure Generation
- Nature Reviews style guide 준수
- 모든 interaction에 evidence reference 포함
- Cell-type specific context 명확히 표시

## Bundled Resources

### scripts/scrna_pipeline.py
scRNA-seq 분석을 위한 comprehensive pipeline:
- Data loading (10X, h5ad, CSV)
- Quality control with visualization
- Normalization and feature selection
- Dimensionality reduction (PCA, UMAP, t-SNE)
- Clustering (Leiden, Louvain)
- Differential expression analysis
- Novel gene discovery

### scripts/pathway_analyzer.py
Deep pathway analysis module:
- Multi-database integration
- KO/treatment effect analysis
- Cell-type specific pathway analysis
- Crosstalk and regulation mapping
- Molecular interaction extraction
- Validation reagent lookup

### scripts/competitive_analysis.py
Competitive landscape analysis:
- Literature search across multiple databases
- Research group identification
- Trend analysis
- Gap identification
- Report generation

### scripts/figure_generator.py
Publication-quality figure generation:
- Nature Reviews style templates
- Pathway schematic generation
- Interaction visualization
- Annotation system
- Multi-format export

### scripts/database_integrator.py
40+ database integration:
- Unified query interface
- Cross-database ID mapping
- Result aggregation
- Evidence scoring

### references/
- `scrna_analysis_pipeline.md`: Detailed scRNA-seq analysis guide
- `pathway_deep_dive.md`: Experimental evidence categories
- `database_integration.md`: All database APIs and queries
- `competitive_analysis.md`: Literature and trend analysis
- `figure_generation.md`: Figure style guide and templates
- `molecular_tools.md`: Antibodies, siRNAs, inhibitors reference

## Example Use Cases

### Use Case 1: Novel Endothelial Marker Discovery

```python
# scRNA-seq of developing vasculature
analyzer = ScRNAAnalyzer()
adata = analyzer.full_pipeline("vascular_dev.h5ad")

# Find novel tip cell markers
novel_markers = analyzer.discover_novel_genes(
    adata,
    cluster="tip_cells",
    known_markers=["DLL4", "VEGFR2", "CXCR4"],
    novelty_criteria="not_in_literature"
)

# Validate through pathway analysis
for gene in novel_markers[:5]:
    pathway.validate_candidate(
        gene=gene,
        expected_pathway="angiogenesis",
        check_databases=["uniprot", "string", "reactome"]
    )
```

### Use Case 2: Drug Target Pathway Analysis

```python
# Analyze anti-VEGF resistance pathway
pathway = PathwayDeepDive(
    target="anti-VEGF resistance",
    disease_context="colorectal_cancer"
)

# Find compensatory pathways
compensatory = pathway.find_compensatory_pathways(
    blocked_pathway="VEGF-VEGFR2",
    search_databases=["kegg", "reactome", "literature"]
)

# Get druggable targets in compensatory pathways
targets = pathway.find_druggable_targets(
    pathways=compensatory,
    criteria=["approved_drug_exists", "clinical_trial"]
)

# Generate combination therapy figure
figure = NatureReviewsFigure()
figure.create_combination_therapy_schematic(
    primary_pathway="VEGF",
    compensatory_pathways=compensatory,
    drug_targets=targets
)
```

## Integration with Other Skills

이 스킬은 다음 스킬들과 함께 사용할 때 더욱 효과적입니다:

- **scanpy**: 상세한 scRNA-seq 분석
- **biopython**: 서열 분석
- **kegg-database**: Pathway 데이터
- **reactome-database**: Reaction 데이터
- **string-database**: Protein interaction
- **pubmed-database**: 문헌 검색
- **scientific-schematics**: Figure 생성
- **generate-image**: 이미지 생성
- **literature-review**: 체계적 문헌 분석

## Summary

이 스킬은 scRNA-seq 데이터 분석부터 publication-ready figure 생성까지 전체 연구 워크플로우를 지원합니다:

1. **Data Analysis**: QC, normalization, clustering, DE analysis
2. **Novel Discovery**: Reference database 비교를 통한 신규 유전자/단백질 발굴
3. **Pathway Deep Dive**: Experimental evidence 기반 상세 pathway 분석
4. **Competitive Analysis**: 연구 동향 및 gap 분석
5. **Figure Generation**: Nature Reviews 수준의 pathway schematic

핵심은 단순한 "up/down regulation"이 아닌, KO mouse, specific cell types, molecular interactions 등 **실험적으로 검증 가능한 수준의 상세 정보**를 제공하는 것입니다.
