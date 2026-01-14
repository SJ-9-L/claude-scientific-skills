# 타겟 발굴 파이프라인 (Target Discovery Pipeline)

## 개요
scRNA-seq 데이터에서 새로운 타겟을 발굴하고, 28개 데이터베이스 + 12개 임상/의료 스킬을 통합 조회하여 Nature Reviews급 pathway figure와 함께 연구에 바로 사용할 수 있는 리포트를 생성합니다.

---

## 입력 정보 확인

사용자에게 다음을 확인하세요:

**필수:**
- **데이터 경로**: h5ad 파일 경로 (예: `./data/organoid.h5ad`)
- **비교 조건**: 어떤 조건끼리 비교할지 (예: Treatment vs Control)
- **조건 컬럼명**: adata.obs에서 조건을 나타내는 컬럼
- **관심 세포 유형**: 특정 세포만 분석할지 (예: endothelial cell)

**선택:**
- **Reference 데이터**: CellxGene 등에서 가져온 비교 데이터
- **관심 Pathway**: 특히 보고 싶은 pathway (예: VEGF signaling)
- **DEG 기준**: FDR threshold, log2FC threshold

---

## PHASE 1: scRNA-seq 데이터 분석

### Step 1.1: 데이터 로드 및 QC

```python
import scanpy as sc
import pandas as pd
import numpy as np

# 데이터 로드
adata = sc.read_h5ad("[data_path]")

# QC 메트릭 계산
adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# QC 필터링
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

print(f"QC 후 세포 수: {adata.n_obs}")
print(f"유전자 수: {adata.n_vars}")
```

### Step 1.2: 정규화 및 고변이 유전자 선택

```python
# 정규화
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Raw counts 보존 (DEG 분석용)
adata.raw = adata

# 고변이 유전자 선택
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
adata_hvg = adata[:, adata.var.highly_variable]

# 스케일링
sc.pp.scale(adata_hvg, max_value=10)
```

### Step 1.3: 차원 축소 및 클러스터링

```python
# PCA
sc.tl.pca(adata_hvg, svd_solver='arpack', n_comps=50)

# Neighbors & UMAP
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_hvg)

# Leiden 클러스터링
sc.tl.leiden(adata_hvg, resolution=0.5)

# 시각화
sc.pl.umap(adata_hvg, color=['leiden', '[condition_column]'], save='_overview.png')
```

### Step 1.4: 세포 유형 어노테이션

```python
# 마커 유전자로 세포 유형 확인
markers = {
    'Endothelial': ['PECAM1', 'CDH5', 'KDR', 'FLT1'],
    'Pericyte': ['PDGFRB', 'RGS5', 'ACTA2'],
    'Tip cell': ['ESM1', 'APLN', 'ANGPT2'],
    'Stalk cell': ['HEY1', 'JAG1', 'DLL4'],
    # 사용자 연구에 맞게 추가
}

sc.pl.dotplot(adata_hvg, markers, groupby='leiden', save='_markers.png')
```

---

## PHASE 2: 차등 발현 분석 (DEG)

### Step 2.1: PyDESeq2를 이용한 DEG 분석

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 관심 세포만 선택 (선택적)
adata_subset = adata[adata.obs['cell_type'] == '[target_cell_type]']

# DESeq2 분석
counts = adata_subset.raw.X.toarray() if hasattr(adata_subset.raw.X, 'toarray') else adata_subset.raw.X
dds = DeseqDataSet(
    counts=pd.DataFrame(counts, columns=adata_subset.var_names),
    metadata=adata_subset.obs,
    design_factors="[condition_column]"
)

dds.deseq2()
stat_res = DeseqStats(dds, contrast=["[condition_column]", "[treatment]", "[control]"])
stat_res.summary()

# 결과 추출
deg_results = stat_res.results_df
deg_results['gene'] = adata_subset.var_names
```

### Step 2.2: 유의한 DEG 필터링

```python
# 기준 설정
FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

# 필터링
significant_up = deg_results[
    (deg_results['padj'] < FDR_THRESHOLD) &
    (deg_results['log2FoldChange'] > LOG2FC_THRESHOLD)
].sort_values('log2FoldChange', ascending=False)

significant_down = deg_results[
    (deg_results['padj'] < FDR_THRESHOLD) &
    (deg_results['log2FoldChange'] < -LOG2FC_THRESHOLD)
].sort_values('log2FoldChange')

print(f"Upregulated genes: {len(significant_up)}")
print(f"Downregulated genes: {len(significant_down)}")

# 상위 타겟 선정
top_targets = pd.concat([
    significant_up.head(20),
    significant_down.head(20)
])
```

### Step 2.3: Volcano Plot 생성

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 8))
ax.scatter(deg_results['log2FoldChange'], -np.log10(deg_results['padj']),
           c='gray', alpha=0.5, s=1)
ax.scatter(significant_up['log2FoldChange'], -np.log10(significant_up['padj']),
           c='red', alpha=0.7, s=10, label='Upregulated')
ax.scatter(significant_down['log2FoldChange'], -np.log10(significant_down['padj']),
           c='blue', alpha=0.7, s=10, label='Downregulated')

# 상위 유전자 라벨링
for _, row in top_targets.head(10).iterrows():
    ax.annotate(row['gene'], (row['log2FoldChange'], -np.log10(row['padj'])))

plt.savefig('volcano_plot.png', dpi=300, bbox_inches='tight')
```

---

## PHASE 3: 28개 데이터베이스 통합 조회

상위 타겟 유전자들에 대해 다음 데이터베이스를 **순차적으로** 조회합니다.

### Step 3.1: 단백질/유전자 기본 정보

| 스킬 | 조회 내용 |
|------|----------|
| `uniprot-database` | 단백질 기능, 도메인, subcellular location, tissue specificity |
| `gene-database` | 유전자 설명, aliases, gene ontology |
| `ensembl-database` | Transcript variants, regulatory regions |
| `pdb-database` | 3D 구조 존재 여부, druggable pocket |

```
각 타겟 유전자에 대해:
1. UniProt에서 단백질 ID, 기능, 도메인 정보 조회
2. 분비 단백질인지, 막 단백질인지, 세포내 단백질인지 분류
3. 알려진 기능과 새로 발견한 발현 패턴의 연관성 분석
```

### Step 3.2: 상호작용 및 네트워크

| 스킬 | 조회 내용 |
|------|----------|
| `string-database` | 단백질-단백질 상호작용 네트워크 (confidence > 700) |
| `kegg-database` | 속한 pathway 목록, pathway position |
| `reactome-database` | 상세 pathway, upstream/downstream |

```
각 타겟에 대해:
1. STRING에서 상위 20개 상호작용 파트너 추출
2. 상호작용 파트너들이 DEG 목록에 있는지 cross-check
3. "새로운 연결" 발견 시 하이라이트
4. 네트워크 이미지 생성 (evidence 색상)
```

### Step 3.3: 임상적 의미

| 스킬 | 조회 내용 |
|------|----------|
| `clinvar-database` | 알려진 병원성 변이 |
| `cosmic-database` | 암에서의 mutation 빈도 |
| `gwas-database` | 관련 GWAS 연구 |
| `opentargets-database` | 질병 연관성, druggability score |

```
각 타겟에 대해:
1. 망막 질환 관련 변이 있는지 확인
2. Druggability score 확인 (타겟 가능성)
3. 기존에 약물 타겟으로 연구된 적 있는지
```

### Step 3.4: 약물/시약 정보

| 스킬 | 조회 내용 |
|------|----------|
| `drugbank-database` | 기존 약물, mechanism of action |
| `chembl-database` | bioactivity data, IC50 정보 |
| `pubchem-database` | 화합물 정보, 구매 가능 여부 |

```
각 타겟에 대해:
1. 이미 개발된 inhibitor/activator 존재 여부
2. 상용 항체 존재 여부
3. 연구용 시약 가용성
```

### Step 3.5: 임상 시험 정보

| 스킬 | 조회 내용 |
|------|----------|
| `clinicaltrials-database` | 진행 중인 임상시험 |
| `fda-database` | 승인된 약물 |

---

## PHASE 4: 경쟁 그룹 & 최신 동향 분석

### Step 4.1: 문헌 검색

```
각 타겟 유전자에 대해 다음 검색 수행:

PubMed 검색:
- "[Gene] AND [tissue/cell type]"
- "[Gene] AND [관심 pathway]"
- "[Gene] AND single-cell"

bioRxiv 검색:
- 최근 6개월 preprint

OpenAlex 검색:
- 저자별 그룹핑
- Citation 네트워크 분석
```

### Step 4.2: 경쟁 분석 매트릭스

| 연구 그룹 | 주요 논문 | 연구 방향 | 우리 연구와 차별점 |
|----------|----------|----------|------------------|
| Lab A | PMID:xxx | VEGF-Notch | 우리는 새로운 X 발견 |
| Lab B | PMID:xxx | Tip cell | 다른 모델 시스템 |

### Step 4.3: Unmet Needs 분석

```
문헌에서 발견한 정보를 바탕으로:
1. 아직 밝혀지지 않은 질문들 정리
2. 우리 데이터가 답할 수 있는 것들 매칭
3. "Novel finding" 포인트 도출
```

---

## PHASE 5: Nature Reviews급 Pathway Figure 생성

### Step 5.1: Figure 구성 요소 정의

`scientific-schematics` 스킬을 사용하여 다음 요소를 포함한 figure 생성:

```
Figure 구성:
┌─────────────────────────────────────────────────────────────┐
│                                                             │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────┐     │
│  │  Ligands    │───▶│  Receptors  │───▶│  Signaling  │     │
│  │  (VEGF-A,   │    │  (VEGFR2,   │    │  cascades   │     │
│  │   Ang-1/2)  │    │   Tie2)     │    │             │     │
│  └─────────────┘    └─────────────┘    └─────────────┘     │
│        │                  │                  │              │
│        ▼                  ▼                  ▼              │
│  ┌─────────────────────────────────────────────────┐       │
│  │              Cellular Responses                 │       │
│  │  • Tip cell migration    • Vascular formation  │       │
│  │  • Proliferation         • Permeability        │       │
│  └─────────────────────────────────────────────────┘       │
│                                                             │
│  ★ = Our novel findings (highlighted in red)               │
│  ? = Unmet needs / Future directions                       │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

### Step 5.2: Figure 생성 프롬프트

```
scientific-schematics 스킬에 전달할 정보:
1. Pathway 구성요소 목록 (DB에서 추출)
2. 각 노드 간 관계 (activation/inhibition)
3. 우리가 발견한 새로운 연결 (하이라이트)
4. 색상 코드 정의
5. 레이아웃 선호도 (horizontal/vertical)
```

---

## PHASE 6: 최종 리포트 생성

### Output 파일들:

**1. Executive Summary** (`target_discovery_summary.md`)
```markdown
# Target Discovery Report

## Key Findings
- 발견된 상위 타겟 N개
- Novel pathway connections
- Druggable targets

## Top 5 Targets
| Rank | Gene | log2FC | Function | Druggability | Novel? |
|------|------|--------|----------|--------------|--------|

## Recommended Next Steps
1. 실험 검증 우선순위
2. 추가 분석 제안
```

**2. 상세 타겟 정보** (`target_details/[GENE].md`)
```markdown
# [GENE NAME] - Detailed Analysis

## Basic Information
- UniProt, Gene description

## Expression Pattern
- scRNA 분석 결과

## Network Context
- STRING interaction partners
- Pathway involvement

## Clinical Relevance
- Disease associations
- Existing drugs

## Literature Summary
- Key papers
- Competing groups

## Experimental Validation Plan
- Suggested experiments
- Available reagents
```

**3. Pathway Figure** (`figures/pathway_figure.png`, `.svg`)

**4. 데이터 파일들**
- `deg_results.csv` - 전체 DEG 결과
- `target_database_query.json` - DB 조회 결과
- `literature_analysis.json` - 문헌 분석 결과

**5. 참고문헌** (`references.bib`)

---

## 다음 단계

타겟 선정 후 `/experiment-design` 명령으로 구체적인 실험 설계를 진행하세요.

```
/experiment-design --target [GENE_NAME] --model "mouse retina" --phenotype "tip cell migration"
```
