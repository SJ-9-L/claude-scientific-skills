# 🔬 Vascular Biology Research Tools

VEGF/VEGFR, retina vasculature, organoid 연구를 위한 Claude Scientific Skills 기반 연구 도구 모음입니다.

---

## 📋 사용 가능한 명령어

| 명령어 | 용도 | 입력 | 출력 |
|--------|------|------|------|
| `/find-scrna-data` | scRNA-seq 공개 데이터 탐색 및 확보 | 키워드, 종, 세포 유형 | 데이터셋 목록, 다운로드 코드 |
| `/target-discovery` | 타겟 발굴 파이프라인 | scRNA 데이터 (h5ad) | DEG, DB 조회 결과, Pathway Figure |
| `/experiment-design` | Deep dive 실험 설계 | 타겟 유전자, 모델, 표현형 | 실험 매트릭스, 시약 목록, Crosstalk 다이어그램 |

---

## 🚀 Quick Start

### 1. 데이터 확보 (새 프로젝트 시작 시)
```
/find-scrna-data

키워드: "mouse retina endothelial VEGF"
종: Mouse
세포: endothelial cell, tip cell
```

### 2. 타겟 발굴 (데이터 있을 때)
```
/target-discovery

데이터: ./data/my_organoid.h5ad
비교: Treatment vs Control
세포: endothelial cell
```

### 3. 실험 설계 (타겟 선정 후)
```
/experiment-design

타겟: VEGFR2
모델: mouse retina, hiPSC organoid
표현형: tip cell migration
비교: VEGFR2 KO vs VEGFR3 KO vs WT
```

---

## 📊 Workflow Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                     Complete Research Pipeline                   │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────────┐                                           │
│  │ /find-scrna-data │  ← 공개 데이터 탐색                        │
│  │                  │    CellxGene, GEO, PubMed                 │
│  └────────┬─────────┘                                           │
│           │                                                      │
│           ▼                                                      │
│  ┌──────────────────┐                                           │
│  │ /target-discovery│  ← 타겟 발굴                               │
│  │                  │    DEG → DB 28개 조회 → Figure             │
│  └────────┬─────────┘                                           │
│           │                                                      │
│           ▼                                                      │
│  ┌──────────────────┐                                           │
│  │/experiment-design│  ← 실험 설계                               │
│  │                  │    시약, KO 비교, Crosstalk                │
│  └──────────────────┘                                           │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## 🗄️ 사용되는 Scientific Skills

### 데이터베이스 (28개)
| Category | Skills |
|----------|--------|
| Protein/Gene | `uniprot-database`, `gene-database`, `ensembl-database`, `pdb-database` |
| Pathway | `kegg-database`, `reactome-database`, `string-database` |
| Clinical | `clinvar-database`, `cosmic-database`, `gwas-database`, `clinicaltrials-database` |
| Drug/Target | `drugbank-database`, `chembl-database`, `opentargets-database`, `pubchem-database` |
| Literature | `pubmed-database`, `openalex-database`, `biorxiv-database` |
| scRNA | `cellxgene-census`, `geo-database` |

### 분석 도구
| Category | Skills |
|----------|--------|
| scRNA Analysis | `scanpy`, `scvi-tools`, `anndata` |
| DEG Analysis | `pydeseq2` |
| Visualization | `matplotlib`, `plotly`, `scientific-schematics` |
| Writing | `scientific-writing`, `literature-review` |

---

## 📁 Output 파일 구조

### /target-discovery 출력
```
target_discovery_output/
├── deg_results.csv              # DEG 전체 결과
├── target_summary.md            # 상위 타겟 요약
├── database_queries/
│   ├── uniprot_results.json
│   ├── string_networks/
│   └── pathway_enrichment.csv
├── competitor_analysis.md       # 경쟁 그룹 분석
├── figures/
│   ├── volcano_plot.png
│   ├── pathway_figure.png
│   └── pathway_figure.svg
└── references.bib
```

### /experiment-design 출력
```
experiment_design_output/
├── executive_summary.md         # 핵심 요약
├── reagents/
│   ├── genetic_tools.csv        # siRNA, CRISPR, KO mice
│   ├── small_molecules.csv      # 억제제
│   ├── antibodies.csv           # 항체
│   └── order_list.xlsx          # 주문 목록
├── experiment_matrix.csv        # 실험 설계
├── ko_comparison.csv            # KO 표현형 비교
├── figures/
│   └── crosstalk_diagram.png
└── references.bib
```

---

## ⚠️ 주의사항

1. **데이터 크기**: CellxGene에서 대용량 데이터 다운로드 시 메모리 확인
2. **API 제한**: 일부 DB는 rate limit 있음 (1초 간격 권장)
3. **버전 관리**: CellxGene Census 버전 명시 권장 (재현성)
4. **시약 확인**: 문헌 기반 정보이므로 최신 vendor 카탈로그 확인 필요

---

## 🔧 Customization

각 명령어의 `.md` 파일을 수정하여 연구 분야에 맞게 커스터마이즈할 수 있습니다:

- `find-scrna-data.md`: 검색 키워드, 데이터 소스 추가
- `target-discovery.md`: DEG 기준, 조회할 DB 변경
- `experiment-design.md`: 관심 pathway, 시약 카테고리 추가

---

## 📞 지원

문제 발생 시 Claude에게 직접 물어보세요:
- "이 에러 어떻게 해결해?"
- "DEG 기준을 바꾸고 싶어"
- "다른 pathway도 추가해줘"
