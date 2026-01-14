# scRNA-seq 데이터 탐색 및 확보

## 개요
이 명령은 특정 연구 주제에 대한 공개 scRNA-seq 데이터를 자동으로 찾고, 다운로드 가능한 형태로 정리합니다.

## 입력 정보 확인
사용자에게 다음을 확인하세요:
- **연구 키워드**: (예: "mouse retina endothelial", "VEGF angiogenesis")
- **종 (Species)**: Human (9606) / Mouse (10090) / 둘 다
- **관심 세포 유형**: (예: endothelial cell, tip cell, pericyte)
- **실험 조건**: (예: normal, diabetic retinopathy, hypoxia)

---

## Step 1: CellxGene Census 검색 (가장 빠름)

CellxGene Census에서 바로 사용 가능한 데이터를 먼저 확인합니다.

```python
import cellxgene_census

with cellxgene_census.open_soma() as census:
    # 1. 어떤 조직/세포 데이터가 있는지 확인
    cell_metadata = cellxgene_census.get_obs(
        census,
        "[species]",  # "homo_sapiens" 또는 "mus_musculus"
        value_filter="tissue_general == '[tissue]' and is_primary_data == True",
        column_names=["cell_type", "tissue", "disease", "assay", "dataset_id"]
    )

    # 2. 세포 유형별 분포 확인
    print(cell_metadata["cell_type"].value_counts())

    # 3. 관련 데이터셋 정보 추출
    datasets = census["census_info"]["datasets"].read().concat().to_pandas()
    relevant_datasets = datasets[datasets["dataset_id"].isin(
        cell_metadata["dataset_id"].unique()
    )]
```

**출력**:
- 사용 가능한 세포 수
- 세포 유형별 분포
- 데이터셋 목록 (논문 정보 포함)

---

## Step 2: GEO Database 검색

`geo-database` 스킬을 사용하여 GEO에서 관련 데이터셋을 검색합니다.

**검색 쿼리 구성**:
```
("[keyword]"[Title/Abstract]) AND ("single cell RNA"[Title/Abstract] OR "scRNA-seq"[Title/Abstract] OR "10x Genomics"[Title/Abstract])
```

**각 데이터셋에서 확인할 정보**:
- GSE accession number
- 샘플 수 및 조건
- 플랫폼 (10X Genomics, Smart-seq2, Drop-seq 등)
- Supplementary files 형식 (h5ad, mtx, raw fastq)
- 관련 논문 (PMID)

**출력 테이블**:
| GSE | Title | Samples | Platform | Processed? | PMID | Download |
|-----|-------|---------|----------|------------|------|----------|

---

## Step 3: PubMed에서 논문 → 데이터 추적

`pubmed-database` 스킬로 최근 관련 논문을 검색하고, Data Availability 섹션에서 데이터 accession을 추출합니다.

**검색 전략**:
```
("[keyword]") AND ("single-cell" OR "scRNA-seq") AND ("2022"[Date - Publication] : "2024"[Date - Publication])
```

**각 논문에서 추출**:
1. 제목, 저자, 저널
2. Data Availability 섹션 텍스트
3. GEO/ArrayExpress/Zenodo accession 번호
4. 코드 repository (GitHub)

---

## Step 4: 데이터 품질 평가 및 우선순위화

찾은 데이터셋들을 다음 기준으로 평가:

| 기준 | 가중치 | 설명 |
|------|--------|------|
| 관심 세포 포함 여부 | 높음 | 해당 세포 유형이 있는지 |
| 처리된 데이터 여부 | 높음 | h5ad/anndata vs raw fastq |
| 세포 수 | 중간 | 충분한 통계적 power |
| 조건 매칭 | 중간 | 실험 조건이 맞는지 |
| 플랫폼 호환성 | 낮음 | 기존 데이터와 통합 가능성 |

**우선순위 점수 계산 후 상위 5개 추천**

---

## Step 5: 데이터 다운로드 코드 생성

선택된 데이터셋에 대해 다운로드 코드를 자동 생성합니다.

### CellxGene Census 데이터:
```python
import cellxgene_census
import scanpy as sc

with cellxgene_census.open_soma() as census:
    adata = cellxgene_census.get_anndata(
        census=census,
        organism="[species]",
        obs_value_filter="[filter_condition]",
        var_value_filter="[gene_filter_if_needed]",
    )

# 저장
adata.write("cellxgene_[keyword].h5ad")
```

### GEO 데이터:
```bash
# Raw count matrix 다운로드
wget -O GSE[XXXXX]_counts.h5ad "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE[XXXXX]&format=file&file=..."

# 또는 GEOquery 사용
```

```python
import GEOparse
import scanpy as sc

# GEO에서 다운로드
gse = GEOparse.get_GEO(geo="GSE[XXXXX]", destdir="./data")

# Supplementary file 처리
adata = sc.read_h5ad("./data/GSE[XXXXX]_*.h5ad")
```

---

## Step 6: 통합 리포트 생성

### 최종 출력물:

**1. 데이터셋 요약 테이블** (`scrna_datasets_summary.csv`)
| Source | ID | Cells | Cell Types | Conditions | Format | Download Status |
|--------|-----|-------|------------|------------|--------|-----------------|

**2. 다운로드 스크립트** (`download_datasets.sh` / `download_datasets.py`)

**3. 메타데이터 정리** (`dataset_metadata.json`)
```json
{
  "datasets": [
    {
      "id": "GSE123456",
      "source": "GEO",
      "cells": 50000,
      "cell_types": ["endothelial", "pericyte", "..."],
      "conditions": ["normal", "diabetic"],
      "paper": "PMID:12345678",
      "download_url": "...",
      "local_path": "./data/GSE123456.h5ad"
    }
  ]
}
```

**4. 추천 분석 전략**
- 어떤 데이터셋을 reference로 쓸지
- 어떤 데이터셋끼리 통합 가능한지
- Batch effect 예상 정도

---

## 다음 단계 안내

데이터 확보 후 `/target-discovery` 명령으로 타겟 발굴 파이프라인을 실행하세요.

```
/target-discovery --data ./data/my_data.h5ad --reference ./data/cellxgene_reference.h5ad
```
