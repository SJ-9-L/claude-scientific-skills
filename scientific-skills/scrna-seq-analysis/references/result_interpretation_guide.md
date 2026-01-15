# Result Interpretation & Publication Guide

## Overview

scRNA-seq 분석 결과를 올바르게 해석하고, Nature/Cell/Science급 저널 publication에 활용하는 가이드

---

## 1. Quality Control 결과 해석

### 1.1 QC Metrics 해석

```
QC Parameter          | Good      | Warning   | Action Required
---------------------|-----------|-----------|------------------
Cells retained       | >80%      | 60-80%    | <60% (re-evaluate)
Genes/cell (median)  | 2000-5000 | 1000-2000 | <1000 (low quality)
UMIs/cell (median)   | 5000-20000| 2000-5000 | <2000 (low depth)
MT% (median)         | <5%       | 5-10%     | >10% (dead cells)
Doublet rate         | <5%       | 5-10%     | >10% (re-filter)
```

### 1.2 Common QC Issues & Solutions

| 문제 | 원인 | 해결방법 |
|------|------|----------|
| 높은 MT% | 죽은/dying 세포 | Stricter MT threshold, FACS sorting |
| 낮은 gene count | Shallow sequencing | Re-sequence or adjust threshold |
| High doublet rate | Overloading | DoubletFinder/Scrublet 적용 |
| Batch effects | Technical variation | Batch correction (Harmony/scVI) |

### 1.3 QC 결과 보고 (논문용)

```markdown
**Quality Control**
After quality control filtering, we retained X cells (Y% of total)
across N samples. Cells were filtered based on the following criteria:
- Minimum genes detected: 200
- Maximum genes detected: 6,000 (doublet exclusion)
- Maximum mitochondrial percentage: 15%
The median genes per cell was X (range: Y-Z), and median UMIs per cell
was X (range: Y-Z), indicating high-quality single-cell capture.
```

---

## 2. Clustering 결과 해석

### 2.1 Resolution 선택 기준

| Resolution | Cluster 수 | 적합한 상황 |
|------------|-----------|-------------|
| 0.2-0.4 | 소수 (5-8) | Major cell types only |
| 0.6-0.8 | 중간 (10-15) | **일반적 권장값** |
| 1.0-1.2 | 다수 (15-25) | Subtle subtypes 필요시 |
| >1.5 | 과다 (>30) | Over-clustering 위험 |

### 2.2 Cluster 품질 평가

```python
# Good cluster indicators
good_cluster = {
    'silhouette_score': > 0.3,        # Well-separated
    'marker_genes': > 10,              # Distinct markers
    'cells_per_cluster': > 50,         # Sufficient cells
    'reproducible': True,              # Consistent across resolutions
}

# Warning signs
bad_cluster = {
    'silhouette_score': < 0.1,         # Poorly separated
    'marker_genes': < 5,               # Few specific markers
    'cells_per_cluster': < 20,         # Too few cells
    'high_mt_pct': True,               # Dying cell cluster
    'doublet_score': high,             # Doublet cluster
}
```

### 2.3 Cluster Annotation 신뢰도

| Confidence Level | Criteria | Action |
|-----------------|----------|--------|
| **High** | 3+ canonical markers, literature support | Report with confidence |
| **Medium** | 1-2 markers, partial literature | Report with caveat |
| **Low** | Novel markers only | Label as "Unknown" or validate experimentally |

---

## 3. Differential Expression 결과 해석

### 3.1 통계적 기준

```
Significance Thresholds (권장):
- Adjusted p-value: < 0.05 (FDR correction)
- Log2 Fold Change: |log2FC| > 0.5 (minimal), > 1.0 (strong)
- Detection rate: >10% in target cluster

Prioritization Formula:
Priority Score = |log2FC| × (-log10(p_adj)) × specificity
```

### 3.2 DE 결과 해석 매트릭스

| log2FC | Adj P-value | Detection | 해석 |
|--------|-------------|-----------|------|
| >2 | <0.001 | >50% | **Strong marker** - 보고 |
| 1-2 | <0.01 | 30-50% | Good marker - 보고 |
| 0.5-1 | <0.05 | 20-30% | Moderate - context 고려 |
| <0.5 | >0.05 | <20% | Weak - 보고하지 않음 |

### 3.3 Common Pitfalls

❌ **Avoid:**
- Reporting thousands of DE genes (top 50-100 sufficient)
- Ignoring detection rates (high FC but low detection = noisy)
- Over-interpreting small fold changes
- Not correcting for multiple testing

✅ **Do:**
- Focus on top markers with biological relevance
- Validate computationally (multiple methods)
- Consider effect size AND statistical significance
- Report detection rates alongside fold changes

---

## 4. Novel Gene Discovery 해석

### 4.1 Novelty 기준

```
Novelty Levels:

Level 1: TRULY NOVEL (Highest Priority)
- Gene not associated with cell type in ANY database
- No PubMed hits for "gene + cell type"
- Strong expression specificity (>0.8)
→ Highest publication value

Level 2: CONTEXT NOVEL (High Priority)
- Gene known in other contexts
- Not associated with THIS specific cell type/condition
- PubMed hits exist but not relevant
→ Good publication value

Level 3: UNDER-APPRECIATED (Medium Priority)
- Gene mentioned but not well-characterized
- Few functional studies
- Mechanism unknown
→ Moderate publication value

Level 4: CONFIRMATORY (Lower Priority)
- Gene already known marker
- Validates existing knowledge
→ Supporting data only
```

### 4.2 Novel Gene Validation Checklist

```markdown
Before claiming novel discovery, verify:

□ Expression specificity (not just high expression)
□ Consistent across biological replicates
□ Literature search confirms novelty
□ Database search confirms novelty (HCA, CellxGene)
□ Protein product exists and is detectable
□ Biologically plausible (pathway context)
□ Reproducible in independent dataset (if available)
```

### 4.3 Novel Finding 보고 방법

```markdown
**Strong claim (validated):**
"We identify GENE_X as a novel marker of [cell type], not previously
associated with this population in published literature or reference
atlases (Methods). GENE_X showed high specificity (score=0.92) with
expression in 78% of [cell type] versus 12% of other cells (log2FC=3.2,
adj.p<0.001)."

**Moderate claim (computational only):**
"Our analysis suggests GENE_X as a candidate novel marker of [cell type],
pending experimental validation. Computational analysis showed high
specificity (score=0.85, log2FC=2.5, adj.p<0.001) and absence from
existing reference atlases."
```

---

## 5. Pathway Analysis 결과 해석

### 5.1 Enrichment 결과 해석

| FDR | Fold Enrichment | Gene Count | 해석 |
|-----|-----------------|------------|------|
| <0.001 | >3x | >10 | **Highly significant** |
| <0.01 | 2-3x | 5-10 | Significant |
| <0.05 | 1.5-2x | 3-5 | Moderate (context dependent) |
| >0.05 | <1.5x | <3 | Not significant |

### 5.2 Pathway Redundancy 처리

```python
# Many pathways are redundant (e.g., "MAPK", "ERK", "RAS" all overlap)
# Solutions:

1. Use pathway clustering (e.g., enrichR, clusterProfiler)
2. Report parent terms, not all children
3. Focus on non-redundant top pathways
4. Use GSEA for better specificity

# Bad: Report 50 redundant GO terms
# Good: Report 5-10 independent biological themes
```

### 5.3 Pathway 결과 보고

```markdown
**Publication format:**
Gene set enrichment analysis revealed significant enrichment of
[pathway 1] (FDR=X, n=Y genes), [pathway 2] (FDR=X, n=Y genes),
and [pathway 3] (FDR=X, n=Y genes) in [condition/cluster].
Complete pathway analysis results are provided in Supplementary Table X.
```

---

## 6. Figure 품질 기준

### 6.1 UMAP/Visualization 체크리스트

```
Essential elements:
☑ Clear cluster separation (not over-mixed)
☑ Appropriate point size (not too large/small)
☑ Color-blind friendly palette
☑ Legend with all labels
☑ Scale bar or axis labels
☑ Title or panel label

Quality checks:
☑ No obvious batch effects visible
☑ Clusters are reproducible
☑ Annotation makes biological sense
☑ Resolution appropriate for message
```

### 6.2 Publication Figure 사양

| Journal Level | Resolution | Format | Size |
|--------------|------------|--------|------|
| Nature/Cell/Science | 300 DPI min | Vector (PDF/EPS) | 1-2 column |
| Good journals | 300 DPI | PDF/PNG | Varies |
| Preprints | 150 DPI | PNG | Flexible |

### 6.3 Color Schemes

```python
# Nature-style color palettes
nature_colors = {
    'categorical': ['#E64B35', '#4DBBD5', '#00A087', '#3C5488',
                   '#F39B7F', '#8491B4', '#91D1C2', '#DC0000'],
    'sequential': 'viridis',  # or 'magma', 'plasma'
    'diverging': 'RdBu_r',    # Red-Blue for differential
}

# Avoid: Rainbow, jet, bright neon colors
# Use: Colorblind-safe, publication-standard palettes
```

---

## 7. Statistical Reporting Standards

### 7.1 필수 보고 항목

```markdown
For each analysis, report:
- Sample sizes (cells per condition/cluster)
- Statistical test used
- Multiple testing correction method
- Significance threshold
- Effect size (fold change)
- Confidence intervals (when applicable)
```

### 7.2 P-value 보고 방식

| P-value | 보고 방식 | 그래프 표시 |
|---------|----------|-------------|
| <0.0001 | p < 0.0001 | **** |
| <0.001 | p < 0.001 | *** |
| <0.01 | p = 0.00X | ** |
| <0.05 | p = 0.0XX | * |
| ≥0.05 | p = 0.XX (NS) | ns |

### 7.3 Sample Size Justification

```markdown
**For discovery studies:**
"We analyzed X cells from Y patients/samples, providing sufficient
statistical power to detect cell populations comprising >Z% of
total cells and differential expression with effect sizes >W
(log2FC) at FDR<0.05."

**For validation studies:**
"Power analysis indicated N samples required to detect [effect]
with 80% power at α=0.05. Our study included M samples, exceeding
this requirement."
```

---

## 8. Supplementary Materials 구성

### 8.1 필수 Supplementary Tables

| Table | Content | Format |
|-------|---------|--------|
| Table S1 | Sample metadata | Excel |
| Table S2 | QC statistics per sample | Excel |
| Table S3 | Cluster marker genes (all) | Excel |
| Table S4 | Pathway enrichment results | Excel |
| Table S5 | Differential expression results | Excel |

### 8.2 필수 Supplementary Figures

| Figure | Content |
|--------|---------|
| Fig S1 | QC plots (violin, scatter) |
| Fig S2 | Clustering at multiple resolutions |
| Fig S3 | Marker gene expression heatmap |
| Fig S4 | Batch correction evaluation |
| Fig S5 | Additional validation plots |

### 8.3 Data Availability Statement

```markdown
**Standard format:**
Raw and processed single-cell RNA sequencing data have been deposited
in the Gene Expression Omnibus (GEO) under accession number GSE[XXXXX].
Processed data objects (h5ad format) and analysis code are available
at [GitHub repository URL]. Interactive data exploration is available
at [cellxgene/UCSC Cell Browser URL].

**Checklist before submission:**
☐ Raw data deposited (GEO/SRA/ENA)
☐ Processed data available (h5ad/Seurat object)
☐ Analysis code deposited (GitHub/Zenodo)
☐ Interactive viewer set up (optional but recommended)
```

---

## 9. Reviewer Response 준비

### 9.1 Common Reviewer Questions

| Question | Preparation |
|----------|-------------|
| "Why this clustering resolution?" | Test multiple, show stability |
| "Batch effects?" | Show UMAP colored by batch |
| "How novel is this gene?" | Literature/database search evidence |
| "Reproducibility?" | Split-half analysis, independent datasets |
| "Statistical power?" | Cell numbers per comparison |

### 9.2 Pre-emptive Analyses

```python
# Prepare these analyses before submission:

1. Robustness analysis
   - Clustering stability across resolutions
   - DE consistency across methods (Wilcoxon, t-test, MAST)
   - Results with/without batch correction

2. Validation analyses
   - Cross-reference with public datasets
   - Orthogonal validation (protein, FISH, flow)
   - Split-half reproducibility

3. Negative controls
   - Housekeeping gene expression (GAPDH, ACTB)
   - Known negative markers for cell types
   - Random gene set enrichment

4. Sensitivity analyses
   - Different QC thresholds
   - Different normalization methods
   - Subset analyses (by patient, batch)
```

---

## 10. Claim 강도별 표현 가이드

### 10.1 강한 Claim (실험 검증 있음)

```markdown
"We demonstrate that..."
"Our data establish..."
"We show that..."
"These results prove..."
```

### 10.2 중간 Claim (강한 computational 증거)

```markdown
"Our analysis reveals..."
"These data indicate..."
"We identify..."
"Our findings suggest..."
```

### 10.3 약한 Claim (추가 검증 필요)

```markdown
"Our data are consistent with..."
"These results may indicate..."
"This suggests the possibility that..."
"Further studies are needed to confirm..."
```

### 10.4 회피해야 할 표현

```markdown
❌ Avoid:
- "We prove..." (without experimental validation)
- "This definitively shows..."
- "Undoubtedly..."
- "For the first time ever..." (unless verified)

✅ Use:
- Hedged language appropriate to evidence level
- Specific quantitative statements
- Clear distinction between observation and interpretation
```

---

## 11. Nature/Cell/Science 수준 체크리스트

### 11.1 논문 투고 전 최종 점검

```markdown
**Data Quality:**
☐ Sufficient cell numbers (>5,000 for discovery)
☐ Appropriate biological replicates (n≥3)
☐ High-quality QC metrics
☐ Proper batch correction (if multiple batches)

**Analysis Rigor:**
☐ Multiple clustering resolutions tested
☐ Markers validated against literature
☐ Pathway analysis non-redundant
☐ Statistics properly corrected

**Novelty:**
☐ Novel finding clearly stated
☐ Comparison to existing knowledge
☐ Biological/clinical significance explained
☐ Mechanistic insight (not just observation)

**Reproducibility:**
☐ Methods detailed and reproducible
☐ Code available
☐ Data deposited
☐ Key findings validated (computationally or experimentally)

**Figures:**
☐ Publication quality (300 DPI, vector)
☐ Clear and informative
☐ Proper statistical annotations
☐ Colorblind accessible

**Writing:**
☐ Claims matched to evidence
☐ Limitations acknowledged
☐ Proper citations
☐ Clear take-home message
```

### 11.2 Journal-Specific Requirements

| Journal | Special Requirements |
|---------|---------------------|
| Nature | Reporting summary, data availability |
| Cell | STAR Methods, Key Resources Table |
| Science | Structured abstract, supplementary refs |
| Nature Methods | Benchmarking, method comparison |
| Nature Communications | Reproducibility checklist |

---

## 12. 결과 해석 Red Flags

### 12.1 의심해야 할 결과

```markdown
⚠️ Red Flags:

1. 모든 cluster가 완벽하게 분리됨
   - 가능성: Over-clustering 또는 batch effect
   - 확인: 다른 resolution, batch correction

2. DE genes가 수천 개
   - 가능성: Threshold 너무 관대
   - 확인: Stricter FDR, higher fold change cutoff

3. 알려진 marker가 발현 안 됨
   - 가능성: Data quality, wrong annotation
   - 확인: QC review, marker literature check

4. Novel gene이 너무 많음
   - 가능성: Database search 불완전
   - 확인: Multiple database search, manual literature review

5. Pathway enrichment이 general terms만
   - 가능성: Gene list too large, redundant pathways
   - 확인: Stricter gene list, pathway clustering
```

### 12.2 Troubleshooting Guide

| Problem | Likely Cause | Solution |
|---------|--------------|----------|
| Poor cluster separation | Technical noise | Better QC, batch correction |
| No DE genes | Low power | More cells, relaxed threshold |
| Too many DE genes | Threshold too loose | Stricter FDR, higher FC |
| Unexpected cell types | Contamination | Re-annotate, check markers |
| Irreproducible clusters | Parameter sensitivity | Multiple resolution, consensus |
