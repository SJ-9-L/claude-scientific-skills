# Deep Dive 실험 설계 도구 (Experiment Design Tool)

## 개요
단순한 "upregulation/downregulation" 수준이 아닌, 구체적인 실험 조건(KO mouse, 특정 세포 유형, 특정 표현형)에서의 pathway crosstalk 정보와 사용 가능한 모든 시약(항체, 억제제, siRNA 등)을 comprehensive하게 제공합니다.

---

## 입력 정보 확인

사용자에게 다음을 확인하세요:

**필수:**
- **타겟 유전자/단백질**: (예: VEGFR2, KDR, Notch1)
- **연구 모델**: (예: mouse retina, hiPSC organoid, mESC organoid)
- **관심 표현형**: (예: tip cell migration, vascular formation, sprouting angiogenesis)
- **관심 세포 유형**: (예: tip cell, stalk cell, endothelial cell, pericyte)

**선택:**
- **비교 조건**: (예: VEGFR2 KO vs WT, Hypoxia vs Normoxia)
- **기존 처리 조건**: (예: DAPT, VEGF-A treatment)
- **관심 crosstalk pathway**: (예: Notch, Wnt, TGF-β)

---

## PHASE 1: 타겟 단백질 심층 분석

### Step 1.1: 기본 정보 수집

```
[타겟]에 대해 다음 정보를 수집:

UniProt 조회:
- Protein function (detailed)
- Domain structure (kinase domain, ligand-binding domain 등)
- Post-translational modifications
- Subcellular location
- Tissue expression pattern

Gene Database 조회:
- Official symbol, aliases
- Gene ontology (BP, MF, CC)
- Orthologs (human, mouse 비교)
- Known phenotypes from knockout studies
```

### Step 1.2: 알려진 KO/Mutant 표현형

**PubMed 심층 검색**:
```
검색 쿼리 목록:
1. "[Gene] knockout mouse phenotype"
2. "[Gene] conditional knockout [tissue]"
3. "[Gene] null mutant [organism]"
4. "[Gene] heterozygous [phenotype]"
5. "[Gene] -/- embryonic lethal"
6. "[Gene] endothelial specific knockout"
7. "[Gene] inducible knockout retina"
```

**출력 테이블**:
| Gene | KO Type | Organism | Tissue | Phenotype | Viable? | PMID |
|------|---------|----------|--------|-----------|---------|------|
| Kdr | Global | Mouse | - | Embryonic lethal E8.5, no blood islands | No | 9697698 |
| Kdr | EC-specific | Mouse | Retina | Reduced vascular density, fewer tip cells | Yes | 15364902 |

---

## PHASE 2: Pathway Crosstalk 분석

### Step 2.1: 직접 상호작용 파트너

```
STRING Database (Confidence > 900, Physical only):

[타겟]과 직접 결합하는 단백질:
- 각 파트너의 기능
- 결합 부위 (어느 도메인끼리)
- 결합의 결과 (activation/inhibition)
```

### Step 2.2: Pathway 맥락 분석

**KEGG/Reactome 통합 분석**:
```
[타겟]이 속한 모든 Pathway:

1. VEGF signaling pathway (hsa04370)
   - Position: Receptor (membrane)
   - Upstream: VEGF-A, VEGF-B, PlGF
   - Downstream: PLCγ, PI3K/Akt, MAPK cascade
   - Crosstalk with: Notch, Hippo, Wnt

2. Angiogenesis (R-HSA-194138)
   - Role: Primary receptor for sprouting
   - Key interactions: DLL4-Notch, Ang-Tie system

각 pathway에 대해:
- [타겟]의 위치 (upstream/midstream/downstream)
- 주요 조절자 (regulators)
- 주요 effector
- 다른 pathway와의 crosstalk point
```

### Step 2.3: Crosstalk 상세 분석

**주요 Crosstalk Pathway별 분석**:

```
[타겟: VEGFR2] × [Pathway: Notch] Crosstalk:

1. 연결 메커니즘:
   - VEGFR2 activation → DLL4 upregulation in tip cells
   - DLL4 → Notch1 activation in adjacent cells
   - Notch1 → VEGFR2 downregulation (lateral inhibition)

2. 세포 수준 결과:
   - Tip cell: High VEGFR2, High DLL4, Low Notch
   - Stalk cell: Low VEGFR2, Low DLL4, High Notch

3. 실험적 증거:
   | Perturbation | Result | Cell Type | Reference |
   |--------------|--------|-----------|-----------|
   | DAPT (γ-secretase inh) | ↑ Tip cells | EC | PMID:xxx |
   | DLL4 blocking Ab | ↑ Sprouting | Retina | PMID:xxx |
   | Notch1 EC-KO | ↑ Tip cell genes | EC | PMID:xxx |
```

---

## PHASE 3: 구체적 실험 조건별 문헌 검색

### Step 3.1: KO Mouse 비교 실험 검색

```
검색 패턴:
"[Gene1] knockout AND [Gene2] knockout comparison"
"[Gene1] -/- AND [Gene2] -/- [tissue]"
"double knockout [Gene1] [Gene2] [phenotype]"

예시 검색:
- "VEGFR2 knockout AND VEGFR3 knockout retina"
- "Kdr -/- AND Flt4 -/- angiogenesis"
- "DLL4 knockout tip cell VEGFR2"
```

**결과 테이블**:
| Genotype 1 | Genotype 2 | Comparison | Tissue | Finding | PMID |
|------------|------------|------------|--------|---------|------|
| Kdr EC-KO | WT | Kdr KO vs WT | Retina | 50% reduction in vascular density | xxx |
| Kdr EC-KO | Flt4 EC-KO | Kdr vs Flt4 | Retina | Kdr: tip cell defect, Flt4: lymphatic defect | xxx |

### Step 3.2: 약물/억제제 처리 실험 검색

```
검색 패턴:
"[Drug/Inhibitor] treatment [Gene] [cell type]"
"[Drug] AND [phenotype] AND [model]"

예시 검색:
- "DAPT treatment tip cell migration"
- "DC101 VEGFR2 blocking retina angiogenesis"
- "Axitinib endothelial cell sprouting"
```

**결과 테이블**:
| Treatment | Target | Dose | Duration | Model | Effect on [Phenotype] | PMID |
|-----------|--------|------|----------|-------|----------------------|------|
| DAPT | γ-secretase | 10μM | 48h | HUVEC | ↑ Tip cell markers | xxx |
| DC101 | VEGFR2 | 40mg/kg | 5 days | Mouse retina | ↓ Vascular density 60% | xxx |

### Step 3.3: 세포 유형별 특이적 효과 검색

```
검색 패턴:
"[Gene] [specific cell type] specific"
"[Gene] expression tip cell vs stalk cell"
"[Gene] endothelial heterogeneity"

예시 검색:
- "VEGFR2 tip cell specific expression"
- "DLL4 tip cell stalk cell difference"
- "ESM1 tip cell marker VEGF"
```

---

## PHASE 4: 사용 가능한 시약 종합

### Step 4.1: Genetic Tools

**4.1.1 KO Mouse Lines (JAX, EMMA, etc.)**
```
검색: Mouse Genome Informatics (MGI)

[Gene] 관련 마우스 라인:
| Allele | Type | Promoter | Availability | Stock # | Notes |
|--------|------|----------|--------------|---------|-------|
| Kdr tm1Jrt | Null | - | JAX | 002938 | Embryonic lethal |
| Kdr tm2.1Jrt | Floxed | - | JAX | 012345 | Conditional ready |
| Cdh5-CreERT2 | Cre | VE-Cadherin | JAX | 012345 | EC-specific, tamoxifen |
```

**4.1.2 siRNA/shRNA**
```
Validated sequences from literature:

[Gene] siRNA:
| Sequence (5'-3') | Target Region | Knockdown % | Cell Type | PMID |
|------------------|---------------|-------------|-----------|------|
| GGAAUAUUCUGUUAUC... | CDS 1234-1254 | 85% | HUVEC | xxx |
| CCUUAAGGAAUCUUGG... | 3'UTR 2345-2365 | 78% | mESC | xxx |

Commercial sources:
- Dharmacon ON-TARGETplus: Cat# L-003xxx
- Sigma MISSION: Cat# SASI_Hs01_xxx
- Thermo Silencer Select: Cat# s1234
```

**4.1.3 CRISPR Guide RNAs**
```
Validated gRNA sequences:

[Gene] CRISPR:
| gRNA Sequence | PAM | Exon | Efficiency | Off-target Score | Source |
|---------------|-----|------|------------|------------------|--------|
| ACGTACGT... | NGG | 3 | 89% | 92 | Broad GPP |
| TGCATGCA... | NGG | 5 | 85% | 95 | CRISPRscan |
```

### Step 4.2: Pharmacological Tools

**4.2.1 Small Molecule Inhibitors**

ChEMBL/DrugBank 조회:
```
[타겟] Inhibitors:

| Compound | Type | IC50 | Selectivity | In vivo dose | Route | PMID |
|----------|------|------|-------------|--------------|-------|------|
| Axitinib | TKI | 0.1nM | VEGFR1/2/3 | 25mg/kg bid | PO | xxx |
| Tivozanib | TKI | 0.21nM | VEGFR1/2/3 | 1mg/kg qd | PO | xxx |
| SU5416 | TKI | 1.2μM | VEGFR2 | 20mg/kg | IP | xxx |

Commonly used in vitro concentrations:
| Compound | Low | Medium | High | Duration | Notes |
|----------|-----|--------|------|----------|-------|
| DAPT | 1μM | 10μM | 25μM | 24-72h | γ-secretase |
| SU5416 | 1μM | 5μM | 10μM | 24-48h | VEGFR2 |
```

**4.2.2 Blocking Antibodies**

```
[타겟] Blocking Antibodies:

| Antibody | Clone | Species | Blocking? | In vivo dose | Vendor | Cat# | PMID |
|----------|-------|---------|-----------|--------------|--------|------|------|
| anti-VEGFR2 | DC101 | Mouse | Yes | 40mg/kg 2x/wk | BioXCell | BE0060 | xxx |
| anti-DLL4 | YW152F | Human | Yes | 10mg/kg | - | - | xxx |
| anti-Notch1 | - | Human | Yes | 5mg/kg | - | - | xxx |

For in vitro:
| Antibody | Working conc. | Application | Vendor | Cat# |
|----------|---------------|-------------|--------|------|
| anti-VEGFR2 | 10μg/mL | Blocking | R&D | MAB3571 |
```

**4.2.3 Recombinant Proteins**

```
[Pathway] Recombinant Proteins:

| Protein | Form | Species | Working conc. | Vendor | Cat# |
|---------|------|---------|---------------|--------|------|
| VEGF-A165 | Carrier-free | Human | 10-50ng/mL | R&D | 293-VE |
| VEGF-A164 | Carrier-free | Mouse | 10-50ng/mL | R&D | 493-MV |
| DLL4-Fc | Fc chimera | Human | 1-5μg/mL | R&D | 1506-D4 |
| Ang-1 | His-tag | Human | 100-500ng/mL | R&D | 923-AN |
```

### Step 4.3: Detection/Analysis Tools

**4.3.1 Antibodies for Detection**

```
[Pathway] Detection Antibodies:

| Target | Clone | Application | Dilution | Vendor | Cat# |
|--------|-------|-------------|----------|--------|------|
| VEGFR2 | 55B11 | WB, IHC | 1:1000, 1:200 | CST | 2479 |
| p-VEGFR2 (Y1175) | 19A10 | WB | 1:1000 | CST | 2478 |
| CD31/PECAM1 | MEC13.3 | IHC, Flow | 1:100 | BD | 553370 |
| ESM1 | - | IHC | 1:200 | R&D | AF1999 |
```

**4.3.2 Reporter Assays**

```
[Pathway] Reporter Systems:

| Reporter | Readout | Application | Source |
|----------|---------|-------------|--------|
| HRE-Luciferase | Hypoxia/HIF | Hypoxia response | Addgene #26731 |
| VEGFR2-GFP | VEGFR2 expression | Live imaging | Published |
| Notch reporter (CBF1) | Notch activity | Notch signaling | Addgene #26905 |
```

---

## PHASE 5: 실험 설계 매트릭스 생성

### Step 5.1: 종합 실험 매트릭스

```
[타겟: VEGFR2] × [표현형: Tip cell migration] 실험 매트릭스:

| Experiment | Perturbation | Control | Model | Readout | Expected Result | Priority | Reference |
|------------|--------------|---------|-------|---------|-----------------|----------|-----------|
| 1 | VEGFR2 siRNA | Scramble | HUVEC | Migration assay | ↓ Migration 50% | High | PMID:xxx |
| 2 | SU5416 5μM | DMSO | HUVEC | Sprouting | ↓ Sprouts 60% | High | PMID:xxx |
| 3 | Kdr EC-KO | WT | Mouse retina P5 | Isolectin | ↓ Tip cells | High | PMID:xxx |
| 4 | DC101 Ab | IgG | Mouse retina | Tip cell count | ↓ Tips 40% | Medium | PMID:xxx |
| 5 | VEGF-A + DAPT | VEGF-A alone | Organoid | Sprouting | ↑ Tips | Medium | PMID:xxx |
| 6 | VEGFR2 OE | Control | mESC EC | Migration | ↑ Migration | Low | - |
```

### Step 5.2: 비교 실험 설계

```
[VEGFR2 KO vs VEGFR3 KO 비교]:

| Parameter | VEGFR2 KO (Kdr-/-) | VEGFR3 KO (Flt4-/-) | Reference |
|-----------|---------------------|---------------------|-----------|
| Viability | Embryonic lethal E8.5 | Embryonic lethal E10.5 | xxx, xxx |
| EC-specific KO viability | Viable | Viable | xxx, xxx |
| Retinal vascular density | ↓↓ 60% | ↓ 20% (lymphatic) | xxx |
| Tip cell number | ↓↓ 70% | → No change | xxx |
| Stalk cell proliferation | ↓ 40% | → | xxx |
| Lymphatic development | → | ↓↓↓ | xxx |
| Primary defect | Angiogenesis | Lymphangiogenesis | - |
```

---

## PHASE 6: Crosstalk Diagram 생성

### Step 6.1: Diagram 요소 정의

```
Nodes (각 노드에 시약 정보 포함):
- VEGF-A [Recomb: R&D 293-VE]
- VEGFR2 [Inh: Axitinib, Ab: DC101, siRNA: available]
- PLCγ [No specific inhibitor]
- MAPK [Inh: U0126, PD98059]
- Notch1 [Inh: DAPT, Ab: anti-Notch1]
- DLL4 [Ab: anti-DLL4, Recomb: DLL4-Fc]
- HIF-1α [Inh: Echinomycin]

Edges:
- VEGF-A → VEGFR2: Activation (ligand binding)
- VEGFR2 → PLCγ: Activation (phosphorylation)
- VEGFR2 → DLL4: Upregulation (transcription)
- DLL4 → Notch1: Activation (juxtacrine)
- Notch1 → VEGFR2: Inhibition (transcription)
- Hypoxia → HIF-1α: Stabilization
- HIF-1α → VEGF-A: Upregulation
```

### Step 6.2: Scientific-schematics 프롬프트

```
Generate a pathway crosstalk diagram with:

1. Layout: Left-to-right signal flow
2. Compartments:
   - Extracellular (ligands)
   - Membrane (receptors)
   - Cytoplasm (signaling)
   - Nucleus (transcription)

3. Node styling:
   - Receptors: Rectangle with rounded corners
   - Ligands: Oval
   - Kinases: Diamond
   - Transcription factors: Hexagon

4. Edge styling:
   - Activation: Green arrow
   - Inhibition: Red bar-headed line
   - Transcriptional: Dashed line

5. Annotations:
   - Small boxes next to each node listing available reagents
   - Star symbol for our experimental targets
   - Question mark for unknown connections

6. Color coding:
   - VEGF pathway: Blue
   - Notch pathway: Orange
   - Hypoxia pathway: Purple
   - Novel findings: Red highlight
```

---

## PHASE 7: 최종 출력물

### Output 파일 구조:

```
experiment_design_output/
├── summary/
│   ├── executive_summary.md        # 1-2페이지 요약
│   └── key_experiments.md          # 우선순위 실험 목록
├── detailed_analysis/
│   ├── target_deep_dive.md         # 타겟 상세 분석
│   ├── pathway_crosstalk.md        # Crosstalk 분석
│   ├── ko_phenotypes.md            # KO 표현형 정리
│   └── literature_evidence.md      # 문헌 근거
├── reagents/
│   ├── genetic_tools.csv           # siRNA, CRISPR, KO lines
│   ├── small_molecules.csv         # 억제제, 활성제
│   ├── antibodies.csv              # 차단, 검출 항체
│   ├── recombinant_proteins.csv    # 재조합 단백질
│   └── vendor_order_list.xlsx      # 주문 목록 (vendor, cat#, price)
├── experiment_matrix/
│   ├── full_matrix.csv             # 전체 실험 매트릭스
│   ├── priority_experiments.csv    # 우선순위 실험
│   └── comparison_table.csv        # KO 비교 테이블
├── figures/
│   ├── crosstalk_diagram.png       # Pathway 다이어그램
│   ├── crosstalk_diagram.svg       # 편집 가능 버전
│   └── experimental_scheme.png     # 실험 설계도
└── references/
    ├── references.bib              # BibTeX
    └── pmid_list.txt               # PMID 목록
```

### Executive Summary 템플릿:

```markdown
# [타겟] Experiment Design Summary

## Target Overview
- **Gene**: [Gene symbol] ([Full name])
- **Function**: [One sentence]
- **Druggability**: [Score/Assessment]

## Key Pathway Interactions
1. [타겟] ↔ [Pathway 1]: [Relationship]
2. [타겟] ↔ [Pathway 2]: [Relationship]

## Recommended Experiments (Priority Order)

### Experiment 1: [Name]
- **Perturbation**: [What to do]
- **Control**: [Control group]
- **Model**: [Model system]
- **Readout**: [Measurement]
- **Expected**: [Expected result]
- **Reagents**: [Specific reagents with cat#]
- **Reference**: PMID:[xxx]

### Experiment 2: [Name]
...

## Available Reagents Summary
- **Inhibitors**: [N] compounds, best: [Name] (IC50: X)
- **Antibodies**: [N] blocking, [M] detection
- **Genetic**: [KO lines], [siRNA sequences]

## Key Literature
1. [Author et al., Year] - [One sentence finding]
2. ...

## Unresolved Questions
1. [Question 1]
2. [Question 2]

## Suggested Timeline
Week 1-2: [Experiments]
Week 3-4: [Experiments]
...
```

---

## 사용 예시

```
/experiment-design

입력:
- 타겟: VEGFR2 (KDR)
- 모델: mouse retina, hiPSC retinal organoid
- 표현형: tip cell migration, vascular sprouting
- 비교: VEGFR2 KO vs VEGFR3 KO vs WT
- 관심 crosstalk: Notch, Hippo, Wnt

출력:
→ 위의 모든 분석 + 파일 생성
```
