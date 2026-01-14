# 🔄 Continuation Notes (다음 세션을 위한 정리)

## ✅ 완료된 작업

### 1. 세 가지 연구 도구 생성 완료
| 파일 | 용도 | 상태 |
|------|------|------|
| `find-scrna-data.md` | scRNA 공개 데이터 탐색 | ✅ 완료 |
| `target-discovery.md` | 타겟 발굴 파이프라인 | ✅ 완료 |
| `experiment-design.md` | Deep dive 실험 설계 | ✅ 완료 |
| `README.md` | 사용 가이드 | ✅ 완료 |

### 2. Git 상태
- Branch: `claude/explain-scientific-skills-dy9Gw`
- Commit: `7c84476` - "Add vascular biology research tools for scRNA-seq analysis"
- Push: ✅ 완료

### 3. 파일 위치
```
/home/user/claude-scientific-skills/
├── .claude/commands/           # 원본 (gitignore됨)
│   ├── find-scrna-data.md
│   ├── target-discovery.md
│   ├── experiment-design.md
│   └── README.md
└── docs/vascular-research-tools/  # Git에 포함된 복사본
    ├── find-scrna-data.md
    ├── target-discovery.md
    ├── experiment-design.md
    └── README.md
```

---

## 📋 다음에 할 작업 (우선순위 순)

### 1. 실제 테스트 필요
- [ ] 샘플 scRNA 데이터로 `/target-discovery` 테스트
- [ ] VEGFR2로 `/experiment-design` 테스트
- [ ] CellxGene Census 연결 테스트

### 2. 개선 사항
- [ ] 에러 핸들링 추가 (API 실패 시)
- [ ] 진행 상황 표시 추가
- [ ] 출력 파일 자동 저장 경로 설정

### 3. 추가 기능 고려
- [ ] `/compare-datasets` - 여러 scRNA 데이터셋 비교
- [ ] `/pathway-explorer` - 특정 pathway 심층 탐색
- [ ] `/reagent-finder` - 시약 검색 특화

### 4. 문서화
- [ ] 실제 사용 예시 추가 (스크린샷 포함)
- [ ] 랩미팅 데모 슬라이드 준비

---

## 🎯 사용자 연구 맥락

**연구 분야:**
- VEGF / VEGFR signaling
- Mouse retina vasculature
- mESC organoid
- hiPSC organoid
- Tip cell / Stalk cell biology
- Angiogenesis

**사용자 요구사항 핵심:**
1. 단순 검색이 아닌 **통합 파이프라인**
2. "upregulation/downregulation" 수준이 아닌 **구체적 실험 조건**
3. KO mouse 비교, 특정 세포 유형, 특정 표현형
4. **시약 정보까지** (항체, 억제제, siRNA 등)
5. **Nature Reviews급 Figure** 생성

**금요일 랩미팅 데모 예정**
- 교수님이 AI에 대해 잘 아심
- 얕은 데모는 소용없음
- 실제 연구에 도움되는 것 필요

---

## 💡 기억해야 할 포인트

1. **scRNA 데이터 위치 문제**: 한 곳에 없음
   - CellxGene Census: 정제된 atlas (API 가능)
   - GEO: 가장 많음 (전처리 필요한 경우 많음)
   - 논문: Data Availability 섹션에서 추출

2. **차별화 포인트**:
   - 여러 DB 통합 조회
   - 자연어로 복잡한 분석 지시
   - 시약 정보까지 한번에

3. **.gitignore 이슈**: `.claude` 폴더가 ignore됨
   - 해결: `docs/vascular-research-tools/`에 복사본 저장

---

## 🔗 관련 리소스

- CellxGene Census: 61M+ cells
- STRING: protein-protein interactions
- KEGG/Reactome: pathway databases
- DrugBank/ChEMBL: drug/compound info
- PubMed: literature search

---

## ⏰ 마지막 업데이트
- 날짜: 2026-01-14
- 세션: claude/explain-scientific-skills-dy9Gw
