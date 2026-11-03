# Presentation

## 1. Deck outline

### Slide 1: Title

**Romarinova AI**
*From Tunisian Rosemary Chemotypes to Evidence-Guided Candidate Prioritization for α-Glucosidase Inhibition*

EthnoHack 2026 · Track 2 — Natural Products and Drug Discovery
Three-student team · 2–8 Nov 2026
*"A prioritization and traceability system, not a claim that rosemary cures diabetes."*

### Slide 2: Problem and motivation

Tunisian rosemary is represented by distinct infraspecific taxa and chemotypes: **var. typicus** (1,8-cineole-rich), **var. troglodytorum** (camphor-rich), and the predominant **1,8-cineole/camphor/α-pinene** chemotype.

α-Glucosidase inhibition is an established research strategy for postprandial hyperglycemia in type 2 diabetes research.

**The fragmentation problem:**
- Essential-oil chemistry exists per population, not linked to enzyme-level bioactivity
- Bioactivity exists for some rosemary compounds, not linked to Tunisian variety/chemotype sources
- Computational rosemary-diabetes studies rarely provide variety-level traceability or evidence grading

**Our contribution:** variety- and chemotype-level traceability, scaffold-aware validation with honest uncertainty, developability filtering, an ML-independent comparison view, and a reproducible open pipeline.

### Slide 3: Data and approach

**Bioactivity data — ChEMBL 37:**
- Target: human Maltase-glucoamylase (MGAM, `tid = 1`); intestinal brush-border α-glucosidase
- 746 raw records → 179 eligible (31 active, 148 inactive); 108 ambiguous (1–10 µM) excluded, 2 activity-conflicts excluded
- Standardization: RDKit LargestFragmentChooser + Uncharger + canonical SMILES

**Rosemary catalogue:**
- 64 cultivar–compound records across 28 distinct compounds
- Every compound–variety link verified against primary LC-MS/GC-MS publications

**Split — document-aware (protocol §5.4):**
- 27 distinct `document_chembl_id` (largest: CHEMBL5244222, 26 compounds, 0 actives)
- Train 128 (14 actives) / Val 30 (10 actives) / Test 21 (7 actives)
- **Zero shared documents** across train↔val, train↔test, val↔test

### Slide 4: ML screening and evaluation

**Model:** Random Forest, `n_estimators=100`, `max_depth=5`, `class_weight='balanced'`, `random_state=42`
**Features:** ECFP4 2048-bit binary (locked); scaffold-grouped CV within training split

| Metric | Test | Validation |
|---|---|---|
| PR-AUC | **0.781** | 0.981 |
| ROC-AUC | **0.760** | 0.990 |
| N actives | 7 | 10 |
| N inactives | 14 | 20 |

**Small-N flags (protocol §5.3):** `ef_unstable = True`, `ef_not_reported = False`, `ci_not_reported = True`
*No 95% CI is reported because active count is below the 20-compound threshold.*

Baseline (class fraction): PR-AUC 0.333. Corrected model is **2.3× above baseline** on a leakage-free 21-compound test set.

**Validation caveat:** PR-AUC 0.981 is inflated by one document (CHEMBL5500361) contributing 9 of 10 validation actives. Test is the honest estimate.

### Slide 5: Similarity, developability, and scoring

**Similarity:** max Tanimoto (ECFP4) to training actives. Bands: ≥0.40 high / 0.20–0.40 medium / <0.20 novel scaffold.

**Developability** (passes among 4 rules ÷ 4):
- Lipinski (MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10)
- Veber (rotatable bonds ≤ 10, PSA ≤ 140 Å²)
- PAINS (RDKit FilterCatalog PAINS_A/B/C)
- Aggregator (AGG-01 long aliphatic chain + AGG-02 amphiphilic large)

**Locked 4-component score:**
final_score = 0.40 × ML_rank_score
              0.25 × similarity_support
              0.25 × developability_score
              0.10 × novelty_bonus


**ML-independent comparison:** `ligand_only_score` with ML weight 0.00.

### Slide 6: Shortlist and compound-level view

**Primary shortlist (N = 10 pairs, U = 5 unique compounds):**

| Rank | Compound | final_score | ML_rank_score | AD | Developability |
|---|---|---|---|---|---|
| 1 | Borneol | 0.707164 | 0.777778 | extrapolation | 1.00 |
| 2 | α-Eudesmol | 0.701927 | 0.777778 | extrapolation | 1.00 |
| 3 | Camphene | 0.699573 | 0.777778 | extrapolation | 1.00 |
| 5 | β-Pinene | 0.698611 | 0.777778 | extrapolation | 1.00 |
| 6 | Camphor | 0.698611 | 0.777778 | extrapolation | 1.00 |

**Diagnostics hold** (all fail ≥1 developability rule): carnosol (PAINS), carnosic acid (PAINS), rosmarinic acid (Veber + PAINS).

`known_active_count = 0` — no rosemary compound matches a ChEMBL active InChIKey.

### Slide 7: Side-by-side locked vs ML-independent comparison

**Pair-level overlap:** 6/10 → `partially_robust`
**Compound-level overlap:** 4/5
**ML-independent weights:** 0.00 / 0.30 / 0.30 / 0.10 (docking weight reserved but not executed)

**The finding:** all five shortlisted compounds share `ML_rank_score = 0.777778`. A Random Forest trained on 14 actives cannot discriminate among monoterpenes that fall into the same "no tree path fires" leaf cluster. **ML discriminates on test but collapses on the extrapolation regime.**

The shortlist ordering is therefore driven entirely by similarity and developability. The comparison view makes this visible rather than hiding it — the domain-shift caveat rendered as a UI element.

### Slide 8: Diagnostics and limitations

**Applicability domain:**
- Threshold: 0.7330 (95th percentile of training LOO distances)
- **28/28 rosemary compounds flagged `extrapolation`** (distances 0.8014–0.9265)
- `high_ml_ad_extrapolation_count = 0` — but this reflects the flag's ≥ 0.80 threshold, not absence of risk: all 5 shortlist compounds sit at 0.777778, just below the boundary

**Known gaps:**
- Docking: protocol specified, not executed (header-only `docking_results.csv`)
- Web prototype: schemas populated, pages not implemented
- Video/deck: pending submission deadline

**Method limitations:**
- Human MGAM-specific; other glycosidases out of scope
- ECFP4 only; no ECFP6 or tautomer standardization
- PAINS/ADMET are heuristic flags, not ground truth

### Slide 9: Impact and responsible research

**What the pipeline delivers:**
- Reproducible prioritization with provenance to source publications
- 5 developable candidates (borneol, α-eudesmol, camphene, β-pinene, camphor) identified for experimental follow-up
- Characteristic diterpenes preserved in diagnostics with failed-rule annotations

**Responsible research:**
- Traditional-knowledge attribution to source databases and, where available, originating communities
- No clinical claims: no compound treats, prevents, or diagnoses diabetes
- PAINS/ADMET are heuristic flags, not truths
- Experimental and computational outputs kept separate
- Full licence attribution (ChEMBL CC BY-SA 4.0; RDKit BSD 3-Clause; Vina/Meeko Apache 2.0)

### Slide 10: Conclusion and future work

**What was delivered:**
- A leakage-free, document-aware ML pipeline with honest small-N reporting
- A defensible shortlist driven by similarity and developability in the regime where ML collapses
- An ML-independent comparison view that makes the collapse visible
- A reproducible open repo with locked rules, decision log, and compliance audit

**Future work:**
- Execute docking (`prepare_receptor.py`, `dock_ligands.py` stubs → real Vina runs)
- Implement the three-page web prototype and demo screencast
- Prospective experimental validation of shortlist candidates
- Extend to α-amylase and other intestinal glycosidases

**Closing:** *Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that reports where its own ML contribution ends.*

## 2. Narration outline

1. **Opening / problem** — "Rosemary has been studied for α-glucosidase inhibition, but the data is fragmented across varieties, chemotypes, and assay papers. We built a system that traces every ranking back to its source."
2. **Data and approach** — "179 eligible ChEMBL compounds against human MGAM, split document-aware so no paper spans train and test. 28 rosemary compounds verified against primary LC-MS and GC-MS publications."
3. **Results / shortlist** — "Test PR-AUC 0.781 on 7 actives — 2.3× baseline. But `ci_not_reported = True`: 7 positives is too few for a confidence interval. Shortlist N = 10 pairs, U = 5 compounds."
4. **Live demo** — "Three pages: Chemical Atlas, α-Glucosidase Screening with side-by-side locked vs ML-independent view, and Evidence Explorer. Toggle between pair-level and compound-level views. The side-by-side comparison makes the ML collapse visible."
5. **Impact / future** — "Five developable candidates for experimental follow-up; characteristic diterpenes preserved in diagnostics with failed-rule annotations. Docking specified but not executed — highest-priority extension."
6. **Closing disclaimer** — "We do not claim rosemary cures diabetes. We identify which documented compounds merit experimental investigation, and we report where the model's contribution ends."

## 3. Narrative text

### 3.1 ML contribution statement

> "ML achieves a test PR-AUC of 0.781 on 7 actives and 14 inactives — strong on ChEMBL-like chemistry. But all five shortlisted rosemary compounds receive the identical `ML_rank_score` of 0.777778: the model cannot discriminate among monoterpenes that all fall into the same leaf cluster. On the extrapolation regime, ML adds **zero** ranking signal to the shortlist. The ML-independent comparison view (pair overlap 6/10, compound overlap 4/5) confirms this: the shortlist is preserved when ML weight is zeroed because ML was already silent."

### 3.2 Empty / short shortlist narrative

- **N = 0:** state N = 0 and coverage gap; sensitivity `not_applicable_NO`; demo focuses on pipeline, data feasibility ledger, and evidence discipline.
- **N = 1–2:** actual N; coverage gap stated.
- **N ≥ 3 (current state, N = 10):** standard narrative — shortlist of 10 pairs, U = 5, driven by similarity and developability.

### 3.3 Modest-shortlist framing

The pipeline separates developable candidates from large, PAINS-prone diterpenes and preserves the latter for experimental follow-up. Carnosol and carnosic acid fail PAINS; rosmarinic acid fails Veber and PAINS. The five shortlisted compounds are all small monoterpenes and sesquiterpenes. This is prioritization discipline, not a novelty failure — the developability filters are doing measurable work on the diterpene class, and the diagnostics data preserves them with their failed-rule annotations for future experimental follow-up.

### 3.4 ML / AD caveat

`ML_rank_score` is a validation-rank percentile, not a calibrated probability. A high ML rank on an AD-extrapolation compound is not strong activity evidence. All 28 rosemary compounds are outside the applicability domain (extrapolation at threshold 0.7330; distances 0.801–0.927), and `high_ml_ad_extrapolation_count = 0` reflects the flag's ≥ 0.80 boundary, not absence of risk — the actual figure that matters is `ad_extrapolation_count = 5`.

### 3.5 Closing disclaimer

> Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that reports where its own ML contribution ends.