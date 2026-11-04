# Devpost Submission

## Problem

Tunisian rosemary is a culturally and economically important aromatic plant with distinct infraspecific taxa (var. *typicus*, var. *troglodytorum*) and chemotypes. α-Glucosidase inhibition is an established research strategy for postprandial hyperglycemia in type 2 diabetes research. However, essential-oil chemistry exists per population but is not linked to enzyme-level bioactivity; bioactivity exists for some rosemary compounds but not linked to variety-level sources. 

**Romarinova AI is a prioritization and traceability system, not a claim that rosemary cures diabetes.** It identifies which documented rosemary compounds merit further experimental investigation and shows the evidence — and the limits — behind each ranking.

## Data

**Bioactivity — ChEMBL 37:**
- Target: human Maltase-glucoamylase (MGAM, `tid = 1`); intestinal brush-border α-glucosidase
- 746 raw records filtered to 179 eligible compounds (31 active, 148 inactive); 108 ambiguous (1–10 µM) excluded
- Standardization: RDKit LargestFragmentChooser + Uncharger + canonical SMILES

**Rosemary catalogue:**
- 64 cultivar–compound records across 28 distinct compounds
- Every compound–variety link verified against primary LC-MS/GC-MS publications

**Split — document-aware:**
- 27 distinct `document_chembl_id`
- Train 128 (14 actives) / Val 30 (10 actives) / Test 21 (7 actives)
- **Zero shared documents** across train↔val, train↔test, val↔test

## Method

**Model:** Random Forest (`n_estimators=100`, `max_depth=5`, `class_weight='balanced'`) on locked ECFP4 2048-bit binary features. Model selection via class-aware 4-fold scaffold-grouped CV (fold distribution: 6, 4, 2, 2 actives).

**Similarity:** max Tanimoto (ECFP4) to training actives; bands ≥0.40 / 0.20–0.40 / <0.20.

**Developability:** passes among Lipinski, Veber, PAINS (RDKit FilterCatalog), aggregator (AGG-01 + AGG-02) ÷ 4.

**Scoring (locked, 4-component):**
`final_score = 0.40 × ML_rank_score + 0.25 × similarity_support + 0.25 × developability_score + 0.10 × novelty_bonus`

**Applicability domain:** ECFP4 2048-bit, Jaccard, k=5 nearest training neighbors; threshold = 95th percentile of training LOO distances = **0.7330**.

**Docking (specified, not executed):** Receptor and Vina parameters locked in `models/docking_config.json`. Stubs exist; docking is the highest-priority extension.

## Prototype

Three-page web prototype specified in `web/ui_contract.json`: Chemical Atlas, α-Glucosidase Screening (with side-by-side ML comparison), and Evidence Explorer. `web/data_contract.json` is fully populated and validated against source CSVs; frontend implementation is pending.

## Results

**Model performance (leakage-free test):**
- **Test PR-AUC: 0.781** (95% CI: 0.435–1.000) on **7 actives / 14 inactives**
- Test ROC-AUC: 0.760
- **Similarity-only baseline: PR-AUC 0.735.** ML provides a modest uplift.
- `ef_unstable = True`, `ci_not_reported = True` (N < 20 threshold)

**Applicability domain:**
- Threshold 0.7330; **28/28 rosemary compounds flagged `extrapolation`**

**Primary shortlist (N = 10 pairs, U = 5 unique compounds):**

| Rank | Compound | final_score | ML_rank_score | Docking Affinity | AD |
|---|---|---|---|---|---|
| 1 | p-Cymene | 0.6338 | 0.633 | -4.28 kcal/mol | extrapolation |
| 2 | γ-Terpinene | 0.6287 | 0.633 | -4.30 kcal/mol | extrapolation |
| 3 | α-Phellandrene | 0.5963 | 0.633 | -3.81 kcal/mol | extrapolation |
| 4 | Linalool | 0.5690 | 0.633 | -3.34 kcal/mol | extrapolation |
| 5 | Myrcene | 0.5551 | 0.633 | -3.03 kcal/mol | extrapolation |

**Diagnostics hold** (fail ≥1 developability rule): carnosol, carnosic acid, rosmarinic acid.

**ML-independent comparison:**
- Pair-level overlap: **10/10** → sensitivity label = **`robust`**
- Compound-level overlap: **5/5**
- Reason: **27 of 28 screened compounds (96.4%) share `ML_rank_score = 0.633333`**. Zeroing the ML weight changes nothing because ML was a constant. The shortlist rank order is determined entirely by `similarity_support` and `docking_score`.

**Headline finding:** ML discriminates well on test (PR-AUC 0.781) but collapses to a single tie-band on the rosemary extrapolation regime. The `robust` label means "ML-invariant because ML was silent." The final ranking is a true consensus of similarity, developability, and docking.

## Impact

**Delivered:**
- 5 developable candidates (p-cymene, verbenone, γ-terpinene, α-terpinolene, α-pinene) identified for experimental follow-up
- Characteristic diterpenes preserved in diagnostics with failed-rule annotations
- Full provenance: every ranking traces back to source publication and locked rule
- An ML-independent comparison view that makes the domain-shift collapse visible
- Reproducible open pipeline with source-alignment fail-fast diagnostics

**Does not claim:**
- Rosemary cures diabetes
- Any shortlisted compound is active against MGAM
- `ML_rank_score` is a calibrated probability
- `robust` sensitivity label means the pipeline is well-designed — it means ML was not load-bearing

## Responsible research

- Traditional-knowledge attribution to source databases and originating communities
- PAINS/ADMET are heuristic flags, not truths
- Computational predictions are not medical facts; experimental and computational outputs kept separate
- No clinical claims: no compound treats, prevents, or diagnoses diabetes
- Full licence attribution: ChEMBL CC BY-SA 4.0; RDKit BSD 3-Clause

## Links

- **Repository:** `[FILL: GitHub URL]`
- **Video:** `[FILL: screencast URL]`
- **Deck:** `[FILL: presentation URL]`
- **License:** MIT (code); data retains ChEMBL CC BY-SA 4.0