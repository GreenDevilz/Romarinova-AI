# Devpost Submission

## Problem

Tunisian rosemary (*Rosmarinus officinalis* L., syn. *Salvia rosmarinus* Spenn.) is a culturally and economically important aromatic plant with distinct infraspecific taxa (var. *typicus*, var. *troglodytorum*) and chemotypes. α-Glucosidase inhibition is an established research strategy for postprandial hyperglycemia in type 2 diabetes research, and rosemary extracts and constituents have been studied in this context.

But the data is fragmented: essential-oil chemistry exists per population but is not linked to enzyme-level bioactivity; bioactivity exists for some rosemary compounds but not linked to Tunisian variety/chemotype sources; and computational rosemary-diabetes studies rarely provide variety-level traceability or explicit evidence grading.

**Romarinova AI is a prioritization and traceability system, not a claim that rosemary cures diabetes.** It identifies which documented rosemary compounds merit further experimental investigation and shows the evidence — and the limits — behind each ranking.

## Data

**Bioactivity — ChEMBL 37:**
- Target: human Maltase-glucoamylase (MGAM, `tid = 1`); intestinal brush-border α-glucosidase
- 746 raw records filtered to 179 eligible compounds (31 active, 148 inactive); 108 ambiguous (1–10 µM) excluded; 2 activity-conflicts excluded; non-human, wrong-assay-type, and low-confidence records counted and excluded
- Standardization: RDKit LargestFragmentChooser + Uncharger + canonical SMILES + InChIKey

**Rosemary catalogue:**
- 64 cultivar–compound records across 28 distinct compounds
- Every compound–variety link verified against primary LC-MS/GC-MS publications
- Structure source: PubChem (priority 1) → ChEMBL; all 28 defined
- `known_active_overlap_n = 0` — no rosemary compound matches a ChEMBL active InChIKey

**Split — document-aware:**
- 27 distinct `document_chembl_id`; largest is CHEMBL5244222 (26 compounds, 14.5%, 0 actives — benign)
- Train 128 (14 actives) / Val 30 (10 actives) / Test 21 (7 actives)
- **Zero shared documents** across train↔val, train↔test, val↔test
- Full ledger in `results/diagnostics.csv` (`section = data_feasibility_ledger`)

## Method

**Model:** Random Forest (`n_estimators=100`, `max_depth=5`, `class_weight='balanced'`, `random_state=42`) on locked ECFP4 2048-bit binary features; model selection by scaffold-grouped CV within training; validation PR-AUC used to define `ML_rank_score`.

**Similarity:** max Tanimoto (ECFP4) to training actives; bands ≥0.40 / 0.20–0.40 / <0.20.

**Developability:** passes among Lipinski, Veber, PAINS (RDKit FilterCatalog PAINS_A/B/C), aggregator (AGG-01 + AGG-02) ÷ 4.

**Scoring (locked):**
final_score = 0.40 × ML_rank_score
              0.25 × similarity_support
              0.25 × developability_score
              0.10 × novelty_bonus


The ML-independent comparison score sets ML weight to 0.00 (0.00 / 0.60 / 0.30 / 0.10 in the current 4-component configuration).

**Applicability domain:** ECFP4 2048-bit, Jaccard, k=5 nearest training neighbors; threshold = 95th percentile of training LOO distances = **0.7330**.

**Docking (specified, not executed):** `models/docking_config.json` locks receptor PDB 3L4W (human MGAM C-terminal catalytic domain), Vina parameters, and the `docking_score` transform. `scripts/prepare_receptor.py` and `scripts/dock_ligands.py` are stubs; `results/docking_results.csv` is header-only. Docking is the highest-priority extension.

## Prototype

Three-page web prototype specified in `web/ui_contract.json`:

1. **Chemical Atlas** — rosemary compound catalogue with structure, variety/varieties, evidence level, provenance
2. **α-Glucosidase Screening** — ranked compound–variety pairs with side-by-side locked vs ML-independent view; toggle between pair-level shortlist and compound-level view; AD flags with tooltips
3. **Evidence Explorer** — static generated summary with evidence drawer

**Current state:** `web/data_contract.json` and `web/ui_contract.json` are populated schemas; `web/api/` and `web/frontend/` contain `.gitkeep` files only. The pages are specified but not yet implemented.

## Results

**Model performance (leakage-free test):**
- **Test PR-AUC: 0.781** on **7 actives / 14 inactives**
- Test ROC-AUC: 0.760
- **`ef_unstable = True`, `ci_not_reported = True`** — no 95% CI reported because active count is below protocol §5.3's 20-compound threshold
- Baseline (class fraction): PR-AUC 0.333 → **2.3× above baseline**

**Applicability domain:**
- Threshold 0.7330; **28/28 rosemary compounds flagged `extrapolation`** (distances 0.8014–0.9265)

**Primary shortlist (N = 10 pairs, U = 5 unique compounds):**

| Rank | Compound | final_score | ML_rank_score | AD | Developability |
|---|---|---|---|---|---|
| 1 | Borneol | 0.707164 | 0.777778 | extrapolation | 1.00 |
| 2 | α-Eudesmol | 0.701927 | 0.777778 | extrapolation | 1.00 |
| 3 | Camphene | 0.699573 | 0.777778 | extrapolation | 1.00 |
| 5 | β-Pinene | 0.698611 | 0.777778 | extrapolation | 1.00 |
| 6 | Camphor | 0.698611 | 0.777778 | extrapolation | 1.00 |

**Diagnostics hold** (all fail ≥1 developability rule): carnosol (PAINS), carnosic acid (PAINS), rosmarinic acid (Veber + PAINS).

**ML-independent comparison:**
- Pair-level overlap: **6/10** → `partially_robust`
- Compound-level overlap: **4/5**
- Reason: all 5 shortlisted compounds share `ML_rank_score = 0.777778` — ML cannot reorder them

**Sensitivity label:** `partially_robust` (not `not_applicable_NO`)

**Headline finding:** ML discriminates well on test (PR-AUC 0.781) but collapses to a single tie-band on the rosemary extrapolation regime. The pipeline reports the collapse and preserves a defensible shortlist driven by similarity and developability.

## Impact

**Delivered:**
- 5 developable candidates (borneol, α-eudesmol, camphene, β-pinene, camphor) identified for experimental follow-up
- Characteristic diterpenes (carnosol, carnosic acid, rosmarinic acid) preserved in diagnostics with failed-rule annotations
- Full provenance: every ranking traces back to source publication, ChEMBL assay, and locked rule
- An ML-independent comparison view that makes the domain-shift collapse visible rather than hiding it
- Reproducible open pipeline with locked rules, decision log, and compliance audit

**Quantified:**
- Eligible rosemary compounds: 28
- Shortlist pairs: 10
- Unique shortlist compounds (U): 5
- Known-active overlaps: 0
- AD extrapolations in shortlist: 5/5
- ML test PR-AUC: 0.781 (2.3× baseline)

**Does not claim:**
- Rosemary cures diabetes
- Any shortlisted compound is active against MGAM
- Validation PR-AUC 0.981 generalizes (inflated by single-document active cluster)
- `ML_rank_score` is a calibrated probability

## Responsible research

- Traditional-knowledge attribution to source databases and, where available, originating communities; no appropriation or decontextualization
- Benefit-sharing and prior informed consent are ethical requirements beyond hackathon scope; this is a computational prioritization tool, not a claim on traditional knowledge
- PAINS/ADMET are heuristic flags, not truths; natural products often trigger PAINS; the diagnostics data preserves follow-up candidates
- Computational predictions are not medical facts; no compound should be used without experimental validation and clinical trials; experimental and computational outputs kept separate
- Known-compound critique: the ranking may naturally mix known-scaffold analogs and novel/PAINS-flagged candidates; the pipeline does not force novelty; the contribution is prioritization and transparency
- No clinical claims: no compound treats, prevents, or diagnoses diabetes or any metabolic disease; all outputs are exploratory hypotheses
- Full licence attribution: ChEMBL CC BY-SA 4.0; RDKit BSD 3-Clause; AutoDock Vina Apache 2.0; Meeko Apache 2.0; PDB public domain; all other libraries attributed in `README.md`

## Links

- **Repository:** `[FILL: GitHub URL]`
- **Video:** `[FILL: screencast URL]`
- **Deck:** `[FILL: presentation URL]`
- **License:** MIT (code); data retains ChEMBL CC BY-SA 4.0