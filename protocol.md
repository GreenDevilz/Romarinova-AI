# Protocol — Romarinova AI (EthnoHack 2026 Track 2)

**Status:** Final locked in-event protocol. Committed before 2 Nov 2026.

**Protocol posture:** Single-path and bounded. No scope fallbacks or contingencies. All conditionals are pre-registered reporting/QC rules fixed before data freeze.

**Canonical artifacts:** `operating_card.md` (daily execution), `protocol.md` (this document: master definitions and full methodology), `scripts/build_artifacts.py` (single generator in-event). This document is authoritative for disputes and pre-registered rule text.

---

# 1. Final protocol

## 1.1 Scope

- **Target:** human α-glucosidase only. Non-human α-glucosidase data excluded.
- **Data:** ChEMBL human α-glucosidase actives/inactives; rosemary catalogue covering var. *typicus*, var. *troglodytorum*, and the predominant 1,8-cineole/camphor/α-pinene chemotype. Data freeze end of Day 1.
- **Model:** XGBoost and Random Forest on locked ECFP4 2048-bit binary; select by scaffold-grouped training CV; model freeze end of Day 2.
- **Computation universe:** all eligible rosemary compounds with `structure_status = defined`. AD, similarity, developability, `final_score`, and `ligand_only_score` computed for all eligible compounds before ranking lock.
- **Evaluation:** test-set headline with real human inactives; PR-AUC, ROC-AUC, EF@10%, top-20/top-N; scaffold-clustered bootstrap CIs with locked small-N behavior; validation PR-AUC reported separately.
- **AD:** training-LOO threshold primary; validation threshold reported. Computed before ranking lock. `ad_missing` flagged; AD does not modify `final_score`.
- **Shortlist:** up to 10 compound-variety pairs by `final_score` from eligible compounds with ≥1 identified pair; duplication cap two per compound; actual N reported.
- **Docking:** receptor-based docking of all eligible rosemary compounds into the locked human MGAM C-terminal catalytic site (PDB receptor, single locked grid box; methodology in §1.2.1). `docking_score` is a locked component of `final_score`.
- **Web prototype:** three pages — Chemical Atlas (dynamic), α-Glucosidase Screening (dynamic), Evidence Explorer (static generated summary with evidence drawer).
- **Out of scope:** no MD simulation; no free-energy calculations; no COCONUT sample; no external natural-product library screening; no automated literature scan; no manual extraction for activity evidence; no knowledge graph; no branch framework; no D overlay; no late-batch handling; no incomplete-input ranking fallback.

## 1.2 Fixed pipeline

1. **Target-based ML screening:** XGBoost and Random Forest on locked ECFP4; select by scaffold-grouped training CV; validate with scaffold splits; screen eligible rosemary catalogue.
2. **Ligand-based similarity annotation:** max Tanimoto to training actives; label known-scaffold analogs, medium-support compounds, novel scaffolds.
3. **Developability and PAINS filtering:** Lipinski/Veber, PAINS, aggregator, ADMET; compute locked priority score.
4. **Receptor-based docking:** AutoDock Vina docking of every eligible rosemary compound into the locked MGAM grid; `docking_score` enters `final_score` with a locked weight (details in §1.2.1 and §4.7).

No MD simulation. No free-energy calculations. No knowledge graph. No external library screening.

### 1.2.1 Docking protocol

- **Receptor.** Human MGAM C-terminal catalytic domain, PDB entry locked pre-event in `models/docking_config.json` (primary: PDB 3L4W; rejected candidates and reasons recorded in `compliance.md` §4). Remove crystallographic waters and non-catalytic heteroatoms; retain the catalytic Zn; add hydrogens at pH 7.4; assign Gasteiger charges; write PDBQT with Meeko.
- **Grid box.** Single box centered on the catalytic site; center and size locked in `models/docking_config.json`; one box for all compounds; no per-compound redocking box.
- **Ligands.** 3D conformer from the locked canonical SMILES: RDKit ETKDGv3 (`randomSeed = 42`), MMFF94-optimized, protonated at pH 7.4, converted to PDBQT with Meeko. One conformer per compound; no conformer enumeration.
- **Docking engine.** AutoDock Vina, `exhaustiveness = 32`, `num_modes = 9`, `seed = 42`. Single run per compound; top-scoring pose retained.
- **Failure handling.** No pose, non-finite affinity, or preparation failure → `docking_missing = True`, affinity left blank, `docking_score = 0`. The compound is **not** dropped from the shortlist pool.
- **Reproducibility.** Vina version, receptor PDB ID, grid box, seed, and exhaustiveness recorded in `models/docking_config.json` and `compliance.md` §4.

## 1.3 Reporting stance

- Every table caption and narrative claim reports actual N from §2 of this document.
- Shortlist is "up to 10"; if fewer than 10 verified identified pairs exist, report actual N and coverage gap without imputation.
- `ML_rank_score` is a validation-rank percentile, not a calibrated probability. Raw output is `ML_raw_score` and marked uncalibrated.
- A high `ML_rank_score` on an AD-extrapolation compound is not strong activity evidence; flag `high_ml_ad_extrapolation`.
- `literature_evidence_level = experimental` only when `known_active_flag = True`; otherwise `none`. This is a known-active evidence flag, not a literature review.
- PAINS and aggregator flags are heuristic, not truths, especially for natural products.
- Known actives are retained but excluded from novelty claims.
- Variety assignment is provenance-based; the project does not rank varieties against each other.
- Docking scores are scoring-function estimates, not binding free energies; they ignore receptor flexibility, explicit water, entropy, and membrane context. They are one axis of evidence, not a verdict.

## 1.4 Dual ranking view

The α-Glucosidase Screening page shows the locked ranking and the ML-independent ranking side by side and simultaneously. A separate toggle switches between the pair-level shortlist and the compound-level view. The ML-independent ranking is the same pipeline with 0.00 ML weight (0.00/0.30/0.30/0.10/0.30), computed for all eligible rosemary compounds. Docking participates in both views identically; the comparison isolates the ML contribution only. It is a comparison view only: no second shortlist, no second compound-level view, no second diagnostic. Divergence is reported in `results/diagnostics.csv`.

## 1.5 Roles

- **P1:** Data & ML/AD Lead — ChEMBL curation, standardization, species check, Day-1 ledger, feature matrix, scaffold splits, CV folds, XGB/RF training and selection, test/validation evaluation with small-N flags, AD for all eligible compounds, ligand-only ranking, first-pass similarity/developability/scoring execution, `build_artifacts.py` ownership, README/reproduction, final video assembly; **docking execution** (receptor/ligand preparation, Vina runs, docking_score generation) and `models/docking_config.json` ownership.
- **P2:** Cheminformatics & Consensus Lead — RDKit QC, PAINS/ADMET/aggregator with locked rules, review/ownership of similarity annotation and scoring outputs, duplication cap, compound-level view, sensitivity, diagnostics appendix; owns `compliance.md` decision-log section (≤24 rows); **docking review** (grid-box sanity, pose inspection for top shortlist entries, docking QC sign-off).
- **P3:** Literature, Evidence & Frontend Lead — rosemary evidence review, provenance, source-publication records, Evidence Explorer content, diagnostics appendix sections, three-page frontend, real-data integration from Day 2, screencast segments, deck assembly.

P1 is single technical decision-maker; P2 is deputy. Decision-log entry promptly; no committee votes.

## 1.6 Freeze and lock schedule

- **Data freeze:** end of Day 1.
- **Model freeze:** end of Day 2.
- **Frontend UI contract freeze:** end of Day 2.
- **Docking freeze:** end of Day 2 — receptor, grid box, ligand-preparation protocol, and Vina parameters fixed.
- **Ranking lock:** Day 3 — compound identities, `final_score` rankings, `docking_score`, duplication cap, AD flags, compound-level view, ML-independent ranking fixed.
- **Final record verification:** Day 4 — locked shortlist, compound-level view, diagnostics data, generated artifacts verified for consistency.

## 1.7 Reporting conditions

- **N = 0:** state N = 0 and coverage gap; sensitivity `not_applicable_NO`; demo focuses on pipeline, data feasibility ledger, and evidence discipline.
- **N = 1–2:** actual N; coverage gap stated.
- **N ≥ 3:** standard narrative.

## 1.8 Ethics and limitations

- Traditional-knowledge attribution to source databases and, where available, originating communities; no appropriation or decontextualization.
- Benefit-sharing and prior informed consent are beyond hackathon scope; this is a computational prioritization tool, not a claim on traditional knowledge.
- PAINS/ADMET are heuristic flags, not truths; natural products often trigger PAINS; diagnostics appendix preserves follow-up candidates.
- Computational predictions are not medical facts; no compound should be used without experimental validation and clinical trials.
- Known-compound critique: the ranking may mix known-scaffold analogs and novel/PAINS-flagged candidates; the pipeline does not force novelty.
- No clinical claims: no compound treats, prevents, or diagnoses diabetes or any metabolic disease.

## 1.9 Bottom line

One target, one feature representation, one receptor structure, one evaluation protocol, one scoring rule (with a side-by-side ML-independent comparison view), one primary shortlist definition (cap 10) plus one compound-level view, one three-student roster, one hybrid ligand- and structure-based pipeline (ligand-based ML screening + receptor-based docking; no MD), one fixed descoped scope per component, one timeline with buffer, one compliance posture, one reporting rule, one master definitions table. Actual N is reported everywhere; coverage gaps are stated honestly; small-N behavior is locked in advance.

---

# 2. Master definitions and consistency table

Canonical source for all N terms and flags. Every table caption and narrative claim reports the relevant N from this table. All header metadata fields are generated by `scripts/build_artifacts.py` from the single results store and are never hand-maintained.

## 2.1 Master definitions

| Term | Definition | Reported in |
|---|---|---|
| **N\*** | Actual primary-shortlist pair count. | Diagnostics data; shortlist header; every relevant narrative claim. |
| **U\*** | Unique compounds in the primary shortlist (= U). | Diagnostics data; compound-level view; shortlist header. |
| **Test N_actives / N_inactives** | Test-split counts after curation and species check. | Every evaluation table caption. |
| **`known_active_count`** | Unique compounds in the primary shortlist with `known_active_flag = True` (pair-level duplicates counted once). Descriptive. | Ranking-lock header block. |
| **`high_ml_ad_extrapolation_count`** | Unique compounds in the primary shortlist with `high_ml_ad_extrapolation = True`. Descriptive. | Ranking-lock header block. |
| **`ad_extrapolation_count`** | Unique compounds in the primary shortlist with `ad_flag = extrapolation` (excluding `ad_missing = True`). Descriptive. | Ranking-lock header block. |
| **`ad_missing_count`** | Unique compounds in the primary shortlist with `ad_missing = True`. Descriptive. | Ranking-lock header block. |
| **`docking_affinity_kcal_per_mol`** | Best Vina affinity (top pose). | Docking results CSV; compound-level view. |
| **`docking_score`** | Locked [0,1] transform of affinity (§4.7). | Scoring; `final_score`; ML-independent score. |
| **`docking_missing`** | `True` when no pose, non-finite affinity, or preparation failure; `docking_score = 0`. | Diagnostics; compound-level view. |
| **`docking_pose_status`** | `top_pose_only`, `failed`, or `skipped`. | Diagnostics. |

## 2.2 Flag glossary

| Flag | Definition |
|---|---|
| `ad_flag` | `within_AD` if `ad_distance ≤ ad_threshold`, else `extrapolation`. |
| `high_ml_ad_extrapolation` | `True` when `ad_flag = extrapolation` and `ML_rank_score ≥ 0.80`. |
| `ml_rank_saturated` | `True` on entries with `ML_rank_score = 1.0` when the saturation fraction across the full scoring universe exceeds 20%. |
| `known_active_flag` | `True` iff InChIKey matches any curated ChEMBL active after species check. |
| `ef_unstable` | As defined in §5.3. |
| `ci_underpowered` | As defined in §5.3. |
| `ci_not_reported` | As defined in §5.3. |
| `cultivar_evidence_level` | `identified` or `suspected`. Only `identified` is shortlist-eligible. |
| `literature_evidence_level` | `experimental` only when `known_active_flag = True`; otherwise `none`. |
| `structure_status` | `defined`, `ambiguous`, or `missing`. Only `defined` is screened. |
| `docking_missing` | `True` when docking produced no usable pose or a non-finite affinity; `docking_score` set to 0. |
| `docking_pose_status` | `top_pose_only` when a pose was retained; `failed` on engine/preparation failure; `skipped` for compounds outside the docking universe. |

## 2.3 Consistency rule

Every table caption and narrative claim reports the relevant N from this table. If a value is wrong, the generator is fixed and the artifact regenerated; hand-editing generated CSVs is prohibited.

---

# 3. Applicability domain (AD)

## 3.1 Descriptor

Exactly locked ECFP4 2048-bit binary matrix; no scaling; metric Jaccard; k = 5 nearest neighbors in deduplicated training set.

## 3.2 AD distance

ad_distance(x) = (1/5) × Σ_{y ∈ NN_5(x)} (1 - Jaccard(x,y))


Mean Jaccard distance to the 5 nearest training neighbors.

## 3.3 Thresholds

- **Primary threshold:** 95th percentile of the training leave-one-out distance distribution.
- **Secondary threshold:** 95th percentile of validation distances.
- Report 90th/95th; report fraction of validation actives within AD under both thresholds.
- `within_AD` if `ad_distance ≤ ad_threshold`, else `extrapolation`.

## 3.4 Flags and scoring

- AD does not modify `final_score`.
- Computed for all eligible rosemary compounds before the ranking lock.
- Compounds with missing AD are ranked and flagged `ad_missing = True`; missing AD does not drop a compound from the shortlist.
- `high_ml_ad_extrapolation = True` when `ad_flag = extrapolation` and `ML_rank_score ≥ 0.80`.

## 3.5 Reporting

Report AD flags, `ad_missing`, ligand-only ranking, and diagnostics data. A high ML score on an AD-extrapolation compound is not strong activity evidence.

---

# 4. Scoring rule

## 4.1 Locked score

known_active = exact_training_match or exact_test_match or known_active_flag
similarity_support = max(0, min(1, max_similarity))
novelty_bonus = 0 if known_active else (1 if max_similarity < 0.20 else 0)

final_score = 0.30 × ML_rank_score
            + 0.20 × similarity_support
            + 0.20 × developability_score
            + 0.10 × novelty_bonus
            + 0.20 × docking_score

- `ML_rank_score` = validation-rank percentile clipped [0,1]; not a calibrated probability.
- `developability_score` = passes among Lipinski/Veber/PAINS/aggregator ÷ 4.
- `docking_score` = locked [0,1] transform of `docking_affinity_kcal_per_mol` (§4.7); `0` when `docking_missing = True`.
- Weights fixed a priori, not tuned on results.

## 4.2 ML-independent comparison score

ligand_only_score = 0.00 × ML_rank_score
                  + 0.30 × similarity_support
                  + 0.30 × developability_score
                  + 0.10 × novelty_bonus
                  + 0.30 × docking_score

Computed by P1 for all eligible rosemary compounds and stored alongside `final_score`. Comparison view only; it does not generate a second shortlist or a second compound-level view. The field name `ligand_only_score` is retained for artifact compatibility; in prose it is described as the ML-independent comparison score.

## 4.3 Computation universe

`final_score` and `ligand_only_score` computed for all eligible rosemary compounds (`structure_status = defined`). The shortlist is drawn from all eligible rosemary compounds with ≥1 identified pair and `structure_status = defined`.

## 4.4 Ties

Higher `similarity_support`, then higher `developability_score`, then lexicographic InChIKey, then lexicographic `cultivar_source`. Set `high_ml_ad_extrapolation = True` on entries with `ad_flag = extrapolation` and `ML_rank_score ≥ 0.80`. Exact known actives retained but zero novelty bonus.

## 4.5 Shortlist

Up to 10 compound-variety pairs by `final_score` from all eligible rosemary compounds with ≥1 identified pair, restricted to `cultivar_evidence_level = identified` and `structure_status = defined`. Actual N reported. Duplication cap: one compound appears at most twice; when a compound has more than two eligible variety pairs, keep the top 2 varieties by the locked tie-break. Because `final_score` has no variety-specific component, the operative selector is lexicographic `cultivar_source`. This is an arbitrary provenance rule, not a biological ranking.

## 4.6 Sensitivity check

Scenarios: locked weights 0.40/0.25/0.25/0.10; ligand-only 0.00/0.60/0.30/0.10. The ligand-only scenario is the sole sensitivity comparison. If N* = 0, report `not_applicable_NO`. Otherwise labels are computed on pair-level overlap using N* thresholds: robust (≥ [0.67 × N]), partially_robust (≥ [0.50 × N]), weight_sensitive (≤ [0.50 × N]). Pair-level overlap = number of (compound, variety) pairs in both primary and ligand-only shortlists. Report raw pair-level overlap, label, N, and U*.

## 4.7 Docking score

docking_affinity_kcal_per_mol = best Vina affinity (top pose; kcal/mol)

docking_score = clip(((-docking_affinity_kcal_per_mol) - 5.0) / 6.0, 0, 1)

Mapping: -5 kcal/mol → 0.0; -11 kcal/mol → 1.0; stronger binders saturate at 1.0; weaker binders floor at 0.0. The transform is fixed pre-event and is not tuned on results. `docking_missing = True` forces `docking_score = 0` and leaves the affinity blank. Docking is computed for the full scoring universe (§4.3) before the ranking lock.

---

# 5. Evaluation protocol

## 5.1 Evaluation set

- **Headline:** test split (15% scaffold-held-out), real human ChEMBL inactives (≥10 µM, confidence ≥7).
- **Validation transparency:** validation split, same criteria; used for model selection and `ML_rank_score` mapping.
- No decoys; real inactives only.

## 5.2 Locked evaluation table

| Model | PR-AUC (test, real inactives, 95% CI or flag) | ROC-AUC | EF@10% or alternative | Top-20/top-N recovery |
|---|---|---|---|---|
| Dummy | ... | ... | ... | ... |
| Similarity-only | ... | ... | ... | ... |
| Selected Stage 1 (XGB or RF) | ... | ... | ... | ... |
| Docking-only | ... | ... | ... | ... |

**Caption requirement:** every evaluation table caption reports test N actives, test N inactives, number of scaffold groups, and largest scaffold group size.

Validation PR-AUC for XGB and RF is reported separately as rows with `split = validation` in `results/evaluation_table.csv`.

## 5.3 Small-N metric behavior (thresholds on test N_actives)

| Metric | Condition | Behavior | Flag |
|---|---|---|---|
| Top-20 recovery | N_actives < 20 | Top-N with N = N_actives; state N. | `top_recovery_mode="top_N"` |
| EF@10% | 5 ≤ N_actives < 20, or top-10% bin < 5 positives | Report raw EF and bin-positive count; unstable. | `ef_unstable=True` |
| EF@10% | N_actives < 5 | Do not report EF; report raw precision and positives. | `ef_not_reported=True` |
| Bootstrap 95% CI | 20 ≤ N_actives < 50 | Report with underpowered warning. | `ci_underpowered=True` |
| Bootstrap 95% CI | N_actives < 20 | Do not report CI; raw metric only. | `ci_not_reported=True` |

All flags appear in the evaluation table, `protocol.md`, frontend, Devpost, and video. Narrative never reports a metric without N and flag state.

## 5.4 Leakage-safe rules

- Scaffold split first; save indices.
- Cross-split InChIKey check.
- Similarity baseline uses training actives only.
- Identical test set for dummy, similarity-only, XGB, RF, selected model.
- Scaffold-grouped CV for selection.
- Exclude compounds sharing `document_chembl_id` with test compounds from training.
- 1000 scaffold-clustered bootstrap resamples; report CI of difference.
- Report random-vs-scaffold Δ PR-AUC; flag >0.15 as scaffold-memorization risk.
- Docking receptor is not derived from any test-split compound; the receptor is the locked PDB structure, not a co-crystal of a rosemary compound.
- Docking-only baseline uses the same test set as all other baselines.

## 5.5 Literature validation

Limited to exact InChIKey known-active flags from curated ChEMBL data and to rosemary source-publication provenance. `literature_evidence_level = experimental` only when `known_active_flag = True`; otherwise `none`. No automated database scan, no manual extraction for activity evidence. The project assesses prioritization plausibility through exact known-active overlap and rosemary provenance, not through a new literature search.