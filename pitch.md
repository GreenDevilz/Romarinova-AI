# Pitch and Impact

## 1. Impact and future work

### 1.1 Prioritized developable compound-variety candidates

Romarinova AI delivers prioritized developable compound-variety candidates with provenance and evidence from three axes: ligand-based ML screening, ligand-based similarity/developability scoring, and a receptor-based docking specification integrated into the scoring formula but not executed in this submission. Impact is quantified by the actual eligible rosemary compound N (28), shortlist N (10 pairs), and U (5 unique compounds).

### 1.2 Modest-shortlist framing

The locked shortlist (N = 10 pairs, U = 5 unique compounds) is composed entirely of small monoterpenes: p-cymene (rank 1), verbenone (rank 2), γ-terpinene (rank 3), α-terpinolene (rank 4), and α-pinene (rank 5). Every shortlisted compound passes all four developability checks (Lipinski, Veber, PAINS, aggregator).

The most rosemary-characteristic diterpenes — carnosic acid, carnosol, and rosmarinic acid — appear only in the diagnostics data. Carnosol and carnosic acid fail the PAINS filter; rosmarinic acid fails both Veber and PAINS. These are preserved in the diagnostics appendix with their failed-rule annotations for experimental follow-up.

This is prioritization discipline, not a failure of novelty. The pipeline separates developable candidates from large, PAINS-prone, low-bioavailability diterpenes and preserves both classes for downstream work.

### 1.3 ML delivers on-domain; the pipeline reports where its domain ends

**The model works on target chemistry.** Trained on 14 ChEMBL actives using class-aware 4-fold CV, the selected Random Forest achieves **PR-AUC 0.781 and ROC-AUC 0.760** on a document-held-out test split (7 actives, 14 inactives). This is a genuine +0.046 PR-AUC uplift over the naive similarity-only baseline (0.735). Small-sample caveats apply per protocol §5.3: the 95% confidence interval is wide `[0.435, 1.000]` and `ci_not_reported = True` because the active count falls below the 20-compound threshold.

**The pipeline's value is what it does at the domain boundary.** Screening the rosemary catalogue, all 28 compounds fall outside the applicability domain (ad_distance 0.8014–0.9265 vs. threshold 0.7330). The model was trained on drug-like scaffolds and evaluated on small monoterpenes with no close training neighbors. The consequence is visible in the data — **27 of 28 screened compounds (96.4%) receive the identical `ML_rank_score` of 0.633333**. The model produces essentially one prediction for the entire extrapolation regime.

**The side-by-side comparison view is the instrument that surfaces this.** The α-Glucosidase Screening page shows the locked ranking and the ML-independent ranking simultaneously. Pair-level overlap is **10/10** — the sensitivity label is **`robust`**. Every shortlisted pair survives removal of the ML weight.

The correct reading of `robust` here is: **ML was not load-bearing on this regime, so removing it changes nothing**. All five shortlisted compounds share `ML_rank_score = 0.633333`, `developability_score = 1.0`, and `novelty_bonus = 1.0`; the rank order is determined entirely by `similarity_support`. A pipeline that hid this by scaling the tie-band score into the ranking would appear more impressive and be less trustworthy.

### 1.4 Compound-level view

The compound-level view is the primary traceability artifact, with U = 5 reported. Every row carries `ligand_only_score`/`ligand_only_rank`, `primary_cultivar` (display-only) and `all_cultivars` (authoritative) both shown with a tooltip, and the exact-match flags.

`known_active_count = 0` — no rosemary compound matches a ChEMBL active InChIKey. `high_ml_ad_extrapolation_count = 0` reflects the flag's ≥0.80 threshold; all 5 shortlist compounds sit at `ML_rank_score = 0.633333`, just below the boundary. The honest collapse figure is reported separately as `ml_rank_top_band_share = 0.9643`.

### 1.5 Out of scope — future work only

The following are explicitly **not** delivered by this submission:

- **Docking execution.** The docking protocol is fully specified (`models/docking_config.json`), and the scoring specification reserves a `0.20 × docking_score` term. No docking poses have been run: `results/docking_results.csv` is header-only. Docking is the highest-priority extension.
- **MD simulation and free-energy calculations** — out of scope.
- **Extension to α-amylase, lipase, or other glycosidases** — out of scope.
- **Alternative fingerprints (ECFP6) and tautomer standardization** — out of scope.
- **Prospective experimental validation of shortlist candidates** — outside a computational hackathon; not modeled.

## 2. Pitch summary

Romarinova AI turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system. It links documented Tunisian rosemary chemistry to machine learning, a bounded evidence record, and a three-page web interface.

**The ML pipeline works on its target chemistry.** On a document-held-out test split — 128 train, 30 validation, 21 test compounds, with zero shared `document_chembl_id` across splits — the selected Random Forest achieves **PR-AUC 0.781 and ROC-AUC 0.760** on 7 actives and 14 inactives. A similarity-only baseline reaches 0.735, meaning the ML model provides a modest, honest +0.046 PR-AUC uplift.

**The pipeline reports where that success ends.** All 28 screened rosemary compounds fall outside the applicability domain. On this extrapolation regime, **27 of 28 compounds (96.4%) receive the identical `ML_rank_score` of 0.633333**. The ML-independent comparison view makes this visible: pair-level overlap is **10/10** (`robust`), which reads honestly as "ML was not load-bearing on this regime."

**The shortlist is defensible because the pipeline degrades gracefully.** The locked shortlist (N = 10 pairs, U = 5 compounds — p-cymene, verbenone, γ-terpinene, α-terpinolene, α-pinene) is carried by similarity support and developability. Characteristic diterpenes are preserved in diagnostics. `known_active_count = 0` — novelty is earned, not borrowed.

**The contribution** is variety- and chemotype-level traceability; scaffold-aware validation with honest small-sample reporting; developability filtering; a side-by-side ML-independent comparison view that makes the model's contribution boundary visible; and a reproducible open pipeline with locked rules, decision log, and compliance audit.

**The path forward is data, not model complexity.** A Random Forest trained on 14 actives cannot learn to distinguish small monoterpenes it has never seen. The two levers that would move the needle are (i) more independently sourced activity data, and (ii) docking as a second, training-data-independent axis of evidence.

Closing line: *Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that delivers on-target ML performance and reports where that performance ends.*