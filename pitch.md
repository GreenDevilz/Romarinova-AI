# Pitch and Impact

## 1. Impact and future work

### 1.1 Prioritized developable compound-variety candidates

Romarinova AI delivers prioritized developable compound-variety candidates with provenance and evidence from three axes: ligand-based ML screening, ligand-based similarity/developability scoring, and a receptor-based docking specification integrated into the scoring formula but not executed in this submission. Impact is quantified by the actual eligible rosemary compound N (28), shortlist N (10 pairs), and U (5 unique compounds).

### 1.2 Modest-shortlist framing

The locked shortlist (N = 10 pairs, U = 5 unique compounds) is composed entirely of small monoterpenes: p-cymene (rank 1), verbenone (rank 2), γ-terpinene (rank 3), α-terpinolene (rank 4), and α-pinene (rank 5). Every shortlisted compound passes all four developability checks (Lipinski, Veber, PAINS, aggregator).

The most rosemary-characteristic diterpenes — carnosic acid, carnosol, and rosmarinic acid — appear only in the diagnostics data. Carnosol and carnosic acid fail the PAINS filter; rosmarinic acid fails both Veber and PAINS. These are not rejected for lack of interest — they are preserved in the diagnostics appendix with their failed-rule annotations, alongside `final_score` and provenance, so experimental follow-up can pick them up on its own terms.

This is prioritization discipline, not a failure of novelty. The pipeline separates developable candidates from large, PAINS-prone, low-bioavailability diterpenes and preserves both classes for downstream work.

### 1.3 ML delivers on-domain; the pipeline reports where its domain ends

**The model works on target chemistry.** Trained on 14 ChEMBL actives, evaluated on a document-held-out test split (7 actives, 14 inactives), the selected Random Forest achieves **PR-AUC 0.781 and ROC-AUC 0.760** — 2.3× the baseline positive-class fraction (0.333). Validation PR-AUC reaches 0.981, though this is inflated by a single-document active cluster (CHEMBL5500361 carries 9 of 10 validation actives) and is not the headline metric. Small-sample caveats apply per protocol §5.3: `ci_not_reported = True`, `ef_unstable = True` — the active count is below the 20-compound threshold for a 95% confidence interval.

**The pipeline's value is what it does at the domain boundary.** Screening the rosemary catalogue, all 28 compounds fall outside the applicability domain (ad_distance 0.8014–0.9265 vs. threshold 0.7330). This is expected: the model was trained on drug-like scaffolds and evaluated on small monoterpenes with no close training neighbors. The consequence is visible in the data — **27 of 28 screened compounds (96.4%) receive the identical `ML_rank_score` of 0.633333**. The model produces essentially one prediction for the entire extrapolation regime. Rather than scale that constant into the ranking and produce a shortlist that *looks* ML-driven but is actually noise, the pipeline detects the collapse and reports it.

**The side-by-side comparison view is the instrument that surfaces this.** The α-Glucosidase Screening page shows the locked ranking and the ML-independent ranking simultaneously. Pair-level overlap is **10/10** — the sensitivity label is **`robust`**. Every shortlisted pair survives removal of the ML weight.

The correct reading of `robust` here is not "the pipeline is robustly designed." It is: **ML was not load-bearing on this regime, so removing it changes nothing**. All five shortlisted compounds share `ML_rank_score = 0.633333`, `developability_score = 1.0`, and `novelty_bonus = 1.0`; the rank order among them is determined entirely by `similarity_support`. The sensitivity test measured this honestly and returned the honest answer. A pipeline that hid this — by scaling the tie-band score into the ranking, or by using a different weight that produced spurious order changes — would appear more impressive and be less trustworthy.

### 1.4 Compound-level view

The compound-level view is the primary traceability artifact, with U = 5 reported. Every row carries `ligand_only_score`/`ligand_only_rank`, `primary_cultivar` (display-only) and `all_cultivars` (authoritative) both shown with a tooltip, and the exact-match flags (`exact_training_match`, `exact_test_match`, `known_active_flag`).

`known_active_count = 0` — no rosemary compound matches a ChEMBL active InChIKey, so `novelty_bonus` is earned legitimately across the shortlist. `high_ml_ad_extrapolation_count = 0` reflects the flag's ≥0.80 threshold; all 5 shortlist compounds sit at `ML_rank_score = 0.633333`, just below the boundary. The honest collapse figure is reported separately as `ml_rank_top_band_share = 0.9643` (27/28 compounds at 0.633333), and `ad_extrapolation_count = 5` reports the true shortlist picture.

### 1.5 Out of scope — future work only

The following are explicitly **not** delivered by this submission:

- **Docking execution.** The docking protocol is fully specified (`models/docking_config.json` locks receptor PDB 3L4W, grid box, Vina parameters, and the `docking_score` transform), and the scoring specification reserves a `0.20 × docking_score` term for future work. No docking poses have been run: `results/docking_results.csv` is header-only, `scripts/prepare_receptor.py` and `scripts/dock_ligands.py` are stubs, and the delivered `final_score` uses the 4-component formula (`0.40 × ML_rank_score + 0.25 × similarity_support + 0.25 × developability_score + 0.10 × novelty_bonus`). Docking is the highest-priority extension.
- **MD simulation and free-energy calculations** — out of scope; stated as future work in `protocol.md` §1.1.
- **Extension to α-amylase, lipase, or other glycosidases** — out of scope; "future work, not modeled" per `protocol.md` §1.1.
- **Alternative fingerprints (ECFP6) and tautomer standardization** — out of scope; "Locked ECFP4 only" per `protocol.md` §10.
- **Prospective experimental validation of shortlist candidates** — outside a computational hackathon; not modeled.

Nothing in this section is delivered by the pipeline, the web prototype, or the shortlist. The pipeline stops at ligand-based prioritization.

## 2. Pitch summary

Romarinova AI turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system. It links documented Tunisian rosemary chemistry to machine learning, a bounded evidence record, and a three-page web interface.

**The ML pipeline works on its target chemistry.** On a document-held-out test split — 128 train, 30 validation, 21 test compounds, with zero shared `document_chembl_id` across splits — the selected Random Forest achieves **PR-AUC 0.781 and ROC-AUC 0.760 on 7 actives and 14 inactives**, 2.3× the baseline positive-class fraction. A 95% confidence interval is not reported because the active count falls below the protocol's 20-compound threshold (`ci_not_reported = True`, `ef_unstable = True`).

**The pipeline reports where that success ends.** All 28 screened rosemary compounds fall outside the applicability domain (extrapolation at threshold 0.7330; distances 0.801–0.927). On this extrapolation regime, **27 of 28 compounds (96.4%) receive the identical `ML_rank_score` of 0.633333** — the model produces essentially one prediction, and cannot discriminate among small monoterpenes that fall into the same leaf cluster. The ML-independent comparison view makes this visible: pair-level overlap between the locked and ML-independent rankings is **10/10** (`robust`), which reads honestly as "ML was not load-bearing on this regime."

**The shortlist is defensible because the pipeline degrades gracefully.** The locked shortlist (N = 10 pairs, U = 5 compounds — p-cymene, verbenone, γ-terpinene, α-terpinolene, α-pinene) is carried by similarity support and developability, with all five compounds passing every developability check. Characteristic diterpenes (carnosol, carnosic acid, rosmarinic acid) are preserved in diagnostics with their failed-rule annotations. `known_active_count = 0` — no shortlisted compound matches a ChEMBL active InChIKey, so novelty is earned, not borrowed.

**The contribution** is variety- and chemotype-level traceability; scaffold-aware validation with honest small-sample reporting; developability filtering that separates real candidates from large PAINS-prone diterpenes; a side-by-side ML-independent comparison view that makes the model's contribution boundary visible to the user; a specified-but-unexecuted docking axis; and a reproducible open pipeline with locked rules, decision log, and compliance audit.

**The path forward is data, not model complexity.** A Random Forest trained on 14 actives cannot learn to distinguish small monoterpenes it has never seen. Deeper trees, richer fingerprints, or an ensemble would not break the tie-band — they would just memorize the training distribution more efficiently. The two levers that would move the needle are (i) more independently sourced activity data on ChEMBL-like scaffolds relevant to MGAM, and (ii) docking as a second, training-data-independent axis of evidence. Both are named in the report's future-work section.

Closing line: *Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that delivers on-target ML performance and reports where that performance ends.*