# Pitch and Impact

## 1. Impact and future work

### 1.1 Prioritized developable compound-variety candidates

Romarinova AI delivers prioritized developable compound-variety candidates with provenance and evidence from three axes: ligand-based ML screening, ligand-based similarity/developability scoring, and a receptor-based docking specification that is integrated into the scoring formula but not executed in this submission. Impact is quantified by the actual eligible rosemary compound N (28), shortlist N (10 pairs), and U (5 unique compounds).

### 1.2 Modest-shortlist framing

The locked shortlist (N = 10 pairs, U = 5 unique compounds) is composed entirely of small monoterpenes and sesquiterpenes: borneol (rank 1), α-eudesmol (rank 2), camphene (rank 3), β-pinene (rank 5), and camphor (rank 6). The most rosemary-characteristic diterpenes — carnosic acid, carnosol, and rosmarinic acid — appear only in the diagnostics data. Carnosol and carnosic acid both fail the PAINS developability filter; rosmarinic acid fails both Veber and PAINS. The developability filters are doing measurable work on the diterpene class.

This is prioritization discipline, not a failure: the pipeline correctly separates developable candidates from large, PAINS-prone, low-bioavailability diterpenes, and the diagnostics data preserves them for experimental follow-up.

### 1.3 Side-by-side comparison view

The α-Glucosidase Screening page shows the locked ranking and the ML-independent ranking side by side and simultaneously. This makes the domain-shift response concrete and lets a viewer see how much of the shortlist is ML-driven versus similarity- and developability-driven.

The result is the clearest finding of the corrected run. On the leakage-free test split the model achieves PR-AUC 0.781 on 7 actives and 14 inactives — a strong signal on ChEMBL-like chemistry. But **all five shortlisted rosemary compounds receive the identical `ML_rank_score` (0.777778)**: a Random Forest trained on 14 actives cannot discriminate among monoterpene and sesquiterpene scaffolds that all fall into the same "no tree path fires" leaf cluster. ML therefore contributes **nothing** to the ordering of the shortlist, even though it discriminates well on the test set. The comparison view is the instrument that shows this collapse — the domain-shift caveat rendered as a UI element, not a footnote.

Pair-level overlap between the locked ranking and the ML-independent ranking is 6/10 (`partially_robust`); compound-level overlap is 4/5. The two rankings agree because ML was already silent.

### 1.4 Compound-level view

The compound-level view is the primary traceability artifact, with U = 5 reported. Every row carries `ligand_only_score`/`ligand_only_rank`, `primary_cultivar` (display-only) and `all_cultivars` (authoritative) both shown with a tooltip. On the corrected shortlist, `high_ml_ad_extrapolation_count = 0` — but this reflects the flag's ≥ 0.80 threshold, not the absence of risk: all 5 shortlist compounds have `ad_flag = extrapolation` and sit at `ML_rank_score = 0.777778`, just below the flag boundary. `ad_extrapolation_count = 5` reports the true picture.

### 1.5 Out of scope — future work only

The following are explicitly **not** delivered by this submission:

- **Docking execution.** The docking protocol is fully specified (`models/docking_config.json` locks receptor PDB 3L4W, grid box, Vina parameters, and the `docking_score` transform), and the scoring formula preserves a `0.20 × docking_score` term. No docking poses have been run in this submission: `results/docking_results.csv` is header-only, `scripts/prepare_receptor.py` and `scripts/dock_ligands.py` are stubs, and the scoring actually delivered uses the 4-component formula (`0.40 × ML_rank_score + 0.25 × similarity_support + 0.25 × developability_score + 0.10 × novelty_bonus`). Docking is the highest-priority extension.
- **MD simulation and free-energy calculations** — out of scope; stated as future work in `protocol.md` §1.1.
- **Extension to α-amylase, lipase, or other glycosidases** — out of scope; "future work, not modeled" per `protocol.md` §1.1.
- **Alternative fingerprints (ECFP6) and tautomer standardization** — out of scope; "Locked ECFP4 only" per `protocol.md` §10.
- **Prospective experimental validation of shortlist candidates** — outside a computational hackathon; not modeled.

Nothing in this section is delivered by the pipeline, the web prototype, or the shortlist. The pipeline stops at ligand-based prioritization.

## 2. Pitch summary

Romarinova AI turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system. It links documented Tunisian rosemary chemistry to machine learning, a bounded evidence record, and a three-page web interface.

On the corrected, document-aware split (train 128, val 30, test 21 compounds; zero shared documents across any pair of splits), test-set performance is **PR-AUC 0.781 and ROC-AUC 0.760 on 7 actives and 14 inactives**; a 95% confidence interval is not reported because the active count falls below the protocol's 20-compound threshold (`ci_not_reported = True`, `ef_unstable = True`).

The locked shortlist (N = 10 pairs, U = 5 compounds) is composed of small monoterpenes and sesquiterpenes; characteristic diterpenes are preserved in diagnostics. The ML-independent comparison view shows pair-level overlap of 6/10 (`partially_robust`) and compound-level overlap of 4/5 — meaning the shortlist ordering is largely preserved under ML-weight removal, because all five shortlisted compounds receive the same ML rank and ML cannot reorder them. All 28 rosemary compounds are outside the applicability domain (extrapolation at threshold 0.7330; distances 0.801–0.927), and the pipeline reports this as a headline caveat rather than suppressing it.

The contribution is variety- and chemotype-level traceability, scaffold-aware validation with honest uncertainty, developability filtering, a side-by-side ML-independent comparison view that makes the ML collapse visible, a specified-but-unexecuted docking axis, and a reproducible open pipeline. It does not claim rosemary cures diabetes; it identifies which documented rosemary compounds merit further experimental investigation and shows the evidence — and the limits — behind each ranking.

Closing line: *Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that reports where its own ML contribution ends.*