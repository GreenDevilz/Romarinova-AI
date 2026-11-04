# Pitch and Impact

## 1. Impact and future work

### 1.1 Prioritized developable compound-variety candidates

Romarinova AI delivers prioritized developable compound-variety candidates with provenance and evidence from three orthogonal axes: ligand-based ML screening, ligand-based similarity/developability scoring, and receptor-based docking. Impact is quantified by the actual eligible rosemary compound N (28), shortlist N (10 pairs), and U (5 unique compounds).

### 1.2 Modest-shortlist framing

The locked shortlist (N = 10 pairs, U = 5 unique compounds) is composed of small monoterpenes and a phenolic monoterpenoid: **p-cymene** (rank 1), **γ-terpinene** (rank 2), **α-phellandrene** (rank 3), **linalool** (rank 4), and **myrcene** (rank 5). 

Every shortlisted compound passes all four developability checks (Lipinski, Veber, PAINS, aggregator) and shows favorable docking affinities (-3.0 to -4.3 kcal/mol) in the MGAM catalytic site. The most rosemary-characteristic diterpenes — carnosic acid, carnosol, and rosmarinic acid — appear only in the diagnostics data. They fail the PAINS filter and exhibit steric clashes in the rigid receptor site (positive docking penalties). 

This is prioritization discipline: the pipeline separates developable, receptor-compatible candidates from large, PAINS-prone, low-bioavailability diterpenes, preserving the latter in diagnostics for experimental follow-up.

### 1.3 ML delivers on-domain; the pipeline reports where its domain ends

**The model works on target chemistry.** On a document-held-out test split (128 train, 30 validation, 21 test compounds), the selected Random Forest achieves **PR-AUC 0.781 and ROC-AUC 0.760** on 7 actives and 14 inactives. A similarity-only baseline reaches 0.735, meaning the ML model provides a modest, honest uplift.

**The pipeline's value is what it does at the domain boundary.** Screening the rosemary catalogue, all 28 compounds fall outside the applicability domain (ad_distance 0.8014–0.9265 vs. threshold 0.7330). On this extrapolation regime, **27 of 28 screened compounds (96.4%) receive the identical `ML_rank_score` of 0.633333**. The ML model produces essentially one prediction for the entire extrapolation regime.

**The side-by-side comparison view is the instrument that surfaces this.** The α-Glucosidase Screening page shows the locked ranking and the ML-independent ranking simultaneously. Pair-level overlap is **10/10** — the sensitivity label is **`robust`**. 

The correct reading of `robust` here is: **ML was not load-bearing on this regime, so removing it changes nothing**. All five shortlisted compounds share `ML_rank_score = 0.633333`, `developability_score = 1.0`, and `novelty_bonus = 1.0`; the rank order is determined entirely by `similarity_support` and `docking_score`. A pipeline that hid this would appear more impressive and be less trustworthy.

### 1.4 Compound-level view

The compound-level view is the primary traceability artifact, with U = 5 reported. Every row carries `ligand_only_score`/`ligand_only_rank`, `docking_affinity_kcal_per_mol`, `docking_score`, `primary_cultivar` (display-only) and `all_cultivars` (authoritative) both shown with a tooltip, and the exact-match flags. `known_active_count = 0` — novelty is earned, not borrowed.

### 1.5 Out of scope — future work only

The following are explicitly **not** delivered by this submission:

- **MD simulation and free-energy calculations** — out of scope; stated as future work in `protocol.md` §1.1. Docking is delivered as a single-receptor, single-conformer scoring axis.
- **Extension to α-amylase, lipase, or other glycosidases** — out of scope.
- **Alternative fingerprints (ECFP6) and tautomer standardization** — out of scope.
- **Prospective experimental validation of shortlist candidates** — outside a computational hackathon; not modeled.

## 2. Pitch summary

Romarinova AI turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system. It links documented Tunisian rosemary chemistry to machine learning, receptor-based docking, a bounded evidence record, and a three-page web interface.

**The ML pipeline works on its target chemistry.** The selected Random Forest achieves **PR-AUC 0.781 and ROC-AUC 0.760** on a leakage-free test split. A 95% confidence interval is not reported because the active count falls below the protocol's 20-compound threshold (`ci_not_reported = True`, `ef_unstable = True`).

**The pipeline reports where that success ends.** All 28 screened rosemary compounds fall outside the applicability domain. On this extrapolation regime, **27 of 28 compounds (96.4%) receive the identical `ML_rank_score` of 0.633333**. The ML-independent comparison view makes this visible: pair-level overlap is **10/10** (`robust`), which reads honestly as "ML was not load-bearing on this regime."

**The shortlist is defensible because the pipeline degrades gracefully.** The locked shortlist (N = 10 pairs, U = 5 compounds — p-cymene, γ-terpinene, α-phellandrene, linalool, myrcene) is carried by similarity support, developability, and receptor docking. All five compounds pass every developability check and dock favorably into the MGAM active site.

**The contribution** is variety- and chemotype-level traceability; scaffold-aware validation with honest small-sample reporting; developability filtering; receptor-based docking; a side-by-side ML-independent comparison view that makes the model's contribution boundary visible; and a reproducible open pipeline with locked rules, decision log, and compliance audit.

Closing line: *Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that delivers on-target ML performance, receptor-based docking, and reports where each component's contribution ends.*