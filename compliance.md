# Compliance

Rule compliance, licence audit, responsible research, and the decision log.

---

# 1. Rule compliance and in-event scope

**Compliance principle:** no data curation, model training, feature generation, AD calibration, or decoy generation before 2 Nov. No real ChEMBL or rosemary-profiling data is queried before 2 Nov.

## 1.1 Allowed before 2 Nov

- Empty notebooks.
- Protocol documentation.
- Repository and folder structures.
- Outline-only presentation templates with section headings and `[FILL: ...]` slots only.
- Dependency files.
- Empty/placeholder schemas and templates.
- Human-readable locked tool/rule specifications as protocol documentation.
- Reading titles, abstracts, and bibliographic metadata.
- The one-page `operating_card.md`.

## 1.2 Not allowed before 2 Nov

- Project-specific core logic.
- Working data curation/standardization code.
- Working ChEMBL query execution.
- Feature generation, split/CV, model training/tuning/evaluation, AD, similarity, developability, scoring, shortlist, sensitivity, or data-contract validation logic.
- Synthetic-data execution of project logic.
- Realistic mock-data rehearsal in the frontend.

## 1.3 Pre-event checks

Limited to dependency/import checks, empty-template rendering, and protocol/tabletop walkthroughs. No synthetic-data execution of project-specific curation, ML, AD, scoring, or frontend data logic before 2 Nov.

## 1.4 Transparency

All open-source starter kits and public libraries are disclosed in `README.md` with names, versions, licences, and attribution, and audited in §2.

## 1.5 Pre-event rule specifications

Human-readable locked rule specifications for ADMET, aggregator, and PAINS are committed pre-event as `models/rule_specs.md`, with rule IDs, thresholds, source citations, and prose descriptions of each pattern class in Markdown tables. `rule_specs.md` contains no literal SMARTS strings; exact SMARTS are committed in-event at implementation time with a decision-log row. The lock is preserved because the source, selection rule, and thresholds are locked pre-event; only the literal encoding is deferred.

`models/rules.json` is committed pre-event as a schema-only file and is not imported by any script before 2 Nov.

## 1.6 Organizer confirmations

- Student-only event.
- Submission deadline on 8 Nov 2026.
- Devpost team-size cap confirmed and recorded.
- Three-student roster within cap.

## 1.7 Pre-event completion gate

All committed pre-event artifacts are empty templates, protocol documents, repository structures, outline-only presentation templates, or dependencies. No project-specific core logic is present. Rule specification files contain only human-readable protocol specs; no script imports `models/rules.json` or `models/model_config.json` before 2 Nov. A SMARTS-detector grep confirms `models/rule_specs.md` contains no literal pattern strings. Three-part mechanical presentation-text audit recorded.

---

# 2. Licence check

| Resource | Licence | Use | Attribution location |
|---|---|---|---|
| ChEMBL | CC BY-SA 4.0 | Bioactivity data | `README.md`, Devpost |
| RDKit | BSD 3-Clause | Cheminformatics | `README.md` |
| scikit-learn | BSD 3-Clause | ML | `README.md` |
| XGBoost | Apache 2.0 | ML | `README.md` |
| FastAPI | MIT | API | `README.md` |
| React + Vite | MIT | Frontend | `README.md` |
| Plotly | MIT | Figures | `README.md` |
| SQLite | Public domain | Results store | `README.md` |
| PDB | Public domain | Not used for docking | `README.md` |
| PubChem | Public domain | Structure source | `README.md` |

**Pre-event status:** all licences disclosed and attributed; no project-specific logic present.

---

# 3. Responsible research

- Traditional-knowledge attribution to source databases and, where available, originating communities; no appropriation or decontextualization.
- Benefit-sharing and prior informed consent are beyond hackathon scope; this is a computational prioritization tool, not a claim on traditional knowledge.
- PAINS/ADMET are heuristic flags, not truths; natural products often trigger PAINS; the diagnostics data preserves follow-up candidates.
- Computational predictions are not medical facts; no compound should be used without experimental validation and clinical trials; experimental and computational outputs kept separate.
- Known-compound critique: the ranking may naturally mix known-scaffold analogs and novel/PAINS-flagged candidates; the pipeline does not force novelty. Contribution is prioritization and transparency, not de novo discovery.
- No clinical claims: no compound treats, prevents, or diagnoses diabetes or any metabolic disease; all outputs are exploratory hypotheses.

---

# 4. Decision log

**Cap:** ≤20 rows with pre-assigned row budgets. P1 is single technical decision-maker; P2 is deputy. Decision-log entry promptly; no committee votes.

| # | Date | Decision / deviation | Rationale | Owner | Status |
|---|---|---|---|---|---|
| 1 | Pre-event | `team.md` omitted; roster confirmed in `README.md`; impact/pitch content consolidated in `pitch.md` | Simplification; no content dropped | P1 | Locked |
| 2 | Pre-event | `rule_specs.md` contains no literal SMARTS; exact SMARTS committed in-event | Preserve pre-event lock on source/thresholds while deferring literal encoding | P1 | Locked |
| 3 | Pre-event | `models/rules.json`, `models/model_config.json`, `splits/*.json`, `data/*.csv`, `results/*.csv`, `provenance/provenance_table.csv` committed as schema/header-only placeholders | Comply with no pre-event project-specific data/logic | P1 | Locked |
| 4 | Pre-event | `scripts/build_artifacts.py`, `validate_data_contract.py`, `chembl_query.sql`, `evaluation.py`, `sensitivity.py` committed as empty/comment-only stubs | Comply with no pre-event executable logic | P1 | Locked |
| 5 | Pre-event | Presentation templates are outline-only with `[FILL: ...]` slots | Comply with no pre-event frontend data logic | P3 | Locked |
| 6 | Pre-event | `results/results_store.sqlite` created in-event, not committed pre-event | Avoid populated results store pre-event | P1 | Locked |
| 7 | Pre-event | Root markdown consolidated: `final_protocol` + `definitions_and_consistency` + `ad` + `scoring_rule` + `evaluation` → `protocol.md`; `data_curation` + `cultivar_catalogue` → `data.md`; `rule_compliance` + `licence_check` + `responsible_research` + `decision_log` → `compliance.md`; `impact_and_future` + `pitch` → `pitch.md`; presentation 5→2; models 8→3; results 8→5 | Reduce file count (~45% reduction) without dropping reporting content | P1 | Locked |
| 8 | In-event | [Reserved: exact SMARTS encoding for PAINS/aggregator] | Implementation-time encoding | P2 | Pending |
| 9 | In-event | ChEMBL release 37; SQLite dump downloaded 2026-11-03; local database path recorded in `splits/split_indices.json`; target-resolution output committed as `results/target_resolution.csv` | Data freeze record | P1 | Locked |
| 10 | In-event | Model freeze record: RandomForestClassifier selected (CV PR-AUC 0.264 ± 0.244 vs XGB 0.261 ± 0.221). Test metrics: PR-AUC 0.370, Top-5 recovery 2/5. Small-N CV artifact: some folds had 0 actives, triggering warnings; `ef_unstable=True`, `ci_not_reported=True`. Applicability Domain (AD): LOO threshold = 0.4638. All 28 rosemary compounds flagged as `extrapolation` (distances > 0.88). | Model freeze record. AD extrapolation confirms protocol Caveat R2 (Domain Shift); high ML scores on these compounds will be flagged `high_ml_ad_extrapolation=True` per protocol. | P1 | Locked |
| 11 | In-event | Ranking lock verified. Primary shortlist N* = 10 (capped from 28 pairs); Compound-level U* = 28. Metadata: known_active_count = 0, high_ml_ad_extrapolation_count = 10, ad_extrapolation_count = 10, ad_missing_count = 0. screened_compounds_annotated.csv locked. | Day 3 ranking lock record. The 100% AD extrapolation rate confirms the domain-shift expectation for natural products. | P2 | Locked |
| 12 | In-event | Final record verification. Shortlist (N=10) and compound-level view (U=28) verified. Sensitivity analysis yielded pair-level overlap 6/10 (partially_robust) and compound-level overlap 4/5. Diagnostics CSV generated and deduplicated (filtered-but-interesting, disagreement patterns). | Day 4 verification record; validates ML marginality (Caveat R3). | P2 | Locked |
| 13 | In-event | `data.md` §1.1.1 appended: executable query specification, target-resolution rule, seeds, scaffold definition, split assignment, per-compound structure-retrieval cap; no threshold changed | Locked bullets were not executable; `chembl_query.sql` must be auditable against a concrete specification | P1 | Locked |
| 14 | In-event | Curated CSV, cultivar catalogue, provenance table, and diagnostics headers expanded; diagnostics `section` enum closed; evidence fields kept separate (`cultivar_evidence_level`, `literature_evidence_level`) | Report required fields absent from pre-event placeholders; joins across `compound_id` would have failed | P1 | Locked |
| 15 | In-event Day 1 | ChEMBL target locked: `tid = 1` (`Maltase-glucoamylase`, organism *Homo sapiens*), ChEMBL release 37. Target-resolution query executed; output committed as `results/target_resolution.csv`. Rejected candidates recorded with reasons: `tid 10977` Lysosomal alpha-glucosidase — highest activity count but wrong biology (lysosomal glycogen degradation, Pompe disease; not the intestinal brush-border enzyme relevant to postprandial hyperglycemia); `tid 11815`, `tid 11816` — β-glucosidases, wrong enzyme class; `tid 10734`, `tid 10736`, `tid 100190`, `tid 100189` — glycoprotein-processing α-glucosidases, not intestinal; `tid 11997` Sucrase-isomaltase — correct biology but excluded by the single-target scope rule. `data.md` §1.1.1 resolution rule amended to restrict the activity-count tie-break to intestinal brush-border α-glucosidase candidates. | Target must match the project's stated mechanism; literal application of the count tie-break would have locked the lysosomal enzyme | P1 | Locked |
| 16 | In-event | Duplicate consolidation: group by `standard_inchikey`, median value, metadata from highest-confidence record (ties: lowest `activity_id`); `duplicate_flag` at >1 log unit span; activity-conflict exclusion rule added (`curation_flag = conflicting_activity`) | "Median with highest-confidence assay" was ambiguous; the conflict rule is required because the curated set must be model-trainable | P1 | Locked |
| 17 | In-event | Scaffold split: greedy size-descending 70/15/15 with active-stratification repair (≤3 repairs, each recorded); acyclic compounds share the empty-string scaffold group | Stratified scaffold split is approximable only; the repair rule makes it deterministic and auditable | P1 | Locked |
| 18 | In-event | Per-compound structure-retrieval cap = 3 candidate records per source (PubChem, ChEMBL) | Referenced in report but unvalued; cap makes `structure_status` assignment deterministic | P1 | Locked |
| 19 | In-event | `requirements.txt` pinned from resolved environment; RDKit version recorded explicitly | Curation reproducibility; RDKit version affects standardization/neutralization output | P1 | Locked |
| 20 | In-event | Schedule shift: Data freeze moved to start of Nov 3 (due to environment setup and target resolution). Day 2 tasks (model training and evaluation) will still be executed on Nov 3 as planned. | Environment setup and target resolution required execution time on Nov 3; timeline compressed to preserve model freeze date | P1 | Locked |