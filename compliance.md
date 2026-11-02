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
| 1 | Pre-event | `team.md` and standalone `pitch.md` omitted; roster confirmed in `README.md`; impact/pitch content merged into `pitch.md` | Simplification; no content dropped | P1 | Locked |
| 2 | Pre-event | `rule_specs.md` contains no literal SMARTS; exact SMARTS committed in-event | Preserve pre-event lock on source/thresholds while deferring literal encoding | P1 | Locked |
| 3 | Pre-event | `models/rules.json`, `models/model_config.json`, `splits/*.json`, `data/*.csv`, `results/*.csv`, `provenance/provenance_table.csv` committed as schema/header-only placeholders | Comply with no pre-event project-specific data/logic | P1 | Locked |
| 4 | Pre-event | `scripts/build_artifacts.py`, `validate_data_contract.py`, `chembl_query.sql`, `evaluation.py`, `sensitivity.py` committed as empty/comment-only stubs | Comply with no pre-event executable logic | P1 | Locked |
| 5 | Pre-event | Presentation templates are outline-only with `[FILL: ...]` slots | Comply with no pre-event frontend data logic | P3 | Locked |
| 6 | Pre-event | `results/results_store.sqlite` created in-event, not committed pre-event | Avoid populated results store pre-event | P1 | Locked |
| 7 | Pre-event | Root markdown consolidated: `final_protocol` + `definitions_and_consistency` + `ad` + `scoring_rule` + `evaluation` → `protocol.md`; `data_curation` + `cultivar_catalogue` → `data.md`; `rule_compliance` + `licence_check` + `responsible_research` + `decision_log` → `compliance.md`; `impact_and_future` + `pitch` → `pitch.md`; presentation 5→2; models 8→3; results 8→5 | Reduce file count (~45% reduction) without dropping reporting content | P1 | Locked |
| 8 | In-event | [Reserved: exact SMARTS encoding for PAINS/aggregator] | Implementation-time encoding | P2 | Pending |
| 9 | In-event | [Reserved: ChEMBL download date/version] | Data freeze record | P1 | Pending |
| 10 | In-event | [Reserved: model selection outcome XGB vs RF] | Model freeze record | P1 | Pending |
| 11 | In-event | [Reserved: ranking lock verification] | Day 3 lock | P2 | Pending |
| 12 | In-event | [Reserved: final record verification] | Day 4 verification | P2 | Pending |
| 13 | In-event | [Reserved] | | | Pending |
| 14 | In-event | [Reserved] | | | Pending |
| 15 | In-event | [Reserved] | | | Pending |
| 16 | In-event | [Reserved] | | | Pending |
| 17 | In-event | [Reserved] | | | Pending |
| 18 | In-event | [Reserved] | | | Pending |
| 19 | In-event | [Reserved] | | | Pending |
| 20 | In-event | [Reserved] | | | Pending |
| 9 | Pre-event | Tree aligned to `git ls-files`: root `data_feasibility_ledger.csv` removed (merged into `results/diagnostics.csv` section `data_feasibility_ledger`); `figures/` and `demo/` use `.gitkeep` placeholders; real media artifacts produced in-event | `.gitkeep` is the conventional marker for directories whose content arrives later; zero-byte media placeholders are intentionally not committed | P1 | Locked |
| 10 | Pre-event | `pitch.md` future-work section reframed as explicit out-of-scope with protocol citations | Prevent appearance of scope creep; protocol already names these as future work | P1 | Locked |