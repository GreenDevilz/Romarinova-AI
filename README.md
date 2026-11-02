# Romarinova AI — EthnoHack 2026 Track 2

**From Tunisian Rosemary Chemotypes to Evidence-Guided Candidate Prioritization for α-Glucosidase Inhibition: A Consensus Virtual Screening Protocol**

Track 2 — Natural Products and Drug Discovery · EthnoHack 2026 · Nov 2–8, 2026 · 100% Online

## What this repository is

Romarinova AI is a prioritization and traceability system for documented Tunisian rosemary compounds against human α-glucosidase. It does not claim to treat, prevent, or diagnose diabetes mellitus or any metabolic disease. The pipeline is ligand-based only; no docking is performed.

The contribution is variety/chemotype-level traceability, scaffold-aware validation, developability filtering, a side-by-side ML-independent comparison view, and a reproducible open pipeline.

## Pre-event scaffolding disclosure

This repository is committed before 2 Nov 2026 with **empty notebooks, protocol documentation, repository structures, outline-only presentation templates, dependency files, and placeholder/schema-only scaffolding only**. No project-specific core logic is present before 2 Nov. All data curation, model training, feature generation, AD calibration, scoring, and artifact generation occur during the hackathon.

Allowed pre-event artifacts include: empty notebooks; protocol documents; repository/folder structures; outline-only presentation templates with `[FILL: ...]` slots; `requirements.txt`; empty/placeholder schemas and templates; human-readable locked tool/rule specifications; and the one-page `operating_card.md`.

`models/rules.json`, `models/model_config.json`, `splits/split_indices.json`, `splits/cv_scaffold_groups.json`, `data/*.csv`, `results/*.csv`, and `provenance/provenance_table.csv` are committed as **schema-only or header-only placeholders**. No real splits, model parameters, curated data, or result values are present.

## Open-source starter kits and public libraries (disclosed and attributed)

| Library / resource | Version (to be pinned in-event) | Licence | Use |
|---|---|---|---|
| RDKit | TBD | BSD 3-Clause | Cheminformatics, standardization, ECFP4, PAINS/ADMET/aggregator rules |
| scikit-learn | TBD | BSD 3-Clause | Random Forest, CV, metrics, bootstrap |
| XGBoost | TBD | Apache 2.0 | Gradient-boosted comparator model |
| FastAPI | TBD | MIT | API layer |
| React + Vite | TBD | MIT | Frontend |
| Plotly | TBD | MIT | Figures |
| SQLite | TBD | Public domain | Results store |
| ChEMBL | TBD | CC BY-SA 4.0 | Human α-glucosidase bioactivity source |
| PDB | TBD | Public domain | Not used for docking; cited only for general structural context if needed |
| PubChem | TBD | Public domain | Rosemary compound structure source (priority 1) |

All starter kits and public libraries are disclosed here and audited in `compliance.md`.

## Roster and event confirmations

- Student-only event.
- Submission deadline: 8 Nov 2026.
- Devpost team-size cap confirmed and recorded.
- Three-student roster within the confirmed cap.

## Repository structure
Romarinova-AI/
├── README.md
├── LICENSE
├── .gitignore
├── requirements.txt
├── protocol.md # methodology: final protocol + definitions + AD + scoring + evaluation
├── data.md # data curation + cultivar catalogue
├── compliance.md # rule compliance + licence + responsible research + decision log
├── pitch.md # impact and future work + pitch
├── operating_card.md # one-page daily execution artifact (canonical per protocol)
├── presentation/
│ ├── presentation.md # deck outline + narration outline + narrative text
│ └── devpost.md # Devpost template
├── data/
│ ├── curated_chembl.csv
│ ├── cultivar_compound_catalogue.csv
│ └── features/
├── scripts/
│ ├── chembl_query.sql
│ ├── build_artifacts.py
│ ├── validate_data_contract.py
│ ├── evaluation.py
│ └── sensitivity.py
├── splits/
│ ├── split_indices.json
│ └── cv_scaffold_groups.json
├── models/
│ ├── rule_specs.md
│ ├── rules.json
│ └── model_config.json
├── provenance/
│ └── provenance_table.csv
├── notebooks/
│ ├── 01_data_curation.ipynb
│ ├── 02_stage1_training.ipynb
│ ├── 03_stage1_screening_ad.ipynb
│ ├── 04_stage2_similarity.ipynb
│ ├── 05_stage3_developability.ipynb
│ └── 06_sensitivity_analysis.ipynb
├── results/
│ ├── evaluation_table.csv
│ ├── screened_compounds_annotated.csv
│ ├── primary_shortlist.csv
│ ├── primary_shortlist_compound_level.csv
│ └── diagnostics.csv
├── web/
│ ├── data_contract.json
│ ├── ui_contract.json
│ ├── api/
│ └── frontend/
├── figures/
│   └── .gitkeep              # pipeline_diagram.png and sensitivity_overlap.png produced in-event
└── demo/
    └── .gitkeep              # final_demo_screencast.mp4 recorded in-event

> `.gitkeep` files mark directories whose artifacts are produced in-event. Zero-byte media placeholders are intentionally not committed.

`scripts/build_artifacts.py` is the single source of truth generator in-event. `results/primary_shortlist_compound_level.csv` defines U and is the primary traceability artifact.

## Reproduction (to be completed in-event)

Install `requirements.txt`; run data curation with the committed ChEMBL query and species check; run training with locked splits/CV/params/model; run screening/AD for all eligible rosemary compounds; run similarity/developability; compute `final_score` and `ligand_only_score`; run sensitivity; run `scripts/build_artifacts.py`; serve the web prototype against the frozen UI contract and validate exports.

## Licence

Project code: MIT (see `LICENSE`). Data and third-party libraries retain their own licences as listed above. ChEMBL data is CC BY-SA 4.0; RDKit is BSD 3-Clause; PDB is public domain.

## Disclaimer

Computational predictions are exploratory hypotheses, not medical facts. No compound should be used without experimental validation and clinical trials. Experimental and computational outputs are kept separate.