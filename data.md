# Data

Curation and provenance specification for ChEMBL human α-glucosidase bioactivity and the rosemary cultivar catalogue.

---

# 1. Data curation specification

## 1.1 ChEMBL bioactivity

**Locked rules:**

- Human α-glucosidase only.
- Active: IC50/Ki/Kd/EC50 ≤ 1 µM.
- Inactive: ≥ 10 µM.
- Ambiguous 1–10 µM excluded.
- Units standardized to nM.
- ChEMBL confidence ≥ 7.
- Assay types: IC50/Ki/Kd/EC50 only; % inhibition and cell-viability-only excluded.
- Duplicates by median with highest-confidence assay; flag if >1 log unit span.
- RDKit standardization: remove salts, neutralize, canonical SMILES, largest fragment.
- Keep defined stereochemistry and flag undefined stereocenters.
- Do not remove PAINS from training.
- Record ChEMBL version, download date, query.
- Save curated CSV to `data/curated_chembl.csv`.

**Species reconciliation:** the query targets the human α-glucosidase ChEMBL target; confidence ≥ 7 applies within that target. Because confidence reflects target-assignment confidence, not species, the curation script performs an explicit species check and excludes any non-human record. Excluded records are counted and documented in `compliance.md` (decision log) and the diagnostics data. The curated CSV contains human records only.

**Pre-registered query specification:** committed here; executable `scripts/chembl_query.sql` is implemented in-event and must match this specification. No ad-hoc variants after data freeze. ChEMBL is the only bioactivity source.

**Split:** Bemis–Murcko scaffold split first, 70/15/15, stratified by activity. Save `splits/split_indices.json`. All baselines and Stage 1 models use the identical split. Cross-split InChIKey leakage check removes duplicates from test/validation and documents.

**Day-1 data feasibility ledger:** by end of Day 1, P1 and P2 populate `results/diagnostics.csv` rows with `section = data_feasibility_ledger` recording: ChEMBL human α-glucosidase active/inactive counts after curation and species check; scaffold groups and largest scaffold group size per split; verified identified and suspected pair counts; defined/ambiguous/missing structure counts by source; structure exclusion count and impact diagnostic summary; exact known-active overlap count; compounds eligible for screening. The ledger is a recording artifact; it does not trigger a scope change.

---

# 2. Cultivar catalogue specification

## 2.1 Scope

Documented Tunisian rosemary varieties and chemotypes: var. *typicus*, var. *troglodytorum*, and the predominant 1,8-cineole/camphor/α-pinene chemotype.

**Source:** published LC-MS/GC-MS profiling of Tunisian rosemary varieties and chemotypes. Every compound-variety link individually verified against the primary publication. `cultivar_evidence_level ∈ {identified, suspected}`; only identified pairs are shortlist-eligible; suspected reported separately (in `results/diagnostics.csv` with `section = suspected_variety_links`) and not scored or ranked.

## 2.2 Eligibility

Only compounds with a usable canonical SMILES from at least one locked structure source — PubChem or ChEMBL — and `structure_status = defined` are screened. `structure_status = ambiguous` and `structure_status = missing` are excluded. Source priority is PubChem → ChEMBL; the chosen source and its provenance are recorded for every screened compound. Compounds with no usable structure in either locked source, or with `structure_status = missing` or `ambiguous`, are recorded in `results/diagnostics.csv` with `section = missing_structures` and excluded. No manual structure retrieval beyond the per-compound cap. `structure_status` is assigned by locked rules from the chosen structure.

| `structure_status` | Definition | Screening | Reporting |
|---|---|---|---|
| **defined** | Single unambiguous structure with defined stereochemistry. | Included. | InChIKey and source. |
| **ambiguous** | Multiple forms, undefined stereochemistry, or conflicting sources. | Excluded; listed in `results/diagnostics.csv` with `section = missing_structures`, reason `ambiguous_only`. | Diagnostics data and Devpost coverage gap. |
| **missing** | No usable structure in any locked source within the per-compound cap, or standardization failure that removes the structure. | Excluded; listed in `results/diagnostics.csv` with `section = missing_structures`. | Diagnostics data and Devpost coverage gap. |

**Source priority:** PubChem → ChEMBL. When multiple sources provide a structure, the highest-priority source is used; the chosen source is recorded. Conflicts use PubChem and set ambiguous when the discrepancy is material. If PubChem has no usable defined structure, ChEMBL is used; if neither yields a defined structure, the compound is missing or ambiguous and excluded.

**Structure exclusion gap:** compounds with no defined structure in either PubChem or ChEMBL are excluded from screening. This is a locked eligibility rule and is reported as an explicit coverage-gap source: counted in `results/diagnostics.csv` (sections `data_feasibility_ledger` and `missing_structures`), named in the coverage-gap narrative, and disclosed in the Devpost page and video.

**Shortlist universe:** all eligible rosemary compounds with ≥1 identified pair and `structure_status = defined`. Similarity, developability, and `final_score` are computed for all eligible rosemary compounds; AD is also computed for all eligible rosemary compounds before the ranking lock. Suspected-only and ambiguous/missing-structure compounds are outside the scoring and ranked universe.

## 2.3 Evidence fields (never collapsed)

| Field | Values | Meaning |
|---|---|---|
| `cultivar_evidence_level` | `identified`, `suspected` | Compound presence in the rosemary variety/chemotype. Only `identified` is shortlist-eligible. |
| `literature_evidence_level` | `experimental`, `none` | Prior known-active support, limited to exact ChEMBL known-active flags. `experimental` only when `known_active_flag = True`; otherwise `none`. This is not a literature review. |
| `graph_edge_evidence_level` | Not used. | No graph. |

**`known_active_flag` definition:** `known_active_flag = True` if and only if the compound's standardized InChIKey matches any curated ChEMBL active after the species check and curation rules in §1.1. `exact_training_match` and `exact_test_match` are subset flags. `literature_evidence_level = experimental` whenever `known_active_flag = True`, regardless of split.

## 2.4 Reporting

Cultivar presence is provenance-based. The project does not rank varieties against each other. `all_cultivars` is authoritative; `primary_cultivar` is display-only.