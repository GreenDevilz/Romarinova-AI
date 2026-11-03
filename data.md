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

**Docking receptor:** PDB (see §3). Not a ChEMBL source; does not affect bioactivity curation.

**Split:** Bemis–Murcko scaffold split first, 70/15/15, stratified by activity (assignment and repair rules in §1.1.1). Save `splits/split_indices.json`. All baselines and Stage 1 models use the identical split. Cross-split InChIKey leakage check removes duplicates from test/validation and documents.

**Day-1 data feasibility ledger:** by end of Day 1, P1 and P2 populate `results/diagnostics.csv` rows with `section = data_feasibility_ledger` recording: ChEMBL human α-glucosidase active/inactive counts after curation and species check; scaffold groups and largest scaffold group size per split; verified identified and suspected pair counts; defined/ambiguous/missing structure counts by source; structure exclusion count and impact diagnostic summary; exact known-active overlap count; compounds eligible for screening. The ledger is a recording artifact; it does not trigger a scope change. Canonical metric keys are listed in §1.1.2.

### §1.1.1 Executable query specification

This operationalizes the locked rules above. It changes no threshold.

**Access method.** ChEMBL SQL dump loaded to a local SQLite instance; no API. Recorded in `compliance.md` §4.

**ChEMBL release.** Release number, download date, and local database path recorded in `compliance.md` §4 and in `splits/split_indices.json`.

**Target resolution (pre-registered).** "Human α-glucosidase" is ambiguous at the UniProt level (GAA/P10253; MGAM/O43451; SI/P14410). The resolution query output is committed as `results/target_resolution.csv`. Rule:

1. Restrict candidates to rows whose `pref_name` denotes an intestinal brush-border α-glucosidase (`Maltase-glucoamylase`, `Sucrase-isomaltase`).
2. If more than one remains, lock the one with the largest qualifying activity count; record the rejected candidates with reasons.
3. If none remains, stop and escalate.

Locked target: `tid = 1` (`Maltase-glucoamylase`, organism *Homo sapiens*), ChEMBL release 37. Rejected candidates and reasons recorded in `compliance.md` §4 row 15. No target change after data freeze without a decision-log row and P1 approval.

The resolution query (executed once):

```sql
-- scripts/target_resolution.sql
-- Output committed to results/target_resolution.csv
SELECT td.tid, td.pref_name, td.organism, td.tax_id,
       COUNT(a.activity_id) AS n_activities
FROM target_dictionary td
LEFT JOIN assays ass     ON ass.tid = td.tid
LEFT JOIN activities a   ON a.assay_id = ass.assay_id
WHERE td.organism = 'Homo sapiens'
  AND td.target_type = 'SINGLE PROTEIN'
  AND (td.pref_name LIKE '%glucosidase%' OR td.pref_name LIKE '%maltase%')
GROUP BY td.tid, td.pref_name, td.organism, td.tax_id
ORDER BY n_activities DESC;

Record classification (per record, before consolidation). Applies to standard_type IN ('IC50','Ki','Kd','EC50') only.

Condition	Class
standard_relation IN ('=','<=','<') and value ≤ 1000 nM	active
standard_relation IN ('=','>','>=') and value ≥ 10000 nM	inactive
value in (1000, 10000) nM, any relation	ambiguous (excluded)
other relation, or unit not convertible	excluded (unparseable)
Unit conversion to nM. nM×1, uM/µM×1e3, mM×1e6, M×1e9, pM×1e-3. Any other unit (%, ug.mL-1, ratios) → excluded as unparseable.

Species check. assays.tax_id = 9606 OR assays.assay_organism = 'Homo sapiens'. Non-human records excluded and counted, even inside a human target.

Duplicate consolidation. Group by standard_inchikey. Representative value = median of standard_value_nM across the compounds qualifying records. Representative metadata (assay_chembl_id, document_chembl_id, confidence_score) taken from the record with the highest confidence_score, ties broken by lowest activity_id. duplicate_flag = True when max/min > 10 (i.e., >1 log unit span); n_records recorded.

Activity conflict rule. If a compound has at least one active and at least one inactive record after classification, curation_flag = conflicting_activity and the compound is excluded from model training and evaluation; counted as activity_conflict_excluded_n.

Seeds. seed = 42 for scaffold split, CV fold assignment, and bootstrap resampling. Recorded in splits/split_indices.json.

Scaffold definition. Bemis–Murcko, generic (includeChirality = False). Scaffold group key = the generic scaffold SMILES; acyclic compounds share the empty-string group.

Split assignment. Sort scaffold groups by size descending, ties by scaffold string ascending. Assign groups greedily to train until ≥70% of compounds are assigned, then to validation until ≥85%, remainder to test. Stratification repair: if validation or test contains zero actives, move the smallest active-containing group of size ≤ 5 from train to the deficient split (ties: lexicographically first scaffold string); repeat up to 3 times; every repair recorded in the ledger (split_repairs_n).

Per-compound structure-retrieval cap. 3 candidate records per source (PubChem, ChEMBL). Beyond the cap, unexamined records are not retrieved; structure_status is assigned from the examined set only.

§1.1.2 Canonical ledger metric keys
results/diagnostics.csv rows with section = data_feasibility_ledger use exactly these metric values:

chembl_actives_n, chembl_inactives_n, chembl_ambiguous_excluded_n, chembl_unparseable_excluded_n, chembl_nonhuman_excluded_n, chembl_confidence_excluded_n, chembl_assaytype_excluded_n, chembl_standardization_failed_n, chembl_activity_conflict_excluded_n, chembl_duplicate_flagged_n, scaffold_groups_train, scaffold_groups_val, scaffold_groups_test, largest_scaffold_group_train, largest_scaffold_group_val, largest_scaffold_group_test, split_repairs_n, cv_folds_n, inchikey_leakage_removed_n, pair_identified_n, pair_suspected_n, structure_defined_n, structure_ambiguous_n, structure_missing_n, structure_source_pubchem_n, structure_source_chembl_n, structure_exclusion_n, known_active_overlap_n, eligible_rosemary_compound_n.

N_star and U_star stay blank on ledger rows (shortlist does not exist until Day 3).

2. Cultivar catalogue specification
2.1 Scope
Documented Tunisian rosemary varieties and chemotypes: var. typicus, var. troglodytorum, and the predominant 1,8-cineole/camphor/α-pinene chemotype.

Source: published LC-MS/GC-MS profiling of Tunisian rosemary varieties and chemotypes. Every compound-variety link individually verified against the primary publication. cultivar_evidence_level ∈ {identified, suspected}; only identified pairs are shortlist-eligible; suspected reported separately (in results/diagnostics.csv with section = suspected_variety_links) and not scored or ranked.

2.2 Eligibility
Only compounds with a usable canonical SMILES from at least one locked structure source — PubChem or ChEMBL — and structure_status = defined are screened. structure_status = ambiguous and structure_status = missing are excluded. Source priority is PubChem → ChEMBL; the chosen source and its provenance are recorded for every screened compound. Compounds with no usable structure in either locked source, or with structure_status = missing or ambiguous, are recorded in results/diagnostics.csv with section = missing_structures and excluded. No manual structure retrieval beyond the per-compound cap (§1.1.1). structure_status is assigned by locked rules from the chosen structure.

structure_status	Definition	Screening	Reporting
defined	Single unambiguous structure with defined stereochemistry.	Included.	InChIKey and source.
ambiguous	Multiple forms, undefined stereochemistry, or conflicting sources.	Excluded; listed in results/diagnostics.csv with section = missing_structures, reason ambiguous_only.	Diagnostics data and Devpost coverage gap.
missing	No usable structure in any locked source within the per-compound cap, or standardization failure that removes the structure.	Excluded; listed in results/diagnostics.csv with section = missing_structures.	Diagnostics data and Devpost coverage gap.
Source priority: PubChem → ChEMBL. When multiple sources provide a structure, the highest-priority source is used; the chosen source is recorded. Conflicts use PubChem and set ambiguous when the discrepancy is material. If PubChem has no usable defined structure, ChEMBL is used; if neither yields a defined structure, the compound is missing or ambiguous and excluded.

Structure exclusion gap: compounds with no defined structure in either PubChem or ChEMBL are excluded from screening. This is a locked eligibility rule and is reported as an explicit coverage-gap source: counted in results/diagnostics.csv (sections data_feasibility_ledger and missing_structures), named in the coverage-gap narrative, and disclosed in the Devpost page and video.

Shortlist universe: all eligible rosemary compounds with ≥1 identified pair and structure_status = defined. Similarity, developability, and final_score are computed for all eligible rosemary compounds; AD is also computed for all eligible rosemary compounds before the ranking lock. Suspected-only and ambiguous/missing-structure compounds are outside the scoring and ranked universe.

2.3 Evidence fields (never collapsed)
Field	Values	Meaning
cultivar_evidence_level	identified, suspected	Compound presence in the rosemary variety/chemotype. Only identified is shortlist-eligible.
literature_evidence_level	experimental, none	Prior known-active support, limited to exact ChEMBL known-active flags. experimental only when known_active_flag = True; otherwise none. This is not a literature review.
graph_edge_evidence_level	Not used.	No graph.
known_active_flag definition: known_active_flag = True if and only if the compounds standardized InChIKey matches any curated ChEMBL active after the species check and curation rules in §1.1. exact_training_match and exact_test_match are subset flags. literature_evidence_level = experimental whenever known_active_flag = True, regardless of split.

2.4 Reporting
Cultivar presence is provenance-based. The project does not rank varieties against each other. all_cultivars is authoritative; primary_cultivar is display-only.


# 3. Docking input preparation

## 3.1 Receptor

- PDB entry, chain, resolution, and missing-residue handling are locked in `models/docking_config.json` and recorded in `provenance/provenance_table.csv`.
- Preparation is mechanical and scripted (`scripts/prepare_receptor.py`): remove crystallographic waters and non-catalytic heteroatoms; retain the catalytic Zn; add hydrogens at pH 7.4; assign Gasteiger charges; write PDBQT with Meeko. No hand-edited receptor.
- Alternative receptor structures reviewed at lock time are recorded as rejected candidates with reasons in `compliance.md` §4.

## 3.2 Ligands

- 3D conformer from the locked canonical SMILES: RDKit ETKDGv3 (`randomSeed = 42`), MMFF94 optimization, protonation at pH 7.4; single lowest-energy conformer per compound.
- PDBQT conversion with Meeko; no manual ligand editing.
- Compounds with `structure_status ∈ {ambiguous, missing}` are not docked; they are recorded in `results/diagnostics.csv` with section = `missing_structures`.

## 3.3 Grid box

- Single grid box centered on the catalytic site; center and size locked in `models/docking_config.json`. One box for all compounds; no per-compound redocking box.

## 3.4 Docking outputs

- Per-compound: `docking_affinity_kcal_per_mol`, `docking_score`, `docking_missing`, `docking_pose_status`.
- Written to `results/docking_results.csv` by `scripts/dock_ligands.py`; consumed by `scripts/build_artifacts.py`.
- Docking failures are retained with `docking_missing = True` and are not imputed.