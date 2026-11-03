-- Locked target: tid = 1 (Maltase-glucoamylase, Homo sapiens), ChEMBL 37.
-- Per data.md §1.1.1. No ad-hoc variants after data freeze.
-- Raw output written to data/raw_chembl_tid1.csv (immutable).

SELECT
    a.activity_id                AS activity_id,
    md.chembl_id                 AS molecule_chembl_id,
    pmd.chembl_id                AS parent_molecule_chembl_id,
    td.chembl_id                 AS target_chembl_id,
    ass.chembl_id                AS assay_chembl_id,
    d.chembl_id                  AS document_chembl_id,
    a.standard_type              AS activity_type,
    a.standard_relation          AS standard_relation,
    a.standard_value             AS standard_value_raw,
    a.standard_units             AS standard_units_raw,
    ass.assay_organism           AS assay_organism,
    ass.assay_tax_id             AS assay_tax_id,
    ass.confidence_score         AS confidence_score,
    cs.canonical_smiles          AS canonical_smiles,
    cs.standard_inchi_key        AS standard_inchikey,
    a.pchembl_value              AS pchembl_value
FROM activities a
JOIN assays ass              ON ass.assay_id = a.assay_id
JOIN target_dictionary td    ON td.tid = ass.tid
JOIN molecule_dictionary md  ON md.molregno = a.molregno
LEFT JOIN molecule_hierarchy mh   ON mh.molregno = a.molregno
LEFT JOIN molecule_dictionary pmd ON pmd.molregno = mh.parent_molregno
LEFT JOIN compound_structures cs  ON cs.molregno = a.molregno
LEFT JOIN docs d             ON d.doc_id = a.doc_id
WHERE td.tid = 1;