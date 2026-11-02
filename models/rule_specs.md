# Human-Readable Locked Rule Specifications

**Pre-event lock:** This file contains rule IDs, thresholds, source citations, and prose descriptions only. It contains **no literal SMARTS strings**. Exact SMARTS are committed in-event at implementation time with a decision-log row. The source, selection rule, and thresholds are locked pre-event; only the literal encoding is deferred.

## 1. ADMET bioavailability rules (six rules)

| Rule ID | Rule | Threshold | Source citation | Prose description |
|---|---|---|---|---|
| ADMET-01 | Molecular weight | ≤ 500 Da | Lipinski et al., Adv. Drug Deliv. Rev. 1997 | Compounds above 500 Da are flagged as violating the molecular-weight rule. |
| ADMET-02 | logP | ≤ 5 | Lipinski et al., 1997 | Compounds with computed logP above 5 are flagged. |
| ADMET-03 | Hydrogen-bond donors | ≤ 5 | Lipinski et al., 1997 | Donor count above 5 is flagged. |
| ADMET-04 | Hydrogen-bond acceptors | ≤ 10 | Lipinski et al., 1997 | Acceptor count above 10 is flagged. |
| ADMET-05 | Rotatable bonds | ≤ 10 | Veber et al., J. Med. Chem. 2002 | Rotatable bond count above 10 is flagged. |
| ADMET-06 | Polar surface area | ≤ 140 Å² | Veber et al., 2002 | PSA above 140 Å² is flagged. |

`V` = number of the six rules violated, so `V ∈ {0,1,2,3,4,5,6}`. `bioavailability_score = max(0, 1 - 0.1 × V)`. `admet_flag = True` iff `bioavailability_score ≥ 0.55`. `admet_flag` is descriptive only and is not a scoring input. `final_score` uses `developability_score` only.

## 2. Developability flags

| Flag | Rule | Source citation | Prose description |
|---|---|---|---|
| Lipinski | MW ≤ 500, logP ≤ 5, HBD ≤ 5, HBA ≤ 10 | Lipinski et al., 1997 | Pass if all four Lipinski thresholds are met. |
| Veber | Rotatable bonds ≤ 10, PSA ≤ 140 Å² | Veber et al., 2002 | Pass if both Veber thresholds are met. |
| PAINS | RDKit FilterCatalog PAINS_A/B/C only | RDKit documentation; Baell & Holloway, J. Med. Chem. 2010 | Flag if the compound matches any PAINS_A, PAINS_B, or PAINS_C catalog entry. Prose classes include frequent-hitter, catechol, quinone, Michael-acceptor, and related assay-interference motifs. Exact SMARTS deferred to in-event. |
| Aggregator | Locked rule class list and source below | Shoichet lab / published aggregator rules; exact source recorded in decision log at implementation | Flag if the compound matches any locked aggregator rule class. Prose classes include long-chain lipophilic, planar aromatic stacking, and dye-like motifs. Exact SMARTS deferred to in-event. |

`developability_score = passes among Lipinski, Veber, PAINS, aggregator / 4`. PAINS compounds retained and flagged; score already penalizes them.

## 3. PAINS catalog description (prose only)

The PAINS catalog is limited to RDKit `FilterCatalog` PAINS_A, PAINS_B, and PAINS_C. No custom PAINS patterns are added. Prose pattern classes include: frequent-hitter scaffolds, catechol/quinone redox cyclers, Michael acceptors, stilbene-like and curcumin-like motifs, and other assay-interference classes defined by the RDKit catalog. The exact SMARTS strings are not stored here; they are committed in-event as `models/rules.json` with a decision-log row.

## 4. Aggregator rule description (prose only)

The aggregator rule class list is locked to published aggregator motifs; the exact source citation and SMARTS are committed in-event in `models/rules.json` with a decision-log row. Prose classes include long-chain lipophilic compounds, planar aromatic stacking motifs, and dye-like scaffolds. No project-specific aggregator rules are invented.