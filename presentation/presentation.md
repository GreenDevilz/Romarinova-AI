# Presentation

## 1. Deck outline

### Slide 1: Title
**Romarinova AI**
*From Tunisian Rosemary Chemotypes to Evidence-Guided Candidate Prioritization for α-Glucosidase Inhibition*
EthnoHack 2026 · Track 2 — Natural Products and Drug Discovery
*"A prioritization and traceability system, not a claim that rosemary cures diabetes."*

### Slide 2: Problem and motivation
- Tunisian rosemary has distinct chemotypes (var. typicus, var. troglodytorum).
- α-Glucosidase inhibition is a known research strategy for type 2 diabetes.
- **Fragmentation:** Essential-oil chemistry is not linked to enzyme-level bioactivity; bioactivity is not linked to variety-level sources.
- **Our contribution:** variety-level traceability, scaffold-aware validation, and a side-by-side ML-independent comparison view.

### Slide 3: Data and approach
- **Bioactivity:** ChEMBL 37 human MGAM (`tid=1`). 746 raw records → 179 eligible (31 active, 148 inactive). 108 ambiguous excluded.
- **Rosemary catalogue:** 64 records across 28 distinct compounds. Linked to primary LC-MS publications.
- **Split:** Document-aware (zero shared documents across train/val/test). Train 128 (14 actives), Val 30 (10 actives), Test 21 (7 actives).

### Slide 4: ML screening and evaluation
**Model:** Random Forest (100 trees, depth 5) on ECFP4. Selected via class-aware 4-fold CV.
| Metric | Test | Validation |
|---|---|---|
| PR-AUC | **0.781** | 0.981 |
| ROC-AUC | **0.760** | 0.990 |
| N actives | 7 | 10 |

- **Small-N flags:** `ci_not_reported = True` (active count < 20).
- **Similarity-only baseline:** PR-AUC 0.735. The ML model adds a modest +0.046 uplift.

### Slide 5: Similarity, developability, and scoring
- **Similarity:** max Tanimoto to training actives.
- **Developability:** Lipinski, Veber, PAINS, Aggregator.
- **Locked 4-component score:** `0.40 × ML + 0.25 × sim + 0.25 × dev + 0.10 × novelty`.
- **Upstream verification:** AD distances cross-checked between `02` and `03` (max diff 1e-16).

### Slide 6: Shortlist and compound-level view
**Primary shortlist (N=10 pairs, U=5 compounds):**
1. p-Cymene (0.640)
2. Verbenone (0.636)
3. γ-Terpinene (0.633)
4. α-Terpinolene (0.633)
5. α-Pinene (0.632)
- **Diagnostics hold:** Carnosol, carnosic acid, rosmarinic acid (failed developability).
- `known_active_count = 0`.

### Slide 7: Side-by-side locked vs ML-independent comparison
- **Pair-level overlap:** **10/10** → sensitivity label = **`robust`**.
- **The finding:** 27 of 28 screened compounds (96.4%) share `ML_rank_score = 0.633333`. 
- The model produces one prediction for the extrapolation regime. Zeroing ML changes nothing because ML was a constant, not a signal.
- The rank order is determined entirely by `similarity_support`.

### Slide 8: Diagnostics and limitations
- **AD Extrapolation:** 28/28 rosemary compounds (threshold 0.7330).
- **Known gaps:** Docking specified but not executed; web prototype schemas populated but pages pending.
- **Method limitations:** Human MGAM-specific; ECFP4 only; PAINS are heuristic flags.

### Slide 9: Impact and responsible research
- **Delivered:** 5 developable candidates; diterpenes preserved in diagnostics; reproducible open pipeline.
- **Responsible research:** No clinical claims. Experimental and computational outputs kept separate. PAINS/ADMET are heuristic flags.

### Slide 10: Conclusion and future work
- A leakage-free, document-aware ML pipeline with honest small-N reporting.
- A shortlist carried entirely by similarity support in the extrapolation regime.
- An ML-independent comparison view that makes the collapse visible.
- **Future work:** Execute docking; implement web prototype; more independent training actives.

## 2. Narration outline
1. **Opening:** Rosemary is fragmented across varieties and papers. We trace every ranking to its source.
2. **Data:** 179 eligible ChEMBL compounds. Document-aware split (zero leakage). 28 rosemary compounds.
3. **Results:** Test PR-AUC 0.781. But CI is wide due to N=7. Similarity alone gets 0.735.
4. **Demo:** Side-by-side view shows 10/10 overlap. Why? 27/28 compounds share the exact same ML score. ML is silent on extrapolation.
5. **Impact:** 5 monoterpenes prioritized. Diterpenes saved in diagnostics. The ranking is carried by similarity and developability.

## 3. Narrative text

### 3.1 ML contribution statement
> "ML achieves PR-AUC 0.781 on test. But 27 of 28 screened rosemary compounds receive the identical `ML_rank_score` of 0.633333. On the ML-independent comparison view, pair overlap is 10/10. The shortlist is invariant to ML weight removal because ML was never load-bearing. This is a domain-shift disclosure rendered as a UI element."

### 3.2 Modest-shortlist framing
The five shortlisted compounds are all small monoterpenes that pass all four developability checks. The characteristic diterpenes fail PAINS/Veber and are preserved in the diagnostics appendix for experimental follow-up. This is prioritization discipline.

### 3.3 ML / AD caveat
`ML_rank_score` is a validation-rank percentile, not a probability. All 28 rosemary compounds are outside the applicability domain. The honest collapse figure is `ml_rank_top_band_share = 0.9643`.

### 3.4 Closing disclaimer
> Romarinova AI does not claim rosemary cures diabetes; it turns fragmented rosemary chemotype chemistry into a reproducible, evidence-traced prioritization system that reports where its own ML contribution ends.