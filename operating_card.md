# Operating Card — Romarinova AI

**One page. Daily execution. Hard checkpoints, 6 stop rules, 3 lock rules.**

## Hard checkpoints

| Checkpoint | When | Owner | Pass condition |
|---|---|---|---|
| Data feasibility ledger complete | End of Day 1 | P1/P2 | `results/diagnostics.csv` `data_feasibility_ledger` rows populated: ChEMBL counts, scaffold groups, pair counts, structure counts, overlap count, eligible compound N |
| Data freeze | End of Day 1 | P1 | Curated CSV, splits, CV folds, InChIKey check saved; no ad-hoc variants after freeze |
| Model freeze | End of Day 2 | P1 | Selected model, `models/model_config.json`, validation PR-AUCs, test table, AD for all eligible compounds saved |
| UI contract freeze | End of Day 2 | P3 | Page layout, views, tooltips, column sets fixed alongside `web/data_contract.json` and `web/ui_contract.json` |
| Ranking lock | Day 3 | P2 | Compound identities, `final_score` rankings, duplication cap, AD flags, compound-level view, ligand-only ranking fixed |
| Final record verification | Day 4 | P2 | Locked shortlist, compound-level view, diagnostics data, generated artifacts verified for consistency |
| Submission cutoff | 8 Nov 2026 | P3 | Devpost final, repo, video, deck submitted |

## 6 stop rules

1. **Data-freeze stop:** no curation changes after end of Day 1 without a decision-log row and P1 approval.
2. **Model-freeze stop:** no model/tuning changes after end of Day 2 without a decision-log row and P1 approval.
3. **UI-freeze stop:** no UI redesign after end of Day 2; Day 3–5 is data wiring only.
4. **Ranking-lock stop:** no score, rank, duplication-cap, or AD-flag changes after the Day 3 ranking lock without a P1/P2 decision-log entry and final verification rerun.
5. **Documentation stop:** documentation is confined to the Day 6 morning slot; review cutoff applies.
6. **Submission stop:** no changes after the 8 Nov deadline; final checks only on Day 7.

## 3 lock rules

1. **Scoring lock:** `final_score` weights 0.40/0.25/0.25/0.10 and ligand-only weights 0.00/0.60/0.30/0.10 are fixed a priori; not tuned on results.
2. **AD lock:** AD descriptor (ECFP4 2048-bit, Jaccard, k=5, training-LOO 95th percentile primary threshold) is fixed before scoring; AD does not modify `final_score`.
3. **Shortlist lock:** up to 10 compound-variety pairs; duplication cap two per compound; identified pairs only; `structure_status = defined`; tie-break higher `similarity_support`, then higher `developability_score`, then lexicographic InChIKey, then lexicographic `cultivar_source`.

## Daily discipline

- P1 single technical decision-maker; P2 deputy.
- P2 protected block before ranking lock.
- P3 daily stop; single-tasking; video handoff to P1; deck assembly from pre-written outline.
- Day 3 morning dry-run with empty templates.
- Decision log ≤20 rows with pre-assigned row budgets.
- `scripts/build_artifacts.py` is the single generator; hand-editing generated CSVs is prohibited.