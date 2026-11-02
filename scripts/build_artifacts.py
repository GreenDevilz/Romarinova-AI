# Placeholder stub.
# No project-specific logic before 2 Nov 2026.
# Implemented in-event.
#
# Generator contract (documentation only; no executable logic pre-event):
#
# Single source of truth: all generated artifacts in this repository are
# produced by this script from one results store. Hand-editing generated
# CSVs is prohibited; if a value is wrong, fix the generator and regenerate.
#
# Artifacts produced in-event:
#   results/evaluation_table.csv
#   results/screened_compounds_annotated.csv
#   results/primary_shortlist.csv                 (includes header metadata block)
#   results/primary_shortlist_compound_level.csv  (defines U)
#   results/diagnostics.csv
#
# results/diagnostics.csv canonical `section` values:
#   data_feasibility_ledger
#   missing_structures
#   scoring_coverage_gap
#   filtered_but_interesting
#   suspected_variety_links
#   sensitivity
#   disagreement
#   missing_ambiguous_structures
#   other
#
# Required header metadata fields (protocol.md section 2.1), all unique-compound counts:
#   N_star, U_star, known_active_count, high_ml_ad_extrapolation_count,
#   ad_extrapolation_count, ad_missing_count
#
# Every N term and flag emitted here must resolve to a definition in
# protocol.md section 2.1 or 2.2. Any new term requires a protocol amendment,
# not a generator change.

# Placeholder stub.