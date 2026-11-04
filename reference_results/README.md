# Aggregate reference results

`revised/` contains aggregate outputs from the current 30-imputation analyses. `legacy/` preserves aggregate outputs from the historical simple-imputation analyses and the unrounded initial screen.

These files contain model summaries and grouped descriptive results, not participant-level records. They allow numerical validation and provide reference inputs for testing export code. Copying them into an output directory is not equivalent to rerunning the statistical models.

The main pipeline fits models from prepared public data and writes new results to `outputs/`. `analysis/compare_reference_results.py` compares those new outputs with these references. The table exporter uses explicitly labelled historical references for the few retained legacy rows.

The historical `legacy/screen_folate_final.json` has no recovered generating script and should not be described as independently reproduced. The `audit-screen-final` stage documents reanalysis under the recoverable screening model. See `docs/REPRODUCIBILITY_STATUS.md`.
