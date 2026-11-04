# Packaging changes and provenance

The original project files, manuscript documents, raw data and figures were not edited. This repository candidate contains copies and portable replacements.

## Shared and historical source

- The original statistical and figure scripts are retained in `code/`. Data locations now honour `NHANES_DATA_DIR`; historical figure locations honour `NHANES_LEGACY_FIGURES`.
- Directory creation was removed from shared-helper import side effects. The top-level runner creates output directories explicitly.
- The original screening model functions were extracted unchanged into `code/screening_model.py`, allowing the screening audit to reuse the same implementation. The 96 original screening tests still reproduce the recorded unrounded results.
- Historical results remain labelled as historical simple-imputation analyses. Original labels such as "primary" in those scripts do not identify the current MICE primary model.
- Manuscript assembly and editorial audit utilities were excluded because they are not needed to reproduce statistical results and can contain draft manuscript material.

## Current revised workflow

- Revised scripts were renamed descriptively and moved to `analysis/`.
- Computer-specific paths and dynamic runtime loading were replaced with imports and paths relative to this repository.
- The previously separate exclusion of all recorded psychotropic medication users was integrated into the medication-model loop.
- Figure font lookup is configurable. Arial is used by default but is not bundled or redistributed.
- The analysis settings, random seeds, statistical formulas and archived aggregate results were preserved. The 30 primary imputations and main coefficient/contrast models were rerun to check numerical equivalence after packaging.

## Added reproducibility utilities

The package adds a stage runner, public-data downloader, raw-file manifest, data preparation wrapper, result-table exporter, focused statistical tests, primary-result validator, reference comparison, screening diagnostic, and English documentation. `config/source_provenance.json` identifies copied and newly written files.

The historical seven-row final-eligibility folate re-screen has only an aggregate reference in the supplied project. Its exact generating script was not recovered. An optional diagnostic implements the recoverable screening specification and reports discrepancies; it does not silently replace the historical table. See `REPRODUCIBILITY_STATUS.md`.

## Distribution boundary

The package contains no raw XPT files, participant identifiers or records, imputed participant values, manuscript drafts, personal credentials, or font binaries. Aggregate model summaries are retained. Temporary validation outputs were removed before creating the distributable ZIP. `FILE_MANIFEST.csv` records the size and SHA-256 hash of each distributed file other than the manifest itself.

The package is published at https://github.com/zhenxi-wang/nhanes-umfa-depression following the corresponding author's explicit upload request. Documentation and citation metadata now include this verified URL. Selecting a software licence and creating a release DOI remain author decisions. The original project files remain unchanged.
