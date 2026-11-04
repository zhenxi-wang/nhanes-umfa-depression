# Code inventory

This package contains 39 Python source files. All paths below are relative to the repository root. Start with `reproduce.py` and `README.md`.

| File | Purpose |
| --- | --- |
| `analysis/analysis_runtime.py` | Shared model fitting, Rubin pooling, Barnard-Rubin degrees of freedom, covariance-aware contrasts, D1 tests, restricted cubic splines and BH adjustment. |
| `analysis/audit_counts_and_multiplicity.py` | Recalculate raw upper-tail counts, depressive-symptom cases, weighted exposure shares, missingness, complete cases, simple imputation and BH checks. |
| `analysis/audit_final_eligibility_screen.py` | Optional diagnostic for the historical seven-row final-eligibility screen; preserves the recorded table and reports discrepancies. |
| `analysis/compare_reference_results.py` | Compare regenerated aggregate revised results with archived reference summaries using explicit tolerances. |
| `analysis/estimate_curves_and_validate.py` | Calculate spline contrasts, adjusted prevalence and absolute differences, E-values, independent optimisation and pooling checks. |
| `analysis/export_result_tables.py` | Export 27 English CSV tables and a combined Markdown document with explicit missing-data provenance. |
| `analysis/fit_primary_models.py` | Fit current primary, single-biomarker, quartile, spline, threshold, standardisation and DAG models. |
| `analysis/fit_replication_and_sensitivity.py` | Fit current discovery/replication, cycle-specific, subgroup, exclusion and medication models; calculate heterogeneity. |
| `analysis/fit_secondary_outcomes.py` | Fit secondary outcomes and recalculate the defined 26-comparison exploratory BH family. |
| `analysis/impute_analysis_samples.py` | Run independent imputation for discovery, earlier-cycle and individual-cycle samples. |
| `analysis/impute_auxiliary_models.py` | Generate medication-specific and secondary-outcome-specific imputations. |
| `analysis/impute_chained_equations.py` | Generate 30 completed datasets using three seeded chained-equations predictive-mean-matching chains. |
| `analysis/make_main_figures.py` | Draw revised Figures 1-3 and export high-resolution PNG/TIFF plus vector PDF/SVG. |
| `analysis/validate_primary_model.py` | Refit all 30 primary models and compare full-precision ORs, confidence intervals, P values and the formal contrast with references. |
| `code/a0_desc.py` | Historical observed-data cycle descriptives and biomarker/outcome summaries. |
| `code/a1_main.py` | Historical single-biomarker, mutually adjusted, quartile and secondary-outcome models under the original specification. |
| `code/a2_mech.py` | Historical spline, phenotype, effect-modification and sensitivity analyses. |
| `code/a3_extra.py` | Historical reverse-causation checks, residual UMFA, B12 interaction, margins and medication outcomes. |
| `code/a4_tables.py` | Historical Table 1, distribution/correlation summaries and additional sensitivity checks. |
| `code/a5_review.py` | Historical DAG/cycle/supplement/fasting checks, 10-imputation sensitivity, FDR and regression validation. |
| `code/a6_review2.py` | Historical weighted-quantile/SD, nonlinear-confounder, outcome-selection and cycle-heterogeneity checks. |
| `code/a7_species.py` | Historical all-six-folate-species detection, discovery, replication and mutually adjusted exploratory analyses. |
| `code/common.py` | Define analytic eligibility, exposure transformations, model covariates, weighted summaries and the survey design. |
| `code/dagfig.py` | Draw the directed acyclic graph for the reduced adjustment set. |
| `code/fig3_replication.py` | Draw the original simple-imputation replication/cycle forest plot. |
| `code/figs.py` | Draw original main and supplementary statistical figures; historical main-figure numbering differs from the revised manuscript. |
| `code/figstyle.py` | Shared typography, nodes and arrows for conceptual supplementary diagrams. |
| `code/flowfig.py` | Recalculate sample-exclusion counts and draw the participant flow diagram. |
| `code/mechfig.py` | Draw the conceptual biological-pathway diagram; this is not an observed-data result. |
| `code/pool.py` | Read and harmonise five NHANES cycles, folate weights, PHQ-9, covariates and medication indicators. |
| `code/prep.py` | Prepare the August 2021-August 2023 wide frame used by the initial exploratory screen. |
| `code/screen.py` | Regenerate the original 96 exploratory screening tests, including the 82-test biomarker/examination family. |
| `code/screening_model.py` | Reusable original age/sex/race/education/income-adjusted screening model. |
| `code/survey.py` | Implement weighted logistic/linear regression, Taylor-linearised covariance, domain inference and Wald contrasts. |
| `code/tables.py` | Export the original main and supplementary tables as Markdown; these retain historical analysis labels. |
| `reproduce.py` | Command-line entry point; runs the complete revised pipeline or one named stage. |
| `scripts/download_public_data.py` | Download 129 public CDC XPT inputs; check SAS XPORT headers and optional SHA-256 fingerprints. |
| `scripts/prepare_analysis_data.py` | Read raw XPT files and build pooled/screening frames in a separate output location; verify analytic eligibility. |
| `tests/test_statistical_helpers.py` | Six focused numerical tests: contrast covariance, imputation variance/df, BH ordering, joint null and spline tails. |

## Supporting files

| File or folder | Purpose |
| --- | --- |
| `requirements.txt` | Exact validated dependency versions. |
| `config/raw_data_manifest.csv` | Public input filenames, source URLs, bytes and SHA-256 hashes. |
| `config/analysis_specification.json` | Analytic definitions, imputation settings and multiplicity-family settings. |
| `config/source_provenance.json` | Origin and packaging status of each source file. |
| `reference_results/legacy/` | Aggregate historical references, including unrounded initial screening P values. |
| `reference_results/revised/` | Aggregate references from the current 30-imputation workflow. |
| `docs/METHODS.md` | Computational definitions and statistical methods. |
| `docs/ANALYSIS_MAP.md` | Manuscript tables and figures mapped to scripts and methods. |
| `docs/REPRODUCIBILITY_STATUS.md` | Validation performed and explicit reproduction limits. |
| `docs/PACKAGING_CHANGES.md` | Changes made to copies for portability and distribution. |
| `docs/AI_ASSISTANCE.md` | AI assistance used for this code package. |
| `docs/validation_summary.json` | Aggregate machine-readable record of packaging checks. |
| `FILE_MANIFEST.csv` | Distributed filenames, bytes and SHA-256 hashes; excludes itself. |
| `CITATION.cff` | Citation metadata; no invented DOI or repository URL. |
| `LICENSE_DECISION.md` | Authors must select a software licence before granting reuse permissions. |
| `.gitignore` | Exclude participant-level data, working outputs, environments and credentials. |

## What is deliberately excluded

Raw/derived participant data, imputed values, Word/PDF manuscript drafts, editorial comments, temporary validation files, credentials and font binaries are excluded. Raw public inputs can be retrieved using the downloader. Temporary packaging utilities are not required to reproduce the analyses.

The historical seven-row final-eligibility re-screen has no recovered original generating script. Its reference is preserved and a diagnostic is provided; this limitation is not presented as successful reproduction.
