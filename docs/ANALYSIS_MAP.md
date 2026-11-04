# Manuscript-to-code map

| Manuscript component | Current or retained approach | Main files |
| --- | --- | --- |
| Analytic eligibility and variable derivation | Shared data preparation | `scripts/prepare_analysis_data.py`, `code/pool.py`, `code/common.py` |
| Initial exploratory biomarker screen, S8 | Historical screening model; exact P values | `code/prep.py`, `code/screen.py`, `analysis/audit_counts_and_multiplicity.py` |
| Table 1, S1, S2 and descriptive comparisons | Observed-data descriptives | `code/a0_desc.py`, `code/a4_tables.py`, `code/flowfig.py` |
| Table 2A, primary three-species models | MICE, 30 datasets | `analysis/impute_chained_equations.py`, `analysis/fit_primary_models.py` |
| Formal 5-MTHF versus UMFA comparison | Covariance-aware MICE contrast | `analysis/analysis_runtime.py`, `analysis/fit_primary_models.py`, `analysis/estimate_curves_and_validate.py` |
| Table 2B and Figure 1a | Single-biomarker MICE models | `analysis/fit_primary_models.py`, `analysis/make_main_figures.py` |
| Figure 1b | Observed survey-weighted prevalence | `analysis/make_main_figures.py` |
| Table 3 | Single-biomarker MICE quartile models | `analysis/fit_primary_models.py` |
| Figure 2, spline points and adjusted prevalence | One shared three-species MICE spline model | `analysis/fit_primary_models.py`, `analysis/estimate_curves_and_validate.py`, `analysis/make_main_figures.py` |
| Figure 3 and S13 | Separately imputed discovery, replication and cycles | `analysis/impute_analysis_samples.py`, `analysis/fit_replication_and_sensitivity.py`, `analysis/make_main_figures.py` |
| S4, most revised sensitivity rows | MICE; subgroup masks based on observed data | `analysis/fit_replication_and_sensitivity.py` |
| S4, medication rows | Medication-specific MICE restricted to earlier cycles | `analysis/impute_auxiliary_models.py`, `analysis/fit_replication_and_sensitivity.py` |
| S4, strict complete case and original simple imputation | Explicit alternative missing-data approaches | `analysis/audit_counts_and_multiplicity.py` |
| S4, hsCRP/MCV and self-rated-health rows retained from the original analysis | Historical simple-imputation sensitivity analyses | `code/a2_mech.py`, `code/a4_tables.py`, `code/a5_review.py` |
| S3 subgroups; S5 residual UMFA; S6 phenotype; S7 B12 | Retained historical sensitivity/exploratory analyses | `code/a2_mech.py`, `code/a3_extra.py`, `code/a5_review.py` |
| S9 missingness | Pre-imputation observed/derived covariate missingness | `analysis/audit_counts_and_multiplicity.py` |
| S10 regression validation under the original specification | Historical simple-imputation estimator checks | `code/a5_review.py` |
| S10 new MICE pooling/optimiser validation | Revised arithmetic and numerical checks | `analysis/estimate_curves_and_validate.py`, `analysis/validate_primary_model.py`, `tests/test_statistical_helpers.py` |
| S11 | Revised 26-test BH family | `analysis/fit_secondary_outcomes.py` |
| S12, S16-S18 descriptive/correlation results and historical design diagnostics | Observed-data comparisons; historical design diagnostics explicitly labelled | `code/a4_tables.py`, `code/a6_review2.py` |
| S14A weighted standardisation | Revised MICE | `analysis/fit_primary_models.py` |
| S14B weighted quartiles and thresholds; S14C thresholds | Revised MICE plus observed descriptives | `analysis/fit_primary_models.py`, `analysis/export_result_tables.py` |
| S14D upper-tail sample and case counts | Observed raw UMFA and PHQ-9 | `analysis/audit_counts_and_multiplicity.py` |
| S14E matched mutually adjusted quartiles | Revised MICE robustness model | `analysis/fit_primary_models.py` |
| S15 non-linear confounder sensitivity | Retained historical simple-imputation model | `code/a6_review2.py` |
| S19 severity thresholds and S22 secondary outcomes | Secondary-outcome-specific MICE | `analysis/impute_auxiliary_models.py`, `analysis/fit_secondary_outcomes.py` |
| S20 variable definitions | Shared derivations and documented labels | `code/pool.py`, `code/common.py`, `code/tables.py` |
| S21 six-species exploration | Retained historical simple-imputation analysis | `code/a7_species.py` |
| Figures S1-S5 | Historical descriptives, sensitivity models and conceptual drawings | `code/flowfig.py`, `code/figs.py`, `code/mechfig.py`, `code/dagfig.py`, `code/figstyle.py` |

Historical labels such as "primary" inside original scripts refer to the original specification. They must not be substituted for the revised MICE primary result. Historical figure numbering in `code/figs.py` predates the manuscript's current Figure 1/Figure 2 order; the three current main figures are generated only by `analysis/make_main_figures.py`.

Manuscript assembly scripts, confidential draft text, Word comment files and editorial revision instructions are not needed to reproduce statistical results and are not part of this package.
