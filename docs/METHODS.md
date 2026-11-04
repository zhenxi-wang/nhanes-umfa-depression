# Computational methods

## Analytic sample

Five cycles are pooled: 2011-2012, 2013-2014, 2015-2016, 2017-March 2020, and August 2021-August 2023. Eligibility requires age at least 18 years, no recorded pregnancy at examination, complete PHQ-9, observed 5-MTHF and UMFA, and a positive appropriate examination/folate/phlebotomy weight. The verified analytic sample is 24,326.

The release's substituted below-LOD values are retained. Exposures are not imputed. Covariate harmonisation follows `code/pool.py` and `code/common.py`; this repository reproduces those operational definitions rather than claiming that missing-at-random or causal assumptions have been established.

## Primary model

The outcome is PHQ-9 >=10. UMFA, 5-MTHF and MeFox enter the model together, each as a log concentration standardised within cycle. `code/common.py::X_model` defines the demographic, socioeconomic and extended covariate blocks. Model 3 is the primary extended model; Model 5 is the DAG-guided reduced adjustment set.

## Multiple imputation

The current implementation uses predictive mean matching with five observed donors and three stochastic chains. After 20 burn-in sweeps, one dataset is retained every five sweeps through sweep 70. Each chain contributes ten datasets, for 30 in total. Seeds are 20261104, 20261105 and 20261106.

Multicategory targets are quantified using the first canonical variate, then imputed from observed category donors. Continuous and binary targets also use observed donors. The implemented procedure is Python chained-equations PMM; it is not an execution of the R `mice` package. Design predictors include log pooled weights and full PSU fixed effects; PSUs are nested within stratum and cycle. Ordinary covariate missingness is distinguished from cycle-level structural non-availability.

The imputation models include the three primary exposures, primary depressive-symptom outcome, all Model 3 covariates and design predictors. Supplemental folic acid is set to zero for confirmed non-users; unknown supplement-use status produces an unknown amount. The logical zero for imputed non-users is enforced.

Discovery, earlier-cycle and single-cycle imputations are performed separately. Medication imputations include medication indicators and are restricted to the four cycles with drug-level data. Secondary-outcome imputations additionally include the secondary outcomes as auxiliary predictors. Subgroup membership is based on observed values; records missing the grouping variable are not assigned a subgroup by imputation.

## Survey estimation and pooling

Weights are WTMEC2YR for 2011-2016, WTFOLPRP for 2017-March 2020, and WTPH2YR for August 2021-August 2023. Pooled weights multiply the cycle weight by its duration divided by 11.2 years. Taylor linearisation aggregates score contributions by masked PSU within stratum. Domain models retain zero score contributions outside the domain.

Within-imputation coefficient estimates and covariance matrices are combined using Rubin's rules. Scalar confidence intervals use Barnard-Rubin degrees of freedom, with the relevant complete-data survey-design degrees of freedom. The main sample has 87 such degrees of freedom. Multi-parameter tests use D1 pooling with Reiter finite-sample correction. `analysis/analysis_runtime.py` contains the explicit formulas.

The formal 5-MTHF-minus-UMFA contrast is formed within each imputation using the full coefficient covariance matrix. Variances include both diagonal terms and the covariance term. Averaged ORs or rounded manuscript ORs are not used to calculate this contrast.

## Non-linear and categorical exposure models

UMFA spline knots are the 5th, 35th, 65th and 95th percentiles of observed log-UMFA. The reference is the 10th percentile, 0.348 nmol/L. Both revised Figure 2 panels use the same mutually adjusted spline model. Marginal prevalences are pooled on the logit scale; the paired absolute difference is pooled on its probability scale using the within-imputation covariance of the contrast.

Cycle-specific quartiles are based on observed concentrations. Both single-biomarker quartiles and a matched UMFA quartile model additionally adjusting for 5-MTHF and MeFox are fitted. Strict greater-than tail counts use raw UMFA values, not transformed values.

## Multiplicity

The historical screen records 96 tests. The manuscript's biomarker family contains the identified 82 biomarker/examination comparisons; the remaining 14 functional or questionnaire comparisons are disclosed separately. Both 82-test and all-96-test BH adjustments are available.

The revised exploratory family retains 26 tests: six secondary UMFA outcomes, four inferences for each of four quartile exposures, and four upper-tail threshold contrasts. Unrounded P values are used. The matched mutually adjusted quartile model is treated as a robustness check of the primary exposure specification.

## Validation scope

Arithmetic validation of pooling, numerical optimiser agreement and comparison with stored reference results do not prove the imputation assumptions, missing-at-random mechanism, absence of residual confounding, or causal identification. The historical 10-imputation sensitivity analysis remains in the legacy code and is not the current primary imputation procedure.

Method references: [PMM and factor quantification](https://amices.org/mice/reference/mice.impute.pmm.html), [Rubin/Barnard-Rubin pooling](https://amices.org/mice/reference/pool.html), and [D1/Reiter reference implementation](https://github.com/cran/mitml/blob/master/R/internal-pool.R).
