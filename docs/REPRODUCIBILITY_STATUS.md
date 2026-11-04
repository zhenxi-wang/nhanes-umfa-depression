# Reproducibility status

## Checks performed on this portable package

Validation used Python 3.12.10 on Windows and the exact versions in `requirements.txt`.

1. The packaged data builders were run against the public XPT files in a separate read-only source directory. Derived files were written into a new validation directory. The analytic sample was 24,326, with cycle counts 4,529; 5,120; 4,820; 4,829; and 5,028.
2. The packaged primary MICE procedure was run afresh from that rebuilt pooled frame, using all 30 completed datasets.
3. All 30 primary survey-weighted models and their formal coefficient contrast were refitted. The primary UMFA OR was 1.1155211654759016, with 95% CI 1.040458407425754-1.195999245855017 and P=0.0024683243624240158. The formal ratio was 0.876071322520171, with P=0.06848228602059105. These matched the archived revised results.
4. The 96 historical screen tests were regenerated from the rebuilt discovery frame. The maximum absolute difference from the recorded unrounded P values was approximately 4.1e-14.
5. Upper-tail counts, raw missingness, strict complete-case sensitivity, original simple-imputation sensitivity, and both screening BH families were recalculated.
6. The six focused statistical-helper tests passed.
7. Every Python file was syntax-checked, all stage paths were checked by the dry-run runner, and source/documentation files were checked for local author paths and non-English editorial material.
8. English table exports and portable figure rendering were exercised using explicitly staged aggregate reference summaries, together with the freshly rebuilt analytic sample. This checks packaging/export compatibility; it is not a second full refit of every model.

The complete revised workflow, including partition-specific and outcome-specific imputations, model fits and figures, was executed during the preceding statistical revision. Packaging preserves that numerical code apart from documented portability changes and integration of an additional psychotropic-exclusion calculation. A complete second execution of every expensive model and every historical sensitivity analysis was not performed during packaging.

## Boundaries

- The downloader's source URL pattern was checked against official CDC documentation, and local input hashes were recorded. The package was not used to re-download all 129 public files during validation. CDC availability and future public-release changes remain external dependencies.
- Statistical validation is not evidence that missing-at-random assumptions, covariate definitions, or causal identification are correct.
- Some supplementary analyses remain historical simple-imputation sensitivity analyses. Their provenance is explicit in `docs/ANALYSIS_MAP.md` and table-export method labels.
- The historical seven-row folate re-screen on final discovery eligibility is retained as `reference_results/legacy/screen_folate_final.json`, but its exact generating script was not found among the available source files. Applying the recoverable screening model to final eligibility does not exactly reproduce that summary. The optional `audit-screen-final` stage records the differences without changing the reference. This historical subtable should not be described as independently reproduced; it requires author verification before submission. The original 96-test screen, including the 82-test biomarker family, is fully reproduced by `code/screen.py` and is unaffected by this limitation.
- Exact floating-point results can vary slightly with numerical libraries and platforms. The comparison script uses documented tolerances; model definitions, sample sizes and substantive conclusions should be checked if differences exceed them.
- Arial is a system dependency for matching the submitted typography. Font binaries are not included.

The source repository is https://github.com/zhenxi-wang/nhanes-umfa-depression. The validation scope described above is unchanged by repository publication. No release DOI has been assigned.
