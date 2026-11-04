# NHANES UMFA and depressive symptoms

Source repository: [https://github.com/zhenxi-wang/nhanes-umfa-depression](https://github.com/zhenxi-wang/nhanes-umfa-depression).

Reproducibility code for **Circulating unmetabolised folic acid and depressive symptoms among US adults: a cross-sectional analysis of NHANES 2011-2023**.

The primary analysis uses 30 multiply imputed datasets and a survey-weighted, mutually adjusted model of UMFA, 5-MTHF and MeFox. The package also preserves the original simple-imputation analyses used for explicitly identified sensitivity and exploratory results.

## Start here

- `CODE_INVENTORY.md`: every code file, its purpose, and the corresponding analysis.
- `docs/ANALYSIS_MAP.md`: manuscript tables and figures mapped to scripts and missing-data approaches.
- `docs/METHODS.md`: model definitions, MICE settings, pooling, thresholds, and multiplicity families.
- `docs/REPRODUCIBILITY_STATUS.md`: checks actually performed and remaining platform-dependent limitations.
- `config/raw_data_manifest.csv`: 129 public XPT inputs, source URLs, source sizes, and SHA-256 fingerprints.
- `reference_results/`: aggregate reference results; no participant-level records.

## Environment

Use Python 3.12.10 to match the validation environment. Install the recorded package versions:

```text
python -m venv .venv
```

Activate the virtual environment, or use its Python executable, then install dependencies:

```text
python -m pip install -r requirements.txt
```

For example, the environment's Python executable is `.venv/Scripts/python.exe` on Windows or `.venv/bin/python` on macOS/Linux. All commands below are run from this repository's top-level folder. Paths are resolved relative to the package; no original author's computer path is required.

## Obtain the public data

```text
python reproduce.py --stage download
```

Alternatively, place the XPT files named in `config/raw_data_manifest.csv` in `data/`. Downloading may take time and depends on CDC availability. The downloader rejects HTML/error pages masquerading as XPT files. Recorded hashes identify the release used for the analyses; a changed public file is reported rather than silently treated as identical.

The source pages include the [2021-2023 folate documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/FOLFMS_L.htm) and the [2017-March 2020 folate documentation](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_FOLFMS.htm). The full file inventory is in the manifest.

## Run the current analysis

```text
python reproduce.py --stage all
```

This prepares the pooled and screening frames, regenerates the historical screening results, performs the current imputations, fits the revised models, validates pooling, exports tables, and draws Figures 1-3. Outputs are written to `outputs/revised/`; the screen is written to `outputs/legacy/`. The full pipeline is substantially longer than an individual stage.

To inspect commands without executing them:

```text
python reproduce.py --stage all --dry-run
```

The figures use Arial by default. Arial font files are not distributed. On a machine without installed Arial, provide a directory containing the font files:

```text
python reproduce.py --stage figures --font-dir /path/to/arial/fonts
```

For a deliberate alternative font, use `--font-family "DejaVu Sans"`. The statistical results are unaffected; typography and layout can differ.

## Read an existing raw-data folder without changing it

```text
python reproduce.py --stage prepare --source-data /path/to/raw_xpt_files --data-dir outputs/local_run/data --output-dir outputs/local_run/revised
python reproduce.py --stage impute-main --data-dir outputs/local_run/data --output-dir outputs/local_run/revised
python reproduce.py --stage validate-primary --data-dir outputs/local_run/data --output-dir outputs/local_run/revised
```

The preparation command reads XPT files from `--source-data` and writes the new pickle files only to `--data-dir`. Use the same `--data-dir` and `--output-dir` in subsequent stages.

## Individual stages

| Stage | Function |
| --- | --- |
| `download` | Download public XPT inputs from the manifest. |
| `prepare` | Build `pooled.pkl` and `wide.pkl` from raw inputs. |
| `screen` | Regenerate the original 96 recorded screening tests, including the 82-test biomarker subset. |
| `impute-main` | Produce the primary 30 completed datasets. |
| `impute-samples` | Independently impute discovery, earlier-cycle and individual-cycle samples. |
| `impute-aux` | Produce medication-specific and secondary-outcome imputations. |
| `core` | Fit primary, single-biomarker, quartile, spline, threshold and standardisation models. |
| `replication` | Fit replication, cycle-specific, medication and other sensitivity models. |
| `secondary` | Fit secondary outcomes and recalculate the defined 26-comparison BH family. |
| `margins` | Calculate spline contrasts, adjusted prevalence, E-values and independent validation. |
| `audit` | Check upper-tail counts, raw missingness, complete cases and multiplicity calculations. |
| `tables` | Export English CSV and Markdown result tables. |
| `figures` | Export revised Figures 1-3 as PNG, TIFF, PDF and SVG. |
| `validate` | Compare regenerated summaries against aggregate revised references. |
| `validate-primary` | Refit all 30 primary completed datasets and check the main ORs and formal contrast. |
| `audit-screen-final` | Diagnose the unrecovered historical seven-row final-eligibility re-screen; does not replace its stored reference. |
| `legacy` | Re-run the historical simple-imputation analyses and supplementary graphics. |

`legacy` is explicit and is not part of `all`. Historical results must retain their original method labels. The supplementary models intentionally retained under simple imputation are mapped in `docs/ANALYSIS_MAP.md`.

## Statistical checks

```text
python -B -m unittest discover -s tests -v
```

The focused tests address covariance-aware contrasts, within/between-imputation variance, finite-sample pooling, BH ordering, and restricted spline tails. Full data-dependent verification uses the pipeline stages above.

## Files intended for GitHub

Upload the package contents, keeping `.gitignore`. Raw XPT files, participant-level pickle files, imputed values, and generated working outputs are excluded. Small aggregate reference JSON/CSV files are included to make numerical checks possible. The package does not include manuscript drafts, temporary Word files, local credentials, or font binaries.

The code repository is https://github.com/zhenxi-wang/nhanes-umfa-depression. Use this verified repository URL in the manuscript code-availability statement. No release DOI or software licence has been assigned. The authors should select an appropriate software licence before granting reuse permissions.
