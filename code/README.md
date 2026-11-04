# Shared helpers and historical analyses

Use `reproduce.py` at the repository root to control working directories and output locations. Directly running historical scripts from an arbitrary working directory can place their outputs there.

`pool.py`, `prep.py`, `common.py`, `survey.py`, and `screening_model.py` support data preparation and estimation. Other scripts in this directory preserve the original screening, descriptive, simple-imputation and exploratory analyses. The current MICE analyses and current main figures are in `analysis/`.

`a5_review.py` contains the historical 10-imputation sensitivity calculation; it is not the current 30-imputation primary procedure. Original figure numbering in `figs.py` predates the current main-figure order. `analysis/make_main_figures.py` is the source for the three revised main figures.

See `CODE_INVENTORY.md` and `docs/ANALYSIS_MAP.md` for individual file roles and manuscript locations.
