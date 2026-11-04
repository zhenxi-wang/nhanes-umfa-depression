"""Audit final-eligibility folate models using the recoverable screening specification.

The historical seven-row summary lacks its generating script. This diagnostic
does not claim to reproduce that summary or replace it in manuscript exports.
"""
from analysis_runtime import *
from screening_model import run

a = get_base()
eligible = a.loc[a.cycle == c.CYCLES[-1], 'SEQN']
d = pd.read_pickle(DATA / 'wide.pkl').reset_index(drop=True)
pooled = pd.read_pickle(DATA / 'pooled.pkl')
latest = pooled.loc[pooled.cycle == c.CYCLES[-1]].set_index('SEQN')
for variable in ['thf', 'fthf', 'methenyl']:
    d[variable] = d.SEQN.map(latest[variable])

variables = [
    ('fol_5mthf', '5-Methyltetrahydrofolate'),
    ('fol_fa', 'Unmetabolised folic acid'),
    ('fol_mefox', 'MeFox'),
    ('rbc_folate', 'Red-cell folate'),
    ('thf', 'Tetrahydrofolate'),
    ('fthf', '5-Formyltetrahydrofolate'),
    ('methenyl', '5,10-Methenyltetrahydrofolate'),
]
results = []
for variable, label in variables:
    positive = d.loc[d[variable] > 0, variable]
    transformed = 'log_' + variable
    d[transformed] = np.log(d[variable].clip(lower=positive.min()))
    mask = (d.SEQN.isin(eligible) & d[transformed].notna()).to_numpy()
    result = run(d, 'wph', transformed, label, mask=mask)
    if result is None:
        raise RuntimeError('Folate screening model failed: ' + label)
    results.append(result)

destination = OUT / 'screening_audit'
destination.mkdir(parents=True, exist_ok=True)
pd.DataFrame(results).to_csv(destination / 'reanalysis_exact.csv', index=False)
formatted = [dict(Biomarker=r['exposure'], n=int(r['n']),
                  OR=f"{r['OR']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})",
                  p=f"{r['p']:.4f}") for r in results]
(destination / 'reanalysis_formatted.json').write_text(json.dumps(formatted, indent=2), encoding='utf-8')
reference = json.loads((ROOT / 'reference_results/legacy/screen_folate_final.json').read_text())
comparison = [dict(biomarker=r['Biomarker'], recorded=r, reanalysed=n,
                   matches_recorded_rounding=(r == n)) for r, n in zip(reference, formatted)]
(destination / 'comparison_to_historical_summary.json').write_text(json.dumps(comparison, indent=2), encoding='utf-8')
print('AUDITED: seven final-eligibility models. Historical summary is unchanged.')
if not all(row['matches_recorded_rounding'] for row in comparison):
    print('NOTICE: the recoverable screening model does not exactly reproduce the historical seven-row summary; see docs/REPRODUCIBILITY_STATUS.md.')
