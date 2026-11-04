"""Recalculate upper-tail counts, raw missingness, complete-case sensitivity and BH adjustments."""
from analysis_runtime import *
legacy=ROOT/'reference_results'/'legacy'
a=get_base();OUT.mkdir(parents=True,exist_ok=True)
counts=[]
for t in [2,4,8]:
 mask=a.umfa>t;counts.append(dict(threshold=t,n=int(mask.sum()),pct=100*float(mask.mean()),cases=int((a.loc[mask,'phq9']>=10).sum()),weighted_pct=100*float(a.loc[mask,'w'].sum()/a.w.sum()),n_equal=int((a.umfa==t).sum())))
dump('upper_tail_counts.json',counts)
cov={'Education':'edu','Income-to-poverty ratio':'pir','Marital status':'mar','Health insurance':'insured','Food security':'foodinsec','Nativity':'usborn','Smoking status':'smoke','Serum cotinine':'cotinine','Alcohol':'drink','Body mass index':'bmi','Sedentary time':'sed','Recreational physical activity':'mvpa_min','Sleep duration':'sleep','Chronic conditions':'ncond','Diabetes':'dm','eGFR':'egfr','Energy intake':'kcal','Dietary folate':'d_folate','Dietary folic acid':'d_fa','Any supplement use':'supp','Supplemental folic acid':'l_sfa','Fasting hours':'fasthr','MEC session':'session','hsCRP':'crp','Mean corpuscular volume':'mcv','Antidepressant use':'ad'}
dump('covariate_missingness.json',[dict(covariate=label,variable=v,missing_n=int(a[v].isna().sum()),missing_pct=100*float(a[v].isna().mean())) for label,v in cov.items()])
allcov=['edu','pir','foodinsec','usborn','mar','insured','smoke','l_cot','egfr','l_kcal','l_dfol','l_dfa','supp','l_sfa','fasthr','session','drink','bmi','sed','l_mvpa','sleep','ncond','dm']
ac=a.copy();ac['ins']=ac[allcov].notna().all(axis=1).astype(int);cc=pooled([fit(ac)],'UMFA')
simple=a.copy();simple['l_sfa']=np.log1p(simple.s_fa.fillna(0).clip(lower=0));simple_result=pooled([fit(simple)],'UMFA')
dump('missing_data_sensitivity.json',dict(strict_complete_case=cc,original_simple_imputation=simple_result))
screen_path=OUT.parent/'legacy'/'screen.csv'
if not screen_path.exists():screen_path=legacy/'screen.csv'
screen=pd.read_csv(screen_path);assert len(screen)==96
included=['log_fol_5mthf','log_fol_fa','log_fol_mefox','log_rbc_folate','logamh','log_fsh','log_lh','log_e2','log_prog','log_tst','log_shbg','log_dheas','log_andro','log_ohp17','log_e1','log_e1s','romberg_fail4','romberg_fail3','log_lsm','log_cap','log_bche','log_vitd']
small=screen.loc[screen.exposure.isin(included)].copy();assert len(small)==82
small['q_exact']=bh(small.p);screen['q_exact']=bh(screen.p)
small.to_csv(OUT/'screen_82_exact_bh.csv',index=False);screen.to_csv(OUT/'screen_96_exact_bh.csv',index=False)
dump('screen_bh_provenance.json',dict(source=str(screen_path),source_type='regenerated screen' if screen_path.parent==OUT.parent/'legacy' else 'recorded historical screen with unrounded P values',biomarker_family_n=82,all_screen_tests_n=96))
# Verify the historical 26-test family against its unrounded stored model outputs.
old={k:json.loads((legacy/f'out_a{k}.json').read_text()) for k in [1,3,4,5]};family=[]
for r in old[1]['outcomes']:
 if r['outcome']!='PHQ-9 >= 10':family.append((f"Outcome: {r['outcome']}",float(r['_UMFA'][3])))
for r in old[1]['quartiles']:
 ex=r['exposure'];family.extend([(f'Quartile Q4 vs Q1: {ex}',float(r['_q4'][3])),(f'Quartile joint test: {ex}',float(r['_p_joint'])),(f'Quartile trend (quartile median): {ex}',float(r['_p_trend'])),(f'Quartile trend (quartile rank): {ex}',float(old[4]['ordinal_trend'][ex]['p']))])
for r in old[3]['thresholds']:family.append((f"Threshold: {r['threshold']}",float(r['_p'])))
assert len(family)==26
oldq={r['Analysis']:float(r['BH q']) for r in old[5]['fdr']};qq=bh([p for _,p in family])
check=[dict(analysis=lab,p_exact=p,q_exact=float(q),matches_reported_rounding=f'{q:.3f}'==f'{oldq[lab]:.3f}') for (lab,p),q in zip(family,qq)]
assert all(x['matches_reported_rounding'] for x in check);dump('historical_s11_bh_check.json',check)
print('AUDIT COMPLETE: analytic sample, missingness, sensitivity models, and both screening families.')
