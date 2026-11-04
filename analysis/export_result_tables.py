"""Export English CSV/Markdown tables, preserving explicit missing-data provenance."""
from analysis_runtime import *
core=json.loads((OUT/'core_results.json').read_text());rep=json.loads((OUT/'replication_sensitivity_results.json').read_text());sec=json.loads((OUT/'secondary_results.json').read_text());cv=json.loads((OUT/'curves_margins_validation.json').read_text());missing=json.loads((OUT/'missing_data_sensitivity.json').read_text())
folder=OUT/'tables';folder.mkdir(exist_ok=True);sections=[]
def pval(p,d=3):return '<0.001' if p<.001 else f'{p:.{d}f}'
def ci(r,d=2):return f"{r['est']:.{d}f} ({r['lo']:.{d}f}-{r['hi']:.{d}f})"
def export(name,rows,note=''):
 frame=pd.DataFrame(rows);frame.to_csv(folder/(name+'.csv'),index=False)
 def clean(v):return str(v).replace('|','/').replace('\n',' ')
 md='| '+' | '.join(frame.columns)+' |\n| '+' | '.join(['---']*len(frame.columns))+' |\n'+'\n'.join('| '+' | '.join(clean(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None))
 sections.append('## '+name+'\n\n'+md+'\n\n'+note+'\n')
for kind,labels in [('Table_2A_mutually_adjusted',['5-MTHF','UMFA','MeFox']),('Table_2B_single_biomarker',['5-MTHF','Total serum folate','Red-cell folate','UMFA','UMFA fraction','MeFox'])]:
 rows=[]
 for lab in labels:
  row={'Biomarker':lab}
  for mo in [1,2,3]:
   r=core[f'joint_M{mo}'][lab] if kind.startswith('Table_2A') else core[f'single_{lab}_M{mo}'];row['n']=r['n'];row[f'Model {mo}, OR (95% CI)']=ci(r);row[f'Model {mo}, P']=pval(r['p'])
  rows.append(row)
 export(kind,rows,'All rows use 30 imputations. Biomarker exposures were not imputed.')
rows=[]
for lab in ['5-MTHF','UMFA','Total serum folate','UMFA fraction']:
 r=core['quart_'+lab];rows.append({'Biomarker':lab,'n':r['Q2']['n'],'Q1':'1.00 (reference)',**{q:ci(r[q]) for q in ['Q2','Q3','Q4']},'P joint':pval(r['joint']['p'],4)})
export('Table_3_single_biomarker_quartiles',rows,'MICE, 30 datasets. D1 joint tests with finite-sample correction.')
for name,terms in [('Formal_contrast',['5-MTHF','UMFA'])]:
 r=core['formal_contrast'];export(name,[{'Comparison':'5-MTHF versus UMFA','Ratio of ORs (95% CI)':ci(r),'P':pval(r['p']),'Log-OR difference':r['beta'],'Standard error':r['se']}])
export('Table_S13_replication',[{'Sample':lab,'n':r['n'],'OR (95% CI)':ci(r),'P':pval(r['p'])} for lab,r in [('Discovery',rep['L']),('Earlier-cycle replication',rep['earlier']),('Pooled',core['joint_M3']['UMFA'])]],'Imputations were performed separately within each analysis sample.')
export('Table_S13_individual_cycles',[{'Cycle':lab,'n':rep[tag]['n'],'OR (95% CI)':ci(rep[tag]),'P':pval(rep[tag]['p'])} for tag,lab in [('G','2011-2012'),('H','2013-2014'),('I','2015-2016'),('P','2017-March 2020'),('L','August 2021-August 2023')]])
export('Between_cycle_heterogeneity',[rep['heterogeneity']])
export('Table_S14A_standardisation',[{'Standardisation':lab,'n':r['n'],'OR (95% CI)':ci(r),'P':pval(r['p'])} for lab,r in [('Unweighted within-cycle SD, primary',core['joint_M3']['UMFA']),('Survey-weighted within-cycle SD',core['weighted_sd']['UMFA'])]])
a=get_base();des=c.design(a);targets={};qw=np.zeros(len(a));weighted_targets={90:np.zeros(len(a),bool),95:np.zeros(len(a),bool)}
for p in [90,95]:
 flag=np.zeros(len(a),bool)
 for cycle in c.CYCLES:
  ix=a.cycle==cycle;flag[ix]=a.loc[ix,'l_umfa']>np.percentile(a.loc[ix,'l_umfa'],p)
 targets['p'+str(p)]=flag
for p in [2,5]:targets['ge'+str(p)]=(a.umfa>=p).to_numpy()
for cycle in c.CYCLES:
 ix=a.cycle==cycle;cut=c.wquantile(a.loc[ix,'l_umfa'],a.loc[ix,'w'],[.25,.5,.75,.9,.95]);x=a.loc[ix,'l_umfa']
 qw[ix]=np.select([x<=cut[0],x<=cut[1],x<=cut[2],x>cut[2]],[1,2,3,4]);weighted_targets[90][ix]=x>cut[3];weighted_targets[95][ix]=x>cut[4]
export('Table_S14B_weighted_quartiles',[{'Quartile':q,'OR (95% CI)':ci(core['weighted_quart_mutual'][q])} for q in ['Q2','Q3','Q4']],f"Q1 is the reference; joint P={core['weighted_quart_mutual']['joint']['p']:.4f}.")
export('Table_S14B_weighted_thresholds',[{'Threshold':f'> cycle-specific weighted {p}th percentile','n exposed':int(flag.sum()),'OR (95% CI)':ci(core['threshold_weighted_p'+str(p)]),'P':pval(core['threshold_weighted_p'+str(p)]['p'])} for p,flag in weighted_targets.items()])
threshold_rows=[]
for key,flag in targets.items():
 r=core['threshold_'+key];we=a.loc[flag,'w'];wn=a.loc[~flag,'w'];threshold_rows.append({'Threshold':key,'n exposed':int(flag.sum()),'Weighted prevalence exposed, %':100*float(we@a.loc[flag,'dep']/we.sum()),'Weighted prevalence unexposed, %':100*float(wn@a.loc[~flag,'dep']/wn.sum()),'OR (95% CI)':ci(r),'P':pval(r['p'])})
export('Table_S14C_threshold_analyses',threshold_rows,'Adjusted estimates use MICE; prevalences describe observed exposure/outcome data.')
export('Table_S14D_upper_tail_counts',json.loads((OUT/'upper_tail_counts.json').read_text()),'Thresholds are strict greater-than comparisons using raw serum UMFA in nmol/L.')
r=core['quart_UMFA_mutual'];export('Table_S14E_mutually_adjusted_quartiles',[{'Quartile':'Q1','OR (95% CI)':'1.00 (reference)','P':'-'}]+[{'Quartile':q,'OR (95% CI)':ci(r[q]),'P':pval(r[q]['p'])} for q in ['Q2','Q3','Q4']],f"Three-species Model 3; joint P={r['joint']['p']:.4f}. Robustness analysis.")
export('Table_S11_exploratory_BH',sec['BH_family_26'],'Exactly 26 comparisons; BH adjustment uses unrounded P values.')
export('Table_S9_missingness',json.loads((OUT/'covariate_missingness.json').read_text()))
export('Missing_data_sensitivity',[{'Analysis':lab,'n':r['n'],'OR (95% CI)':ci(r),'P':pval(r['p'])} for lab,r in missing.items()])
export('Table_S4_MICE_sensitivities',[{'Analysis':lab,'n':r['n'],'OR (95% CI)':ci(r),'P':pval(r['p']),'Missing-data approach':'MICE, 30 datasets'} for lab,r in rep['sensitivity'].items()])
export('Table_S4_medication_analyses',[{'Analysis':lab,'n':r['n'],'OR (95% CI)':ci(r),'P':pval(r['p']),'Missing-data approach':'MICE restricted to medication-data cycles'} for lab,r in rep['medication'].items()])
legacy5=json.loads((ROOT/'reference_results/legacy/out_a5.json').read_text());legacy4=json.loads((ROOT/'reference_results/legacy/out_a4.json').read_text())
model4=next(r for r in legacy5['specifications'] if r['analysis'].startswith('Model 4'))
rows=[{'Specification':f'Model {m}','n':core[f'joint_M{m}']['UMFA']['n'],'OR (95% CI)':ci(core[f'joint_M{m}']['UMFA']),'P':pval(core[f'joint_M{m}']['UMFA']['p']),'Missing-data approach':'MICE, 30 datasets'} for m in [1,2,3]]
rows.append({'Specification':'DAG-guided Model 5','n':core['DAG']['UMFA']['n'],'OR (95% CI)':ci(core['DAG']['UMFA']),'P':pval(core['DAG']['UMFA']['p']),'Missing-data approach':'MICE, 30 datasets'})
rows.append({'Specification':model4['analysis'],'n':model4['n'],'OR (95% CI)':model4['OR'],'P':model4['p'],'Missing-data approach':'Historical simple imputation with indicators'})
export('Table_S4A_specifications',rows)
srh=next(r for r in legacy4['robust_extra'] if 'self-rated' in r['analysis'])
export('Table_S4_retained_legacy_rows',[{'Analysis':r['analysis'],'n':r['n'],'OR (95% CI)':r['OR'],'P':r['p'],'Missing-data approach':'Historical simple imputation with indicators'} for r in [model4,srh]],'These retained rows were not reclassified as new MICE analyses.')
rows=[]
for label,y in [('PHQ-9 >=15','dep15'),('PHQ-9 >=5','dep5'),('Suicidal ideation','si'),('PHQ-9 total','phq9'),('Somatic subscale','som'),('Cognitive-affective subscale','cog')]:
 row={'Outcome':label,'Estimate type':'OR' if y in ['dep15','dep5','si'] else 'Linear regression coefficient'}
 for term in ['UMFA','5-MTHF','MeFox']:row[term+' (95% CI)']=ci(sec[y][term],2 if row['Estimate type']=='OR' else 3);row[term+' P']=pval(sec[y][term]['p'])
 rows.append(row)
export('Table_S22_secondary_outcomes',rows,'Secondary outcomes were included as auxiliary predictors in the secondary-outcome imputations.')
export('Figure_2_spline_points',[{'UMFA nmol/L':r['umfa'],'OR (95% CI)':ci(r),'P':pval(r['p'])} for r in cv['at_concentrations']])
export('Figure_2_adjusted_prevalence',[{k:r[k] for k in ['percentile','umfa','prev','lo','hi']} for r in cv['prevalence']])
export('Figure_2_absolute_difference',[cv['risk_difference']])
export('E_values',[cv['evalue']])
export('Table_S10_MICE_pooling_validation',[{'Quantity':k,'Value':v} for k,v in cv['validation'].items() if isinstance(v,(float,int))])
summary={'single_biomarker_quartile_joint_P':core['quart_UMFA']['joint']['p'],'mutually_adjusted_quartile_joint_P':core['quart_UMFA_mutual']['joint']['p'],'spline_nonlinearity_P':core['spline']['nonlinearity']['p'],'formal_contrast_P':core['formal_contrast']['p']}
export('Exposure_specification_checks',[summary])
(folder/'revised_tables.md').write_text('# Revised analysis tables\n\n'+'\n'.join(sections),encoding='utf-8')
print('EXPORTED',len(sections),'tables to',folder)
