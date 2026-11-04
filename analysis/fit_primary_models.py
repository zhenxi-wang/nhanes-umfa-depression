from analysis_runtime import *
import time,pickle
base=pd.read_pickle(OUT/'analytic_base.pkl');t0=time.time();records={}
def savefit(key,z,**kw):records.setdefault(key,[]).append(fit(z,**kw))
biomarkers={'5-MTHF':'l_mthf_z','Total serum folate':'l_tfol_z','Red-cell folate':'l_rbcfol_z','UMFA':'l_umfa_z','UMFA fraction':'l_pct_umfa_z','MeFox':'l_mefox_z'}
quart={'5-MTHF':'l_mthf','UMFA':'l_umfa','Total serum folate':'l_tfol','UMFA fraction':'l_pct_umfa'}
kn=np.percentile(base.l_umfa,[5,35,65,95]);ref=float(np.percentile(base.l_umfa,10))
for i,z in enumerate(completed()):
 for mo in [1,2,3]:
  savefit(f'joint_M{mo}',z,level=mo)
  for label,col in biomarkers.items():savefit(f'single_{label}_M{mo}',z,level=mo,exposures=(col,),names=('Exposure',))
 savefit('DAG',z,level=5)
 for label,lv in quart.items():
  q=z[lv+'_q'];ex=[(q==k).astype(float).to_numpy() for k in [2,3,4]]
  savefit('quart_'+label,z,exposures=(),names=(),extra=ex,extra_names=['Q2','Q3','Q4'],mask=q.notna())
  savefit('trend_median_'+label,z,exposures=(),names=(),extra=[np.log(z[lv+'_qmed'].astype(float)).fillna(0).to_numpy()],extra_names=['Trend'],mask=q.notna())
  savefit('trend_rank_'+label,z,exposures=(),names=(),extra=[q.fillna(0).to_numpy()],extra_names=['Trend'],mask=q.notna())
 q=z.l_umfa_q
 savefit('quart_UMFA_mutual',z,exposures=('l_mthf_z','l_mefox_z'),names=('5-MTHF','MeFox'),extra=[(q==k).astype(float).to_numpy() for k in [2,3,4]],extra_names=['Q2','Q3','Q4'])
 x=z.l_umfa.to_numpy()
 savefit('spline',z,exposures=('l_mthf_z','l_mefox_z'),names=('5-MTHF','MeFox'),extra=[x-kn[0]]+rcs(x,kn),extra_names=['umfa_lin','umfa_s1','umfa_s2'])
 # Linear-log model for marginal prevalence, matching the original panel b.
 savefit('prevalence_loglinear',z,exposures=('l_umfa','l_mthf_z','l_mefox_z'),names=('UMFA raw log','5-MTHF','MeFox'))
 for v in ['l_umfa','l_mthf','l_mefox']:
  wz=np.zeros(len(z))
  for cy in c.CYCLES:
   sel=z.cycle==cy;mu,sd=c.wmean_sd(z.loc[sel,v],z.loc[sel,'w']);wz[sel]=(z.loc[sel,v]-mu)/sd
  z[v+'_wz']=wz
 savefit('weighted_sd',z,exposures=('l_umfa_wz','l_mthf_wz','l_mefox_wz'))
 for p in [90,95]:
  exp=np.zeros(len(z))
  for cy in c.CYCLES:
   ix=z.cycle==cy;cut=np.percentile(z.loc[ix,'l_umfa'],p);exp[ix]=z.loc[ix,'l_umfa']>cut
  savefit('threshold_p'+str(p),z,exposures=('l_mthf_z','l_mefox_z'),names=('5-MTHF','MeFox'),extra=[exp],extra_names=['Threshold'])
 for p in [2,5]:savefit('threshold_ge'+str(p),z,exposures=('l_mthf_z','l_mefox_z'),names=('5-MTHF','MeFox'),extra=[(z.umfa>=p).astype(float).to_numpy()],extra_names=['Threshold'])
 q_w=np.zeros(len(z));p90w=np.zeros(len(z));p95w=np.zeros(len(z))
 for cy in c.CYCLES:
  ix=z.cycle==cy;cut=c.wquantile(z.loc[ix,'l_umfa'],z.loc[ix,'w'],[.25,.5,.75,.90,.95]);xx=z.loc[ix,'l_umfa']
  q_w[ix]=np.select([xx<=cut[0],xx<=cut[1],xx<=cut[2],xx>cut[2]],[1,2,3,4]);p90w[ix]=xx>cut[3];p95w[ix]=xx>cut[4]
 savefit('weighted_quart_mutual',z,exposures=('l_mthf_z','l_mefox_z'),names=('5-MTHF','MeFox'),extra=[(q_w==k).astype(float) for k in [2,3,4]],extra_names=['Q2','Q3','Q4'])
 for p,xx in [(90,p90w),(95,p95w)]:savefit('threshold_weighted_p'+str(p),z,exposures=('l_mthf_z','l_mefox_z'),names=('5-MTHF','MeFox'),extra=[xx],extra_names=['Threshold'])
 print('CORE IMPUTATION',i+1,'MODELS',len(records),'SECONDS',round(time.time()-t0,1),flush=True)
with open(OUT/'core_model_records.pkl','wb') as f:pickle.dump(records,f)
res={}
for key,recs in records.items():
 if key.startswith('joint') or key in ['DAG','weighted_sd']:res[key]={t:pooled(recs,t) for t in ['UMFA','5-MTHF','MeFox']}
 elif 'quart' in key:
  res[key]={t:pooled(recs,t) for t in ['Q2','Q3','Q4']};res[key]['joint']=joint(recs,['Q2','Q3','Q4'])
 elif key=='spline':res[key]={'nonlinearity':joint(recs,['umfa_s1','umfa_s2']),'overall':joint(recs,['umfa_lin','umfa_s1','umfa_s2'])}
 elif key=='prevalence_loglinear':res[key]=pooled(recs,'UMFA raw log')
 elif key.startswith('single'):res[key]=pooled(recs,'Exposure')
 elif key.startswith('trend'):res[key]=pooled(recs,'Trend')
 else:res[key]=pooled(recs,'Threshold')
res['formal_contrast']=contrast(records['joint_M3'],['5-MTHF','UMFA'],[1,-1])
res['spline_knots_log']=kn;res['spline_reference']=np.exp(ref)
dump('core_results.json',res)
print('PRIMARY',json.dumps(res['joint_M3']),flush=True);print('CONTRAST',json.dumps(res['formal_contrast']),flush=True)
print('SPLINE',json.dumps(res['spline']),flush=True);print('QUART_MUTUAL',json.dumps(res['quart_UMFA_mutual']),flush=True)
