from analysis_runtime import *
import time,pickle
t0=time.time();allr={};res={}
for tag in ['earlier','L','G','H','I','P']:
 recs=[fit(z) for z in completed(OUT/('mi_'+tag))];allr[tag]=recs;res[tag]=pooled(recs,'UMFA')
 print('REPLICATION',tag,json.dumps(res[tag]),flush=True)
bs=np.array([res[tag]['beta'] for tag in ['G','H','I','P','L']]);vs=np.array([res[tag]['se']**2 for tag in ['G','H','I','P','L']]);wt=1/vs;mu=np.sum(wt*bs)/sum(wt);Q=float(np.sum(wt*(bs-mu)**2))
res['heterogeneity']=dict(Q=Q,df=4,p=float(stats.chi2.sf(Q,4)),I2=max(0.,100*(Q-4)/Q))
med={};mbase=pd.read_pickle(OUT/'mi_medication'/'analytic_base.pkl')
for i,z in enumerate(completed(OUT/'mi_medication')):
 def add(key,**kw):med.setdefault(key,[]).append(fit(z,**kw))
 add('med_cycles_baseline');add('antidepressant_adjusted',extra=[z.ad.to_numpy()],extra_names=['Antidepressant'])
 add('psychotropic_adjusted',extra=[z[v].to_numpy() for v in ['ad','anx','ap']],extra_names=['Antidepressant','Anxiolytic','Antipsychotic'])
 add('antidepressant_users',mask=mbase.ad==1);add('antidepressant_nonusers',mask=mbase.ad==0)
 add('no_psychotropic',mask=(mbase.ad==0)&(mbase.anx==0)&(mbase.ap==0))
 add('no_antidepressant_no_supplement',mask=(mbase.ad==0)&(mbase.supp==0))
 add('no_antidepressant_no_FA_supplement',mask=(mbase.ad==0)&(mbase.fa_supp==0))
 add('antidepressant_outcome',y='ad')
 add('antidepressant_interaction',extra=[z.ad.to_numpy(),z.ad.to_numpy()*z.l_umfa_z.to_numpy()],extra_names=['Antidepressant','UMFA x antidepressant'])
 if (i+1)%10==0:print('MEDICATION',i+1,'SECONDS',round(time.time()-t0,1),flush=True)
res['medication']={key:pooled(rr,'UMFA x antidepressant' if key=='antidepressant_interaction' else 'UMFA') for key,rr in med.items()}
sens={};a=pd.read_pickle(OUT/'analytic_base.pkl')
for i,z in enumerate(completed()):
 def add(key,**kw):sens.setdefault(key,[]).append(fit(z,**kw))
 add('excluding_2011_2012',mask=a.cyc!='G')
 add('no_dietary_supplement_covariates',drop_terms=['log energy intake','log dietary folate','log dietary folic acid','Any supplement use','log supplemental folic acid'])
 add('FA_supp_users',mask=a.fa_supp==1);add('FA_supp_nonusers',mask=a.fa_supp==0);add('any_supp_nonusers',mask=a.supp==0)
 for hr in [8,10,12]:add('fasting_'+str(hr),mask=a.fasthr>=hr)
 add('morning',mask=a.session==1);add('excluding_egfr_lt60',mask=a.egfr>=60);add('excluding_mcv_ge100',mask=a.mcv<100);add('excluding_bmi_ge40',mask=a.bmi<40)
 add('excluding_umfa_top1pct',mask=a.l_umfa<=a.l_umfa.quantile(.99));add('add_RBC',mask=a.rbcfol.notna(),extra=[z.l_rbcfol_z.to_numpy()],extra_names=['RBC folate'])
 add('age_ge20',mask=a.age>=20);add('NH_White',mask=a.race=='NH White')
 if (i+1)%10==0:print('SENSITIVITY',i+1,'SECONDS',round(time.time()-t0,1),flush=True)
res['sensitivity']={key:pooled(rr,'UMFA') for key,rr in sens.items()}
dump('replication_sensitivity_results.json',res)
with open(OUT/'replication_sensitivity_records.pkl','wb') as f:pickle.dump(dict(replication=allr,medication=med,sensitivity=sens),f)
print('MED_FINAL',json.dumps(res['medication']),flush=True);print('HETEROGENEITY',json.dumps(res['heterogeneity']),flush=True)
