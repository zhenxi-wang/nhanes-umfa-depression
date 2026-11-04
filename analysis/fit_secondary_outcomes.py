from analysis_runtime import *
import time,pickle
t0=time.time();rr={}
for i,z in enumerate(completed(OUT/'mi_secondary')):
 for y in ['dep15','dep5','si','phq9','som','cog']:
  rr.setdefault(y,[]).append(fit(z,y=y,linear=y in ['phq9','som','cog']))
 if (i+1)%10==0:print('SECONDARY',i+1,'SECONDS',round(time.time()-t0,1),flush=True)
res={y:{t:pooled(rec,t) for t in ['UMFA','5-MTHF','MeFox']} for y,rec in rr.items()}
core=json.loads((OUT/'core_results.json').read_text());fam=[]
labels={'dep15':'PHQ-9 >= 15','dep5':'PHQ-9 >= 5','si':'Suicidal ideation','phq9':'PHQ-9 total score','som':'Somatic subscale','cog':'Cognitive-affective subscale'}
for y in rr:fam.append(dict(analysis='Outcome: '+labels[y],p=res[y]['UMFA']['p']))
for label in ['5-MTHF','UMFA','Total serum folate','UMFA fraction']:
 fam.extend([dict(analysis='Quartile Q4 vs Q1: '+label,p=core['quart_'+label]['Q4']['p']),dict(analysis='Quartile joint test: '+label,p=core['quart_'+label]['joint']['p']),dict(analysis='Quartile trend (quartile median): '+label,p=core['trend_median_'+label]['p']),dict(analysis='Quartile trend (quartile rank): '+label,p=core['trend_rank_'+label]['p'])])
for label in ['p90','p95','ge2','ge5']:fam.append(dict(analysis='Threshold: '+label,p=core['threshold_'+label]['p']))
assert len(fam)==26
for r,q in zip(fam,bh([r['p'] for r in fam])):r['q']=float(q)
res['BH_family_26']=fam;dump('secondary_results.json',res)
with open(OUT/'secondary_model_records.pkl','wb') as f:pickle.dump(rr,f)
print('SECONDARY_UMFA',json.dumps({y:res[y]['UMFA'] for y in rr}),flush=True);print('BH_MIN',min(r['q'] for r in fam),flush=True)
