"""Re-fit all 30 completed primary models and check the key manuscript results."""
from analysis_runtime import *
recs=[fit(z) for z in completed()]
result={t:pooled(recs,t) for t in ['UMFA','5-MTHF','MeFox']};result['formal_contrast']=contrast(recs,['5-MTHF','UMFA'],[1,-1])
reference=json.loads((ROOT/'reference_results/revised/core_results.json').read_text())
for term in ['UMFA','5-MTHF','MeFox']:
 for key in ['est','lo','hi','p','beta','se']:
  assert np.isclose(result[term][key],reference['joint_M3'][term][key],rtol=1e-6,atol=1e-9),(term,key,result[term][key])
for key in ['est','lo','hi','p']:assert np.isclose(result['formal_contrast'][key],reference['formal_contrast'][key],rtol=1e-6,atol=1e-9)
assert np.isclose(result['formal_contrast']['est'],result['5-MTHF']['est']/result['UMFA']['est'])
dump('primary_validation.json',dict(status='passed',results=result));print(json.dumps(result,indent=2))
