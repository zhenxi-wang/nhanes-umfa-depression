"""Compare regenerated summaries with the archived revised-analysis results."""
from analysis_runtime import *
files=['core_results.json','replication_sensitivity_results.json','secondary_results.json','curves_margins_validation.json','upper_tail_counts.json']
checks=[]
def walk(got,ref,path):
 if isinstance(ref,dict):
  for key,value in ref.items():
   if key=='independent_optimisation':continue
   if key not in got:raise AssertionError('Missing result: '+path+'/'+key)
   walk(got[key],value,path+'/'+key)
 elif isinstance(ref,list):
  assert len(got)==len(ref),path
  for i,(g,r) in enumerate(zip(got,ref)):walk(g,r,path+'/'+str(i))
 elif isinstance(ref,(int,float)) and not isinstance(ref,bool):
  assert np.isclose(got,ref,rtol=2e-5,atol=1e-8),f'{path}: regenerated={got}; reference={ref}'
 else:assert got==ref,path
for filename in files:
 path=OUT/filename
 if not path.exists():raise FileNotFoundError('Run the required pipeline stages first: '+str(path))
 walk(json.loads(path.read_text()),json.loads((ROOT/'reference_results/revised'/filename).read_text()),filename);checks.append(filename)
dump('reference_comparison.json',dict(status='passed',files=checks,relative_tolerance=2e-5,absolute_tolerance=1e-8));print('REFERENCE COMPARISON PASSED:',', '.join(checks))
