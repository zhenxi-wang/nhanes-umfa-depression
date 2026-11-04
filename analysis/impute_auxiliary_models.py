import os,sys,subprocess,concurrent.futures
from pathlib import Path
out=Path(os.environ.get('NHANES_OUTPUT',str(Path(__file__).resolve().parents[1]/'outputs'/'revised')))
out.mkdir(parents=True,exist_ok=True)
def run(spec):
 tag,aux=spec;env=os.environ.copy();env['NHANES_SUBSET']=tag;env['NHANES_AUX']=aux;env['NHANES_OUTPUT']=str(out/('mi_'+aux))
 p=subprocess.run([sys.executable,'-B',str(Path(__file__).parent/'impute_chained_equations.py')],env=env,capture_output=True,text=True)
 (out/('imputation_'+aux+'.log')).write_text(p.stdout+'\n'+p.stderr,encoding='utf-8');print('AUX',aux,'EXIT',p.returncode,p.stderr[-600:] if p.returncode else p.stdout.splitlines()[-1:],flush=True)
 if p.returncode:raise RuntimeError(aux)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[('all','secondary'),('earlier','medication')]))
