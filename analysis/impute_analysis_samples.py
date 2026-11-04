import os,sys,subprocess,concurrent.futures
from pathlib import Path
out=Path(os.environ.get('NHANES_OUTPUT',str(Path(__file__).resolve().parents[1]/'outputs'/'revised')))
out.mkdir(parents=True,exist_ok=True)
def run(tag):
 env=os.environ.copy();env['NHANES_SUBSET']=tag;env['NHANES_OUTPUT']=str(out/('mi_'+tag))
 proc=subprocess.run([sys.executable,'-B',str(Path(__file__).parent/'impute_chained_equations.py')],env=env,capture_output=True,text=True)
 (out/('imputation_'+tag+'.log')).write_text(proc.stdout+'\n'+proc.stderr,encoding='utf-8')
 print('PARTITION',tag,'EXIT',proc.returncode,proc.stdout.splitlines()[-1:] if proc.returncode==0 else proc.stderr[-1200:],flush=True)
 if proc.returncode:raise RuntimeError(tag)
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,['earlier','L','G','H','I','P']))
