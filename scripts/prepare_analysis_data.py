"""Build both analytic frames from public XPT files; keep source files read-only."""
import argparse,hashlib,importlib.util,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source-data',type=Path,required=True);p.add_argument('--output-data',type=Path,required=True);a=p.parse_args()
 src=a.source_data.resolve();dst=a.output_data.resolve();dst.mkdir(parents=True,exist_ok=True)
 # pool.py constructs its cross-cycle drug dictionary at import time.
 import os
 os.environ['NHANES_DATA_DIR']=str(src)
 pool=load('nhanes_pool_builder',ROOT/'code'/'pool.py');pool.DATA=str(src);pooled=pool.build();pooled.to_pickle(dst/'pooled.pkl')
 prep=load('nhanes_screen_builder',ROOT/'code'/'prep.py');prep.DATA=str(src);wide=prep.build();wide.to_pickle(dst/'wide.pkl')
 m=(pooled.age>=18)&(pooled.preg!=1)&pooled.phq9.notna()&pooled.mthf.notna()&pooled.umfa.notna()&(pooled.w>0)
 assert int(m.sum())==24326,f'Unexpected analytic count: {int(m.sum())}'
 assert not pooled.loc[m,['cyc','SEQN']].duplicated().any()
 result={'pooled_rows':len(pooled),'discovery_frame_rows':len(wide),'analytic_n':int(m.sum()),'analytic_n_by_cycle':pooled.loc[m].groupby('cycle',sort=False).size().to_dict(),'source_xpt_directory':str(src),'output_directory':str(dst)}
 (dst/'preparation_summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
