"""Run documented analysis stages from any working directory."""
import argparse,json,os,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
STAGES=['prepare','screen','impute-main','impute-samples','impute-aux','core','replication','secondary','margins','audit','tables','figures','validate']
SCRIPTS={'impute-main':'impute_chained_equations.py','impute-samples':'impute_analysis_samples.py','impute-aux':'impute_auxiliary_models.py','core':'fit_primary_models.py','replication':'fit_replication_and_sensitivity.py','secondary':'fit_secondary_outcomes.py','margins':'estimate_curves_and_validate.py','audit':'audit_counts_and_multiplicity.py','tables':'export_result_tables.py','figures':'make_main_figures.py','validate':'compare_reference_results.py','validate-primary':'validate_primary_model.py'}
def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--stage',choices=['all','download','legacy','validate-primary','audit-screen-final']+STAGES,default='all')
 p.add_argument('--data-dir',type=Path,default=ROOT/'data',help='Public XPT inputs and generated pooled.pkl/wide.pkl.')
 p.add_argument('--source-data',type=Path,help='Read XPT inputs here during preparation; write derived files only to --data-dir.')
 p.add_argument('--output-dir',type=Path,default=ROOT/'outputs'/'revised')
 p.add_argument('--font-family',default='Arial')
 p.add_argument('--font-dir',type=Path)
 p.add_argument('--dry-run',action='store_true',help='Print the stage commands without executing them.')
 p.add_argument('--download',action='store_true',help='Download manifest inputs before running the requested stages.')
 args=p.parse_args();data=args.data_dir.resolve();out=args.output_dir.resolve();legacy=out.parent/'legacy'
 env=os.environ.copy();env.update(NHANES_DATA_DIR=str(data),NHANES_OUTPUT=str(out),NHANES_FONT_FAMILY=args.font_family,NHANES_SUBSET='all',NHANES_AUX='',PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NHANES_LEGACY_FIGURES=str(legacy/'figures'))
 if args.font_dir:env['NHANES_FONT_DIR']=str(args.font_dir.resolve())
 def call(script,extra=(),cwd=None):
  command=[sys.executable,'-B',str(script),*map(str,extra)];print('RUN',subprocess.list2cmdline(command),flush=True)
  if not args.dry_run:subprocess.run(command,check=True,cwd=cwd or ROOT,env=env)
 selected=STAGES if args.stage=='all' else [args.stage]
 if args.download and 'download' not in selected:selected=['download']+selected
 if not args.dry_run:
  out.mkdir(parents=True,exist_ok=True);legacy.mkdir(parents=True,exist_ok=True);(legacy/'figures').mkdir(exist_ok=True)
 t0=time.time()
 for stage in selected:
  print('\nSTAGE',stage,flush=True)
  if stage=='download':call(ROOT/'scripts'/'download_public_data.py',['--data-dir',data])
  elif stage=='prepare':call(ROOT/'scripts'/'prepare_analysis_data.py',['--source-data',args.source_data.resolve() if args.source_data else data,'--output-data',data])
  elif stage=='screen':call(ROOT/'code'/'screen.py',cwd=legacy)
  elif stage=='audit-screen-final':call(ROOT/'analysis'/'audit_final_eligibility_screen.py')
  elif stage=='legacy':
   call(ROOT/'code'/'screen.py',cwd=legacy)
   if not args.dry_run:shutil.copyfile(ROOT/'reference_results/legacy/screen_folate_final.json',legacy/'screen_folate_final.json')
   for script in ['a0_desc.py','a1_main.py','a2_mech.py','a3_extra.py','a4_tables.py','a5_review.py','a6_review2.py','a7_species.py','flowfig.py','tables.py','figs.py','fig3_replication.py','mechfig.py','dagfig.py']:call(ROOT/'code'/script,cwd=legacy)
  else:call(ROOT/'analysis'/SCRIPTS[stage])
 if not args.dry_run:
  (out/'last_run.json').write_text(json.dumps({'stages':selected,'python':sys.version,'data_directory':str(data),'output_directory':str(out),'elapsed_seconds':time.time()-t0},indent=2),encoding='utf-8')
 print('COMPLETED',', '.join(selected),flush=True)
if __name__=='__main__':main()
