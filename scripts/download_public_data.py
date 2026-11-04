"""Download public CDC files in the input manifest, checking XPT headers and optional hashes."""
import argparse,csv,hashlib,urllib.request,warnings
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data-dir',type=Path,default=ROOT/'data');p.add_argument('--strict-hashes',action='store_true');p.add_argument('--list-only',action='store_true');a=p.parse_args()
 folder=a.data_dir.resolve();rows=list(csv.DictReader((ROOT/'config/raw_data_manifest.csv').open(encoding='utf-8')))
 if not a.list_only:folder.mkdir(parents=True,exist_ok=True)
 failures=[]
 for row in rows:
  path=folder/row['filename'];assert path.parent==folder
  urls=[row['download_url']]
  if row['filename']=='RXQ_DRUG.xpt':urls.extend(row['download_url'].replace('/2017/','/'+year+'/') for year in ['2019','2015'])
  if a.list_only:print(row['filename'],urls[0]);continue
  if path.exists():data=path.read_bytes()
  else:
   data=None
   for url in urls:
    try:
     req=urllib.request.Request(url,headers={'User-Agent':'NHANES-reproducibility/1.0'})
     with urllib.request.urlopen(req,timeout=180) as response:candidate=response.read()
     if not candidate.startswith(b'HEADER RECORD*******LIBRARY'):raise ValueError('Response is not a SAS XPORT file.')
     data=candidate;path.write_bytes(data);break
    except Exception as exc:last_error=str(exc)
   if data is None:failures.append((row['filename'],last_error));continue
  if not data.startswith(b'HEADER RECORD*******LIBRARY'):failures.append((row['filename'],'Invalid SAS XPORT header'));continue
  digest=hashlib.sha256(data).hexdigest()
  if digest!=row['sha256']:
   message=row['filename']+': current bytes differ from the recorded source release.'
   if a.strict_hashes:failures.append((row['filename'],message))
   else:warnings.warn(message)
  print('READY',row['filename'],len(data),flush=True)
 if failures:raise SystemExit('Input failures:\n'+'\n'.join(name+': '+error for name,error in failures))
if __name__=='__main__':main()
