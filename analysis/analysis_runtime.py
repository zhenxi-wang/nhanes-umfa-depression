import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
import sys,json,warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg,stats
warnings.filterwarnings('ignore',category=pd.errors.PerformanceWarning)
ROOT=Path(__file__).resolve().parents[1]
DATA=Path(os.environ.get('NHANES_DATA_DIR',str(ROOT/'data'))).resolve()
OUT=Path(os.environ.get('NHANES_OUTPUT',str(ROOT/'outputs'/'revised'))).resolve()
sys.path.insert(0,str(ROOT/'code'))
import survey as s
import common as c
c.DATA=str(DATA)

def get_base():
 d=c.get();a=d.loc[d.ins==1].copy().reset_index(drop=True)
 group=pd.factorize(pd.MultiIndex.from_frame(a[['strata','psu']]),sort=True)[0]
 order=np.argsort(group,kind='stable');a=a.iloc[order].reset_index(drop=True)
 a.loc[a.supp.isna()&a.s_fa.isna(),'l_sfa']=np.nan
 return a

def completed(folder=None):
 folder=OUT if folder is None else Path(folder)
 a=pd.read_pickle(folder/'analytic_base.pkl')
 vals=np.load(folder/'imputed_missing_values.npz')
 meta=json.loads((folder/'imputation_metadata.json').read_text())
 for i in range(30):
  z=a.copy()
  for v in meta['variables']:
   ix=np.flatnonzero(a[v].isna());vv=vals[v][i]
   if v in meta['levels']:
    z.loc[ix,v]=[meta['levels'][v][int(q)] for q in vv]
   else:z.loc[ix,v]=vv
  yield z

def fit(z,level=3,exposures=('l_umfa_z','l_mthf_z','l_mefox_z'),names=('UMFA','5-MTHF','MeFox'),mask=None,y='dep',linear=False,extra=None,extra_names=None,drop_terms=()):
 X,nm=c.X_model(z,level)
 if drop_terms:
  keep=[i for i,t in enumerate(nm) if t not in drop_terms];X=X[:,keep];nm=[nm[i] for i in keep]
 if exposures:X,nm=c.add(X,nm,[z[v].to_numpy(float) for v in exposures],list(names))
 if extra is not None:X,nm=c.add(X,nm,list(extra),list(extra_names))
 mk=(z.ins==1).to_numpy() if mask is None else (z.ins==1).to_numpy()&np.asarray(mask)
 rr,V,n,df=(s.svylin if linear else s.svylogit)(c.design(z),z[y].to_numpy(float),X,mask=mk,names=nm)
 included_cycles=z.loc[mk,'cycle'].unique()
 zd=z.loc[z.cycle.isin(included_cycles)]
 df=zd[['strata','psu']].drop_duplicates().shape[0]-zd.strata.nunique()
 return dict(b=rr.beta.to_numpy(),V=V,n=n,df=int(df),names=rr.term.tolist(),linear=linear)

def scalar(q,u,dfcom=87,exponentiate=True):
 q=np.asarray(q);u=np.asarray(u);m=len(q);b=q.mean();U=u.mean();B=q.var(ddof=1) if m>1 else 0.
 T=U+(1+1/m)*B;se=np.sqrt(max(T,0));lam=(1+1/m)*B/T if T>0 else 0.
 old=(m-1)/max(lam,1e-12)**2 if m>1 else np.inf
 obs=(dfcom+1)/(dfcom+3)*dfcom*(1-lam)
 df=1/(1/old+1/obs) if m>1 else dfcom
 crit=stats.t.ppf(.975,df);p=2*stats.t.sf(abs(b/se),df) if se else 1.
 trans=np.exp if exponentiate else lambda x:x
 return dict(est=float(trans(b)),lo=float(trans(b-crit*se)),hi=float(trans(b+crit*se)),p=float(p),beta=float(b),se=float(se),within=float(U),between=float(B),df=float(df),m=m)

def pooled(recs,term,dfcom=None):
 q=[];u=[]
 for r in recs:
  j=r['names'].index(term);q.append(r['b'][j]);u.append(r['V'][j,j])
 out=scalar(q,u,min(r['df'] for r in recs) if dfcom is None else dfcom,not recs[0]['linear'])
 out['n']=recs[0]['n'];return out

def contrast(recs,terms,weights,dfcom=None,exponentiate=True):
 q=[];u=[]
 for r in recs:
  v=np.zeros(len(r['b']))
  for t,w in zip(terms,weights):v[r['names'].index(t)]=w
  q.append(v@r['b']);u.append(v@r['V']@v)
 return scalar(q,u,min(r['df'] for r in recs) if dfcom is None else dfcom,exponentiate)

def joint(recs,terms,dfcom=None):
 Q=[];U=[]
 for r in recs:
  ix=[r['names'].index(t) for t in terms];Q.append(r['b'][ix]);U.append(r['V'][np.ix_(ix,ix)])
 Q=np.array(Q);U=np.array(U);m,k=Q.shape;qb=Q.mean(axis=0);ub=U.mean(axis=0);B=np.cov(Q,rowvar=False,ddof=1)
 B=np.atleast_2d(B);r=(1+1/m)*np.trace(B@np.linalg.inv(ub))/k
 F=float(qb@np.linalg.solve((1+r)*ub,qb)/k)
 t=k*(m-1);dc=min(x['df'] for x in recs) if dfcom is None else dfcom
 aa=r*t/(t-2);vs=(dc+1)/(dc+3)*dc;c0=1/(t-4);c1=vs-2*(1+aa);c2=vs-4*(1+aa)
 zz=1/c2+c0*(aa*aa*c1/((1+aa)**2*c2))+c0*(8*aa*aa*c1/((1+aa)*c2*c2)+4*aa*aa/((1+aa)*c2))+c0*(4*aa*aa/(c2*c1)+16*aa*aa*c1/c2**3)+c0*(8*aa*aa/c2**2)
 df=4+1/zz
 return dict(F=F,df1=k,df2=float(df),p=float(stats.f.sf(F,k,df)),r=float(r))

def rcs(x,knots):
 t=np.asarray(knots);out=[]
 for j in range(len(t)-2):
  cu=lambda kn:np.maximum(np.asarray(x)-kn,0)**3
  out.append((cu(t[j])-cu(t[-2])*(t[-1]-t[j])/(t[-1]-t[-2])+cu(t[-1])*(t[-2]-t[j])/(t[-1]-t[-2]))/(t[-1]-t[0])**2)
 return out

def bh(p):
 p=np.array(p,float);o=np.argsort(p);q=np.empty(len(p));q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.clip(q,0,1)

def dump(name,data):
 def default(x):
  if isinstance(x,np.ndarray):return x.tolist()
  if isinstance(x,np.generic):return x.item()
  raise TypeError(type(x).__name__)
 (OUT/name).write_text(json.dumps(data,ensure_ascii=False,indent=2,default=default),encoding='utf-8')
