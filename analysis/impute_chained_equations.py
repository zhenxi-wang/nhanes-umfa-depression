from analysis_runtime import *
import time
a=get_base()
tag=os.environ.get('NHANES_SUBSET','all')
if tag=='earlier':a=a.loc[a.cyc!='L'].copy().reset_index(drop=True)
elif tag!='all':a=a.loc[a.cyc==tag].copy().reset_index(drop=True)
OUT.mkdir(parents=True,exist_ok=True)
group=pd.factorize(pd.MultiIndex.from_frame(a[['strata','psu']]),sort=True)[0]
N=len(a);G=int(group.max()+1)
num=['pir','insured','foodinsec','usborn','l_cot','bmi','sed','l_mvpa','sleep','ncond','dm','egfr','l_kcal','l_dfol','l_dfa','supp','l_sfa','fasthr']
cat=['edu','mar','smoke','drink','session','race']
cols=[];blocks={};levels={};original={};scales={}
for v in num:
 y=a[v].to_numpy(float);original[v]=y.copy();mu=np.nanmean(y);sd=np.nanstd(y);sd=sd if sd>1e-10 else 1.
 scales[v]=(mu,sd);blocks[v]=[len(cols)];cols.append((y-mu)/sd)
for v in cat:
 lv=sorted(a[v].dropna().astype(str).unique());levels[v]=lv
 y=a[v].astype(str).map({z:i for i,z in enumerate(lv)}).to_numpy(float);original[v]=y.copy();blocks[v]=[]
 for k in range(1,len(lv)):
  blocks[v].append(len(cols));cols.append(np.where(np.isnan(y),np.nan,(y==k).astype(float)))
for v in ['age','female','l_umfa_z','l_mthf_z','l_mefox_z','dep']:
 y=a[v].to_numpy(float);cols.append((y-y.mean())/y.std())
aux=os.environ.get('NHANES_AUX','')
auxcols={'secondary':['dep15','dep5','si','phq9','som','cog'],'medication':['ad','anx','ap']}.get(aux,[])
for v in auxcols:
 y=a[v].to_numpy(float);assert np.isfinite(y).all(),v
 if y.std()>0:cols.append((y-y.mean())/y.std())
y=((a.age-a.age.mean())/10)**2;cols.append(((y-y.mean())/y.std()).to_numpy())
y=np.log(a.w.to_numpy());cols.append((y-y.mean())/y.std())
rawmat=np.column_stack(cols);variables=[v for v in num+cat if np.isnan(original[v]).any()]
variables.sort(key=lambda v:np.isnan(original[v]).sum());infos={}
for v in variables:
 mis=np.flatnonzero(np.isnan(original[v]));obs=np.flatnonzero(~np.isnan(original[v]));go=group[obs];counts=np.bincount(go,minlength=G)
 assert (counts>0).all();starts=np.r_[0,np.cumsum(counts)[:-1]];keep=np.array([j for j in range(rawmat.shape[1]) if j not in blocks[v]])
 infos[v]=(mis,obs,go,counts,starts,keep)
def put(mat,v,rows,values):
 if v in num:
  mu,sd=scales[v];mat[rows,blocks[v][0]]=(values-mu)/sd
 else:mat[np.ix_(rows,blocks[v])]=(values[:,None]==np.arange(1,len(levels[v]))[None,:]).astype(float)
def values(mat,v):
 if v in num:
  mu,sd=scales[v];return mat[:,blocks[v][0]]*sd+mu
 ar=np.column_stack([1-mat[:,blocks[v]].sum(axis=1),mat[:,blocks[v]]]);return ar.argmax(axis=1).astype(float)
def initialize(rng):
 mat=rawmat.copy()
 for v in variables:
  mis,obs,*_=infos[v];out=np.empty(len(mis))
  for cy in a.cyc.unique():
   mm=a.cyc.to_numpy()[mis]==cy;don=obs[a.cyc.to_numpy()[obs]==cy];out[mm]=rng.choice(original[v][don],int(mm.sum()))
  put(mat,v,mis,out)
 rows=np.flatnonzero(np.isnan(original['l_sfa'])&(values(mat,'supp')==0));put(mat,'l_sfa',rows,np.zeros(len(rows)))
 return mat
def update(mat,v,rng):
 mis,obs,go,counts,starts,keep=infos[v];X=mat[np.ix_(obs,keep)];xm=np.add.reduceat(X,starts,axis=0)/counts[:,None];X-=xm[go]
 S=X.T@X;ridge=1e-5*np.maximum(np.diag(S),1e-8);ch=linalg.cho_factor(S+np.diag(ridge),lower=True,check_finite=False);yo=original[v][obs].copy()
 if v in cat and len(levels[v])>2:
  Y=(yo[:,None]==np.arange(len(levels[v])-1)[None,:]).astype(float);ym=np.add.reduceat(Y,starts,axis=0)/counts[:,None];yc=Y-ym[go];xy=X.T@yc;overall=Y.mean(axis=0)
  H=ym.T@(ym*counts[:,None])+xy.T@linalg.cho_solve(ch,xy,check_finite=False)-len(obs)*np.outer(overall,overall)
  C=Y.T@Y-len(obs)*np.outer(overall,overall);eig,vec=linalg.eigh((H+H.T)/2,(C+C.T)/2,check_finite=False);quant=np.r_[vec[:,-1],0.]
  if quant[np.argmax(np.abs(quant))]<0:quant=-quant
  yo=quant[yo.astype(int)]
 yo=(yo-yo.mean())/(yo.std() or 1.);ym=np.add.reduceat(yo,starts)/counts;yc=yo-ym[go]
 bhat=linalg.cho_solve(ch,X.T@yc,check_finite=False);res=yc-X@bhat;sigma=np.sqrt((res@res)/rng.chisquare(len(obs)-len(keep)-G))
 bdraw=bhat+sigma*linalg.solve_triangular(ch[0].T,rng.normal(size=len(keep)),lower=False,check_finite=False)
 alphad=ym-xm@bdraw+sigma*rng.normal(size=G)/np.sqrt(counts);predobs=X@bhat+ym[go];predmis=mat[np.ix_(mis,keep)]@bdraw+alphad[group[mis]];donorobs=obs
 if v=='l_sfa':
  active=values(mat,'supp')[mis]!=0;put(mat,v,mis[~active],np.zeros((~active).sum()))
  if not active.any():return
  donors=a.supp.to_numpy()[obs]==1;donorobs=obs[donors];predobs=predobs[donors];mis=mis[active];predmis=predmis[active]
 idx=np.argsort(predobs);po=predobs[idx];donorobs=donorobs[idx];loc=np.searchsorted(po,predmis);cand=loc[:,None]+np.arange(-5,5)[None,:]
 valid=(cand>=0)&(cand<len(po));cand=np.clip(cand,0,len(po)-1);dist=np.abs(po[cand]-predmis[:,None]);dist[~valid]=np.inf
 near=np.argpartition(dist,4,axis=1)[:,:5];choice=near[np.arange(len(mis)),rng.integers(0,5,size=len(mis))];selected=donorobs[cand[np.arange(len(mis)),choice]];put(mat,v,mis,original[v][selected])
 if v=='supp':
  rows=np.flatnonzero(np.isnan(original['l_sfa'])&(values(mat,'supp')==0));put(mat,'l_sfa',rows,np.zeros(len(rows)))
saved={v:[] for v in variables};trace=[];t0=time.time()
for chain in range(3):
 rng=np.random.default_rng(20261104+chain);mat=initialize(rng)
 for sweep in range(1,71):
  for v in variables:update(mat,v,rng)
  if sweep%10==0:
   trace.append(dict(chain=chain+1,sweep=sweep,means={v:float(values(mat,v)[infos[v][0]].mean()) for v in variables}))
   print('CHAIN',chain+1,'SWEEP',sweep,'SECONDS',round(time.time()-t0,1),flush=True)
  if sweep>20 and sweep%5==0:
   for v in variables:
    vv=values(mat,v);known=~np.isnan(original[v]);assert np.allclose(vv[known],original[v][known]);saved[v].append(vv[infos[v][0]].copy())
a.to_pickle(OUT/'analytic_base.pkl');np.savez_compressed(OUT/'imputed_missing_values.npz',**{v:np.array(x) for v,x in saved.items()})
dump('imputation_metadata.json',dict(variables=variables,levels=levels,seed=20261104,chains=3,burn_in=20,interval=5,m=30,trace=trace))
print('COMPLETE_30_MICE',round(time.time()-t0,1),flush=True)
