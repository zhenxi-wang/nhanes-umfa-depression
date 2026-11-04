from analysis_runtime import *
import pickle,time
from scipy.special import expit,logit
from scipy.optimize import minimize
records=pickle.load(open(OUT/'core_model_records.pkl','rb'));base=pd.read_pickle(OUT/'analytic_base.pkl');core=json.loads((OUT/'core_results.json').read_text())
kn=np.array(core['spline_knots_log']);ref=np.log(core['spline_reference'])
def bas(g):return np.array([g-kn[0]]+[x[0] for x in rcs(np.array([g]),kn)])
def curvepoint(g):
 dif=bas(g)-bas(ref);q=[];u=[]
 for rec in records['spline']:
  ix=[rec['names'].index(t) for t in ['umfa_lin','umfa_s1','umfa_s2']];q.append(dif@rec['b'][ix]);u.append(dif@rec['V'][np.ix_(ix,ix)]@dif)
 return dict(umfa=float(np.exp(g)),**scalar(q,u))
grid=np.unique(np.r_[np.linspace(np.percentile(base.l_umfa,2),np.percentile(base.l_umfa,99),180),ref,np.log([.5,1,2,4,8])])
curve=[curvepoint(g) for g in grid]
pcts=[10,25,50,75,90,95,99];qs=np.percentile(base.l_umfa,pcts);marg={pc:[] for pc in pcts};diff=[]
t0=__import__('time').time()
for i,z in enumerate(completed()):
 rec=records['spline'][i];X,nm=c.X_model(z,3);X,nm=c.add(X,nm,[z.l_mthf_z.to_numpy(),z.l_mefox_z.to_numpy()],['5-MTHF','MeFox']);X,nm=c.add(X,nm,[z.l_umfa.to_numpy()-kn[0]]+rcs(z.l_umfa.to_numpy(),kn),['umfa_lin','umfa_s1','umfa_s2'])
 assert nm==rec['names'];ix=[nm.index(t) for t in ['umfa_lin','umfa_s1','umfa_s2']]
 b=rec['b'];V=rec['V'];w=z.w.to_numpy();w=w/w.sum();eta0=X@b-X[:,ix]@b[ix];est=[];grads=[]
 for pc,q in zip(pcts,qs):
  vv=bas(q);mu=expit(eta0+vv@b[ix]);pr=float(w@mu);v=w*mu*(1-mu);g=X.T@v;g[ix]=vv*v.sum();var=float(g@V@g)
  marg[pc].append((logit(pr),var/(pr*(1-pr))**2));est.append(pr);grads.append(g)
 gd=grads[-1]-grads[0];diff.append((est[-1]-est[0],float(gd@V@gd)))
 if i==0:
  # Independent numerical optimisation of a completed-data likelihood.
  wr=z.w.to_numpy()/z.w.mean();yy=z.dep.to_numpy()
  def f(bb):
   eta=X@bb;return np.sum(wr*(np.logaddexp(0,eta)-yy*eta))
  def jac(bb):return X.T@(wr*(expit(X@bb)-yy))
  opt=minimize(f,np.zeros(X.shape[1]),jac=jac,method='BFGS',options={'gtol':1e-7,'maxiter':3000})
  js=[nm.index(t) for t in ['umfa_lin','umfa_s1','umfa_s2','5-MTHF','MeFox']]
  err=float(np.max(np.abs(opt.x[js]-b[js])));assert err<1e-5,err
  optimcheck=dict(max_absolute_substantive_coefficient_difference=err,gradient_max=float(np.abs(jac(opt.x)).max()),completed_model='first MICE spline model')
prevalence=[]
for pc,q in zip(pcts,qs):
 qq,uu=np.array(marg[pc]).T;po=scalar(qq,uu,exponentiate=False)
 prevalence.append(dict(percentile=pc,umfa=float(np.exp(q)),prev=100*float(expit(po['est'])),lo=100*float(expit(po['lo'])),hi=100*float(expit(po['hi'])),logit_pooled=po))
qq,uu=np.array(diff).T;rd=scalar(qq,uu,exponentiate=False)
rd={k:(100*v if k in ['est','lo','hi','beta','se'] else v) for k,v in rd.items()}
r=core['joint_M3']['UMFA'];evalue=lambda x:np.sqrt(x)+np.sqrt(np.sqrt(x)*(np.sqrt(x)-1))
ev=dict(point=float(evalue(r['est'])),lower_confidence_limit=float(evalue(r['lo'])))
# Independent matrix construction for scalar pooling and formal contrast.
Q=np.array([x['b'][-3:] for x in records['joint_M3']]);names=records['joint_M3'][0]['names'][-3:]
UU=np.array([x['V'][-3:,-3:] for x in records['joint_M3']]);qb=Q.mean(axis=0);Ub=UU.mean(axis=0);B=np.cov(Q,rowvar=False);T=Ub+(31/30)*B
ju=names.index('UMFA');jm=names.index('5-MTHF');vec=np.zeros(3);vec[jm]=1;vec[ju]=-1
assert np.isclose(qb[ju],r['beta']);assert np.isclose(T[ju,ju],r['se']**2)
con=core['formal_contrast'];assert np.isclose(vec@qb,con['beta']);assert np.isclose(vec@T@vec,con['se']**2)
validation=dict(primary_beta=float(qb[ju]),within_variance=float(Ub[ju,ju]),between_variance=float(B[ju,ju]),total_variance=float(T[ju,ju]),df=r['df'],p=r['p'],formal_contrast_beta=float(vec@qb),formal_contrast_total_variance=float(vec@T@vec),pooled_coefficient_covariance_matrix=T.tolist(),terms=names,independent_optimisation=optimcheck)
res=dict(curve=curve,at_concentrations=[curvepoint(np.log(g)) for g in [.5,1,2,4,8]],prevalence=prevalence,risk_difference=rd,evalue=ev,validation=validation,knots_nmol=np.exp(kn).tolist(),reference=np.exp(ref),nonlinearity=core['spline']['nonlinearity'])
dump('curves_margins_validation.json',res)
print('SPLINE_POINTS',json.dumps(res['at_concentrations']),flush=True);print('PREVALENCE',json.dumps(prevalence),flush=True);print('RISK_DIFF',json.dumps(rd),'EVALUES',json.dumps(ev),flush=True);print('VALIDATION',json.dumps(validation),flush=True)
