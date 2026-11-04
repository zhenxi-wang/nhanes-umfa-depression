"""Additional analyses requested at review: DAG-minimal model, cycle/mode sensitivity,
folic-acid-specific supplement strata, stricter fasting, multiple imputation, FDR control,
covariate missingness audit and a survey-estimation validation check."""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from common import get, design, X_model, add, cont, dum, fmt, pfmt, wquantile, CYCLES
from survey import SvyDesign, svylogit

d = get()
des = design(d)
ins = (d.ins == 1).to_numpy()
OUT = {}
pd.set_option("display.width", 250)
EX = [("UMFA", "l_umfa_z"), ("5-MTHF", "l_mthf_z")]


def fit(y, cols, names, level, mask=None, dd=None, dsg=None):
    dd = d if dd is None else dd
    dsg = des if dsg is None else dsg
    X, nm = X_model(dd, level)
    X, nm = add(X, nm, cols, names)
    mk = ins if mask is None else (ins & mask)
    return svylogit(dsg, dd[y].to_numpy(float), X, mask=mk, names=nm)


def umfa_row(label, mask=None, level=3, extra=(), extranm=(), dd=None, dsg=None):
    dd = d if dd is None else dd
    cols = [dd.l_umfa_z.to_numpy(), dd.l_mthf_z.to_numpy(), dd.l_mefox_z.to_numpy()] + list(extra)
    nms = ["UMFA", "5-MTHF", "MeFox"] + list(extranm)
    res, V, n, df = fit("dep", cols, nms, level, mask=mask, dd=dd, dsg=dsg)
    r = res[res.term == "UMFA"].iloc[0]
    return dict(analysis=label, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                _or=float(r.OR), _lo=float(r.lo), _hi=float(r.hi), _p=float(r.p))


# ------------------------------------------------ 1. covariate missingness audit
print("=" * 100)
print("COVARIATE MISSINGNESS IN THE ANALYTIC SAMPLE (n = 24,326)")
cov = {"Education": "edu", "Income-to-poverty ratio": "pir", "Marital status": "mar",
       "Health insurance": "insured", "Food security": "foodinsec", "Nativity": "usborn",
       "Smoking status": "smoke", "Serum cotinine": "cotinine", "Alcohol": "drink",
       "Body mass index": "bmi", "Sedentary time": "sed", "Recreational MVPA": "mvpa",
       "Sleep duration": "sleep", "Chronic conditions": "ncond", "Diabetes": "dm",
       "eGFR": "egfr", "Energy intake": "kcal", "Dietary folate": "d_folate",
       "Dietary folic acid": "d_fa", "Any supplement use": "supp",
       "Fasting hours": "fasthr", "MEC session": "session", "hsCRP": "crp",
       "Mean corpuscular volume": "mcv", "Antidepressant use": "ad"}
TYPE = {"edu": "categorical", "mar": "categorical", "smoke": "categorical",
        "drink": "categorical", "session": "categorical", "race": "categorical",
        "insured": "binary", "foodinsec": "binary", "usborn": "binary",
        "dm": "binary", "supp": "binary",
        "crp": "structurally unavailable in 2011-2014",
        "ad": "structurally unavailable in Aug 2021-Aug 2023"}
a = d[d.ins == 1]
mi = pd.DataFrame([{"Covariate": k, "Type": TYPE.get(v, "continuous"),
                    "Missing, n": int(a[v].isna().sum()),
                    "Missing, %": round(100 * a[v].isna().mean(), 1)} for k, v in cov.items()])
print(mi.to_string(index=False))
OUT["missingness"] = mi.to_dict("records")

# ------------------------------------------------ 2. DAG-minimal model
print("\n" + "=" * 100)
print("MODEL SPECIFICATIONS")
rows = [umfa_row("Model 1: demographic", level=1),
        umfa_row("Model 2: + socioeconomic", level=2),
        umfa_row("Model 3: extended (primary)", level=3),
        umfa_row("Model 5: DAG-guided reduced adjustment set", level=5),
        umfa_row("Model 4: model 3 + hsCRP and MCV",
                 mask=d.crp.notna().to_numpy(), level=4)]
t = pd.DataFrame(rows)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
OUT["specifications"] = t.to_dict("records")

# ------------------------------------------------ 3. cycle / measurement-mode sensitivity
print("\n" + "=" * 100)
print("SURVEY-CYCLE AND MEASUREMENT-MODE SENSITIVITY")
no_supp_cov = [d.l_kcal.to_numpy()]
rows = [
    umfa_row("All cycles (primary)"),
    umfa_row("Excluding Aug 2021-Aug 2023 (supplement/medication collection changed)",
             mask=(d.cycle != "Aug 2021-Aug 2023").to_numpy()),
    umfa_row("Aug 2021-Aug 2023 only", mask=(d.cycle == "Aug 2021-Aug 2023").to_numpy()),
    umfa_row("Excluding 2011-2012", mask=(d.cycle != "2011-2012").to_numpy()),
]
# reduced model that uses no dietary-supplement covariates at all
X, nm = X_model(d, 1)
for c_, n_ in [(d.l_cot.fillna(d.l_cot[ins].mean()).to_numpy(), "log cotinine"),
               (d.egfr.fillna(d.egfr[ins].mean()).to_numpy(), "eGFR"),
               (d.fasthr.fillna(d.fasthr[ins].mean()).to_numpy(), "Fasting hours")]:
    X, nm = add(X, nm, [c_], [n_])
X2, nm2 = add(X, nm, [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy()], ["UMFA", "5-MTHF"])
res, V, n, df = svylogit(des, d.dep.to_numpy(float), X2, mask=ins, names=nm2)
r = res[res.term == "UMFA"].iloc[0]
rows.append(dict(analysis="Reduced model with no dietary or supplement covariates",
                 n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p), _or=float(r.OR),
                 _lo=float(r.lo), _hi=float(r.hi), _p=float(r.p)))
t = pd.DataFrame(rows)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
OUT["cycle_mode"] = t.to_dict("records")

# ------------------------------------------------ 4. supplement / fasting strata
print("\n" + "=" * 100)
print("FOLIC-ACID-SPECIFIC SUPPLEMENT AND FASTING STRATA")
rows = [
    umfa_row("Folic-acid-containing supplement users", mask=(d.fa_supp == 1).to_numpy()),
    umfa_row("Non-users of folic-acid-containing supplements", mask=(d.fa_supp == 0).to_numpy()),
    umfa_row("Non-users of any dietary supplement", mask=(d.supp == 0).to_numpy()),
    umfa_row("Fasted >= 8 h", mask=(d.fasthr >= 8).to_numpy()),
    umfa_row("Fasted >= 10 h", mask=(d.fasthr >= 10).to_numpy()),
    umfa_row("Fasted >= 12 h", mask=(d.fasthr >= 12).to_numpy()),
    umfa_row("Morning examination session only", mask=(d.session == 1).to_numpy()),
]
t = pd.DataFrame(rows)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
OUT["supp_fast"] = t.to_dict("records")
# interaction with folic-acid supplement use
cc = d.fa_supp.fillna(0).to_numpy(float)
res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), cc,
                            d.l_umfa_z.to_numpy() * cc],
                    ["UMFA", "5-MTHF", "FA supplement", "UMFA x FA supplement"], 3)
r = res[res.term == "UMFA x FA supplement"].iloc[0]
print(f"  interaction UMFA x folic-acid supplement: OR {fmt(r.OR, r.lo, r.hi)}, p={pfmt(r.p)}")
OUT["fa_supp_interaction"] = dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), p=float(r.p))

# ------------------------------------------------ 5. multiple imputation (MICE)
print("\n" + "=" * 100)
print("MULTIPLE IMPUTATION (chained equations, m = 10)")
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import BayesianRidge

num = ["pir", "insured", "foodinsec", "usborn", "l_cot", "bmi", "sed", "mvpa", "sleep",
       "ncond", "dm", "egfr", "l_kcal", "l_dfol", "l_dfa", "supp", "l_sfa", "fasthr"]
catv = {"edu": ["<High school", "High school/GED", "Some college", "College graduate+"],
        "mar": ["Married/partnered", "Widowed/divorced/separated", "Never married"],
        "smoke": ["Never", "Former", "Current"],
        "drink": ["None", "Occasional", "Frequent"]}
base = d[num].copy()
for c_, lv in catv.items():
    for k, level in enumerate(lv):
        s = d[c_]
        base[f"{c_}__{k}"] = np.where(s.isna(), np.nan, (s.astype(str) == level).astype(float))
base["age"] = d.age
base["female"] = d.female
base["l_umfa_z"] = d.l_umfa_z
base["l_mthf_z"] = d.l_mthf_z
base["l_mefox_z"] = d.l_mefox_z
base["dep"] = d.dep
ests, ses = [], []
for m in range(10):
    imp = IterativeImputer(estimator=BayesianRidge(), max_iter=10, sample_posterior=True,
                           random_state=100 + m)
    filled = pd.DataFrame(imp.fit_transform(base.loc[ins]), columns=base.columns,
                          index=base.index[ins])
    dm_ = d.copy()
    for c_ in num:
        dm_.loc[ins, c_] = filled[c_].to_numpy()
    for c_, lv in catv.items():
        arr = filled[[f"{c_}__{k}" for k in range(len(lv))]].to_numpy()
        dm_.loc[ins, c_] = [lv[j] for j in arr.argmax(axis=1)]
    dm_["l_cot"] = dm_["l_cot"]
    dsg = SvyDesign(dm_, "w", "strata", "psu")
    X, nm = X_model(dm_, 3)
    X, nm = add(X, nm, [dm_.l_umfa_z.to_numpy(), dm_.l_mthf_z.to_numpy(),
                        dm_.l_mefox_z.to_numpy()], ["UMFA", "5-MTHF", "MeFox"])
    res, V, n, df = svylogit(dsg, dm_.dep.to_numpy(float), X, mask=ins, names=nm)
    r = res[res.term == "UMFA"].iloc[0]
    ests.append(r.beta)
    ses.append(r.se)
ests, ses = np.array(ests), np.array(ses)
qbar = ests.mean()
U = (ses ** 2).mean()
B = ests.var(ddof=1)
T = U + (1 + 1 / len(ests)) * B
sepool = np.sqrt(T)
nu = (len(ests) - 1) * (1 + U / ((1 + 1 / len(ests)) * B)) ** 2
tq = st.t.ppf(0.975, nu)
print(f"  pooled OR {np.exp(qbar):.3f} ({np.exp(qbar - tq*sepool):.3f}-"
      f"{np.exp(qbar + tq*sepool):.3f}), p={2*st.t.sf(abs(qbar/sepool), nu):.4f}  "
      f"(m=10, FMI-adjusted df={nu:.0f})")
OUT["mice"] = dict(OR=float(np.exp(qbar)), lo=float(np.exp(qbar - tq * sepool)),
                   hi=float(np.exp(qbar + tq * sepool)),
                   p=float(2 * st.t.sf(abs(qbar / sepool), nu)), m=len(ests))

# ------------------------------------------------ 6. FDR across secondary analyses
print("\n" + "=" * 100)
print("BENJAMINI-HOCHBERG FDR: secondary outcomes and alternative exposure parameterisations")
A1 = json.load(open("out_a1.json"))
A2 = json.load(open("out_a2.json"))
A3 = json.load(open("out_a3.json"))
sec = []
for r in A1["outcomes"]:
    if r["outcome"] != "PHQ-9 >= 10":
        sec.append((f"Outcome: {r['outcome']}", float(r["_UMFA"][3])))
A4j = json.load(open("out_a4.json"))
# Every quartile inference that Table 3 reports enters the family, not only the extreme
# contrast: the joint test and both trend tests are reported estimates too, and excluding
# the nominally significant ones would make the correction self-serving.
for r in A1["quartiles"]:
    ex = r["exposure"]
    sec.append((f"Quartile Q4 vs Q1: {ex}", float(r["_q4"][3])))
    sec.append((f"Quartile joint test: {ex}", float(r["_p_joint"])))
    sec.append((f"Quartile trend (quartile median): {ex}", float(r["_p_trend"])))
    sec.append((f"Quartile trend (quartile rank): {ex}", float(A4j["ordinal_trend"][ex]["p"])))
for r in A3["thresholds"]:
    sec.append((f"Threshold: {r['threshold']}", float(r["_p"])))
def _p(rec):
    if "_p" in rec:
        return float(rec["_p"])
    s = str(rec.get("p", "1"))
    return 0.0005 if s.startswith("<") else float(s)


# Subgroup strata and robustness checks are deliberately NOT included: subgroup inference
# belongs to the interaction tests, and sensitivity analyses re-examine one hypothesis
# rather than testing new ones.
lab = [s[0] for s in sec]
p = np.array([s[1] for s in sec])
o = np.argsort(p)
q = np.empty_like(p)
mtot = len(p)
prev = 1.0
for rank, idx in enumerate(o[::-1]):
    k = mtot - rank
    prev = min(prev, p[idx] * mtot / k)
    q[idx] = prev
fdr = pd.DataFrame({"Analysis": lab, "p": [f"{x:.3f}" for x in p],
                    "BH q": [f"{x:.3f}" for x in q]})
print(fdr.to_string(index=False))
OUT["fdr"] = fdr.to_dict("records")

# ------------------------------------------------ 7. survey estimation validation
print("\n" + "=" * 100)
print("VALIDATION OF THE DESIGN-BASED ESTIMATOR AGAINST PUBLISHED NCHS FIGURES")
val = []
mkl = ins & (d.cyc == "L").to_numpy()
p1 = des.subset(mkl).prop(d.dep.to_numpy())
val.append(dict(estimate="PHQ-9 >= 10, analytic sample, Aug 2021-Aug 2023",
                ours=f"{100*p1['est']:.1f}% ({100*p1['lo']:.1f}-{100*p1['hi']:.1f})",
                published="same as Supplementary Table S2 by construction"))
# wider denominator: all examined adults with a complete PHQ-9, MEC weight, this cycle
dl = d.copy()
dl["w"] = np.where((d.cyc == "L") & (d.age >= 18) & d.phq9.notna(),
                   d.WTMEC2YR.fillna(0.0) if "WTMEC2YR" in d else d.wraw.fillna(0.0), 0.0)
desl = SvyDesign(dl, "w", "strata", "psu")
mk2 = ((d.cyc == "L") & (d.age >= 18) & d.phq9.notna() & (dl.w > 0)).to_numpy()
p2 = desl.subset(mk2).prop(d.dep.to_numpy())
val.append(dict(estimate="PHQ-9 >= 10, all examined adults 18+ with complete PHQ-9, MEC weight",
                ours=f"{100*p2['est']:.1f}% ({100*p2['lo']:.1f}-{100*p2['hi']:.1f})",
                published="13.1% for ages 12+ (NCHS Data Brief 527)"))
# Independent replication of the primary model. The point estimates are re-derived by a
# quasi-Newton maximisation of the same weighted likelihood, written independently of the
# IRLS routine; the Taylor-linearised standard error is checked against a delete-one-PSU
# jackknife. Statsmodels' GLM is not used because raw NHANES weights as frequency weights
# make its optimiser diverge.
from scipy.optimize import minimize
X, nmv = X_model(d, 3)
X, nmv = add(X, nmv, [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()],
             ["UMFA", "5-MTHF", "MeFox"])
res0, V0, n0, df0 = svylogit(des, d.dep.to_numpy(float), X, mask=ins, names=nmv)
r0 = res0[res0.term == "UMFA"].iloc[0]
Xk = np.nan_to_num(X)[:, [nmv.index(t) for t in res0.term]][ins]
yv = d.dep.to_numpy(float)[ins]
wv = des.w[ins] / des.w[ins].mean()
RIDGE = 1e-8


def _nll(b):
    z = np.clip(Xk @ b, -500, 500)
    return -np.sum(wv * (yv * z - np.logaddexp(0, z))) + RIDGE * np.sum(b[1:] ** 2)


def _grad(b):
    z = np.clip(Xk @ b, -500, 500)
    mu = 1 / (1 + np.exp(-z))
    g = -(Xk.T @ (wv * (yv - mu)))
    g[1:] += 2 * RIDGE * b[1:]
    return g


opt = minimize(_nll, np.zeros(Xk.shape[1]), jac=_grad, method="L-BFGS-B",
               options=dict(maxiter=20000, maxfun=40000, ftol=1e-16, gtol=1e-12))
jU = list(res0.term).index("UMFA")
key_idx = [list(res0.term).index(t) for t in ["UMFA", "5-MTHF", "MeFox", "Age, per 10 y", "Female"]]
maxdiff = float(np.max(np.abs(res0.beta.to_numpy()[key_idx] - opt.x[key_idx])))
rank = int(np.linalg.matrix_rank(Xk))
val.append(dict(estimate="Design matrix rank",
                ours=f"{rank} of {Xk.shape[1]} columns",
                published="full rank required for a unique maximum likelihood solution"))
val.append(dict(estimate="Primary UMFA log-odds: IRLS vs independent quasi-Newton maximisation "
                         "of the same weighted likelihood",
                ours=f"{r0.beta:.6f} vs {opt.x[jU]:.6f} (largest absolute difference across the "
                     f"substantive coefficients {maxdiff:.0e})",
                published="agreement to numerical precision expected"))
# delete-one-PSU jackknife
psus = pd.Series(list(zip(d.strata, d.psu)))[ins].unique()
nh = pd.Series([p_[0] for p_ in psus]).value_counts().to_dict()
reps, reps_h = [], []
for p_ in psus:
    keep = ins & ~((d.strata == p_[0]) & (d.psu == p_[1])).to_numpy()
    f = nh[p_[0]] / (nh[p_[0]] - 1.0) if nh[p_[0]] > 1 else 1.0
    dw = d.copy()
    dw["w"] = np.where(keep & (d.strata == p_[0]).to_numpy(), d.w * f,
                       np.where(keep, d.w, 0.0))
    try:
        rr, _, _, _ = svylogit(SvyDesign(dw, "w", "strata", "psu"), d.dep.to_numpy(float),
                               X, mask=keep, names=nmv)
        reps.append(float(rr[rr.term == "UMFA"].iloc[0].beta))
        reps_h.append(p_[0])
    except Exception:
        pass
reps = np.array(reps)
reps_h = np.array(reps_h)
jk_var = 0.0
for h in np.unique(reps_h):
    b_h = reps[reps_h == h]
    if len(b_h) > 1:
        jk_var += (len(b_h) - 1) / len(b_h) * np.sum((b_h - b_h.mean()) ** 2)
jk_se = float(np.sqrt(jk_var))
val.append(dict(estimate="Primary UMFA standard error: Taylor linearisation vs delete-one-PSU "
                         "jackknife",
                ours=f"{r0.se:.5f} vs {jk_se:.5f}",
                published="close agreement supports the variance estimator"))
# antidepressant use as an outcome, so the figure quoted in the text has a home
_rx = ins & d.ad.notna().to_numpy()
_res, _V, _n, _df = fit("ad", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                               d.l_mefox_z.to_numpy()], ["UMFA", "5-MTHF", "MeFox"], 3, mask=_rx)
_r = _res[_res.term == "UMFA"].iloc[0]
OUT["ad_as_outcome"] = dict(n=int(_n), OR=fmt(_r.OR, _r.lo, _r.hi), p=pfmt(_r.p))
print("  UMFA -> antidepressant use:", OUT["ad_as_outcome"])
for v in val:
    print(" ", v)
OUT["validation"] = val

with open("out_a5.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a5.json")
