"""Second review round: survey-weighted quantiles and SDs, spline-adjusted confounders,
selection check for incomplete PHQ-9, discovery/replication split, Model-3 heterogeneity."""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from common import (get, design, X_model, add, fmt, pfmt, wquantile, wmean_sd, CYCLES)
from survey import svylogit, wald, contrast

d = get()
des = design(d)
ins = (d.ins == 1).to_numpy()
w = d.w.to_numpy()
OUT = {}
pd.set_option("display.width", 250)


def fit(y, cols, names, level, mask=None):
    X, nm = X_model(d, level)
    X, nm = add(X, nm, cols, names)
    mk = ins if mask is None else (ins & mask)
    return svylogit(des, d[y].to_numpy(float), X, mask=mk, names=nm)


def prim(cols, names, level=3, mask=None, term="UMFA"):
    res, V, n, df = fit("dep", cols, names, level, mask=mask)
    r = res[res.term == term].iloc[0]
    return dict(n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p), _or=float(r.OR),
                _lo=float(r.lo), _hi=float(r.hi), _p=float(r.p), _b=float(r.beta),
                _se=float(r.se))


TRIO = ([d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()],
        ["UMFA", "5-MTHF", "MeFox"])

# ---------------------------------------------- 1. survey-weighted standardisation
print("=" * 100)
print("SURVEY-WEIGHTED vs UNWEIGHTED WITHIN-CYCLE STANDARDISATION")
for v in ["l_umfa", "l_mthf", "l_mefox"]:
    zc = pd.Series(np.nan, index=d.index)
    for cy in CYCLES:
        sel = (d.cycle == cy).to_numpy()
        mu, sd = wmean_sd(d.loc[sel & ins, v], w[sel & ins])
        zc[sel] = (d.loc[sel, v] - mu) / sd
    d[v + "_wz"] = zc
rows = [dict(analysis="Unweighted within-cycle SD (primary)", **prim(*TRIO)),
        dict(analysis="Survey-weighted within-cycle mean and SD",
             **prim([d.l_umfa_wz.to_numpy(), d.l_mthf_wz.to_numpy(), d.l_mefox_wz.to_numpy()],
                    ["UMFA", "5-MTHF", "MeFox"]))]
sdtab = []
for cy in CYCLES:
    sel = (d.cycle == cy).to_numpy() & ins
    mu_u, sd_u = d.loc[sel, "l_umfa"].mean(), d.loc[sel, "l_umfa"].std()
    mu_w, sd_w = wmean_sd(d.loc[sel, "l_umfa"], w[sel])
    sdtab.append(dict(cycle=cy, sd_unweighted=round(sd_u, 4), sd_weighted=round(sd_w, 4)))
t = pd.DataFrame(rows)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
print(pd.DataFrame(sdtab).to_string(index=False))
OUT["weighted_sd"] = dict(models=t.to_dict("records"), sds=sdtab)

# ---------------------------------------------- 2. survey-weighted quantiles
print("\n" + "=" * 100)
print("SURVEY-WEIGHTED CYCLE-SPECIFIC QUARTILES AND UPPER-TAIL CUT-POINTS")
q_w = pd.Series(np.nan, index=d.index)
p90w = pd.Series(np.nan, index=d.index)
p95w = pd.Series(np.nan, index=d.index)
for cy in CYCLES:
    sel = (d.cycle == cy).to_numpy()
    cut = wquantile(d.loc[sel & ins, "l_umfa"], w[sel & ins], [.25, .5, .75, .90, .95])
    x = d.loc[sel, "l_umfa"]
    q_w[sel] = np.select([x <= cut[0], x <= cut[1], x <= cut[2], x > cut[2]], [1, 2, 3, 4],
                         default=np.nan)
    p90w[sel] = (x > cut[3]).astype(float)
    p95w[sel] = (x > cut[4]).astype(float)
d["umfa_qw"] = q_w
cols = [(d.umfa_qw == k).astype(float).to_numpy() for k in [2, 3, 4]]
nms = [f"Q{k}" for k in [2, 3, 4]]
res, V, n, df = fit("dep", cols + [d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()],
                    nms + ["5-MTHF", "MeFox"], 3, mask=d.umfa_qw.notna().to_numpy())
idx = [list(res.term).index(x) for x in nms]
F, k_, d2, pj = wald(res.beta.to_numpy(), V, idx, df)
row = {"quartiles": "survey-weighted", "Q1": "1.00 (ref)"}
for kk, x in zip([2, 3, 4], nms):
    r = res[res.term == x].iloc[0]
    row[f"Q{kk}"] = fmt(r.OR, r.lo, r.hi)
row["p_joint"] = pfmt(pj)
print(pd.DataFrame([row]).to_string(index=False))
thr = []
for v, lab in [("p90w", "UMFA > weighted 90th pct"), ("p95w", "UMFA > weighted 95th pct")]:
    d[v] = p90w if v == "p90w" else p95w
    r = prim([d[v].to_numpy(float), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()],
             [lab, "5-MTHF", "MeFox"], term=lab)
    r["threshold"] = lab
    r["n_exposed"] = int((ins & (d[v] == 1).to_numpy()).sum())
    thr.append(r)
tt = pd.DataFrame(thr)
print(tt[["threshold", "n_exposed", "OR", "p"]].to_string(index=False))
OUT["weighted_quantiles"] = dict(quartiles=row, thresholds=tt.to_dict("records"))

# ---------------------------------------------- 3. spline-adjusted confounders
print("\n" + "=" * 100)
print("NON-LINEAR ADJUSTMENT FOR KEY CONFOUNDERS")


def rcs(x, knots):
    k = len(knots)
    t = knots
    out = []
    for j in range(k - 2):
        def cu(v, kn):
            return np.maximum(v - kn, 0) ** 3
        num = (cu(x, t[j]) - cu(x, t[k - 2]) * (t[k - 1] - t[j]) / (t[k - 1] - t[k - 2])
               + cu(x, t[k - 1]) * (t[k - 2] - t[j]) / (t[k - 1] - t[k - 2]))
        out.append(num / (t[k - 1] - t[0]) ** 2)
    return out


spl_cols, spl_nms = [], []
for var, lab in [("age", "age"), ("egfr", "eGFR"), ("l_dfa", "dietary folic acid"),
                 ("l_sfa", "supplemental folic acid"), ("fasthr", "fasting hours")]:
    x = d[var].astype(float)
    x = x.fillna(x[ins].mean()).to_numpy()
    kn = np.nanpercentile(x[ins], [10, 50, 90])
    for j, b in enumerate(rcs(x, kn)):
        spl_cols.append(np.nan_to_num(b))
        spl_nms.append(f"{lab} spline {j+1}")
r = prim(TRIO[0] + spl_cols, TRIO[1] + spl_nms)
print(f"  Model 3 + restricted cubic splines for age, eGFR, dietary and supplemental "
      f"folic acid and fasting time: OR {r['OR']}, p={r['p']} (n={r['n']})")
OUT["spline_confounders"] = r

# ---------------------------------------------- 4. selection: incomplete PHQ-9
print("\n" + "=" * 100)
print("PARTICIPANTS EXCLUDED FOR INCOMPLETE PHQ-9 vs INCLUDED")
elig = ((d.age >= 18) & (d.preg != 1) & (d.w > 0) & d.mthf.notna() & d.umfa.notna()).to_numpy()
excl = elig & d.phq9.isna().to_numpy()
incl = elig & d.phq9.notna().to_numpy()
rows = []
for lab, var, kind in [("Age, years", "age", "m"), ("Female, %", "female", "p"),
                       ("NH White, %", "race", "race"), ("UMFA, nmol/L (GM)", "l_umfa", "g"),
                       ("5-MTHF, nmol/L (GM)", "l_mthf", "g"),
                       ("Any supplement, %", "supp", "p"),
                       ("Income-to-poverty ratio", "pir", "m")]:
    def val(mk):
        if kind == "m":
            return f"{np.nanmean(d.loc[mk, var]):.2f}"
        if kind == "p":
            return f"{100*np.nanmean(d.loc[mk, var]):.1f}"
        if kind == "g":
            return f"{np.exp(np.nanmean(d.loc[mk, var])):.2f}"
        return f"{100*(d.loc[mk, var] == 'NH White').mean():.1f}"
    rows.append({"Characteristic": lab, "Included (n=24,326)": val(incl),
                 "Excluded for incomplete PHQ-9 (n=2,677)": val(excl)})
by_cy = pd.crosstab(d.loc[elig, "cycle"], d.loc[elig, "phq9"].isna(), normalize="index")
sel = pd.DataFrame(rows)
print(sel.to_string(index=False))
print("\nProportion excluded by cycle:")
print((100 * by_cy[True]).round(1).to_string())
OUT["selection"] = dict(table=sel.to_dict("records"),
                        by_cycle={k: round(100 * v, 1) for k, v in by_cy[True].items()})
# UMFA distribution difference
tt_ = st.ttest_ind(d.loc[incl, "l_umfa"].dropna(), d.loc[excl, "l_umfa"].dropna(),
                   equal_var=False)
print(f"  log-UMFA included vs excluded: t = {tt_.statistic:.2f}, p = {tt_.pvalue:.3f}")
OUT["selection"]["umfa_t_p"] = float(tt_.pvalue)

# ---------------------------------------------- 5. discovery / replication split
print("\n" + "=" * 100)
print("DISCOVERY (screened cycle) vs OUT-OF-SCREEN REPLICATION vs POOLED")
rows = [dict(sample="Discovery: Aug 2021-Aug 2023 (cycle used for the exposure screen)",
             **prim(*TRIO, mask=(d.cycle == "Aug 2021-Aug 2023").to_numpy())),
        dict(sample="Out-of-screen replication: 2011-March 2020 (four independent cycles)",
             **prim(*TRIO, mask=(d.cycle != "Aug 2021-Aug 2023").to_numpy())),
        dict(sample="Pooled: 2011-2023", **prim(*TRIO))]
t = pd.DataFrame(rows)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
OUT["discovery_replication"] = t.to_dict("records")

# ---------------------------------------------- 6. Model-3 heterogeneity across cycles
print("\n" + "=" * 100)
print("BETWEEN-CYCLE HETEROGENEITY, PRIMARY MODEL 3 SPECIFICATION")
bs, ses, labs = [], [], []
for cy in CYCLES:
    r = prim(*TRIO, mask=(d.cycle == cy).to_numpy())
    bs.append(r["_b"]); ses.append(r["_se"]); labs.append(cy)
    print(f"  {cy:20s} OR {r['OR']}  p={r['p']}")
bs, ses = np.array(bs), np.array(ses)
wt = 1 / ses ** 2
bbar = (wt * bs).sum() / wt.sum()
Q = (wt * (bs - bbar) ** 2).sum()
I2 = max(0.0, 100 * (Q - 4) / Q)
print(f"  pooled OR {np.exp(bbar):.2f}; Cochran Q = {Q:.2f} (df=4), p = {st.chi2.sf(Q,4):.3f}, "
      f"I2 = {I2:.0f}%")
OUT["heterogeneity_m3"] = dict(pooled=float(np.exp(bbar)), Q=float(Q),
                               p=float(st.chi2.sf(Q, 4)), I2=float(I2),
                               per_cycle=[dict(cycle=c, b=float(b), se=float(s))
                                          for c, b, s in zip(labs, bs, ses)])

# ---------------------------------------------- 7. per-biomarker n for Table 2
print("\n" + "=" * 100)
print("ANALYTIC n BY BIOMARKER")
ns = {}
for v, lab in [("l_mthf", "5-Methyltetrahydrofolate"), ("l_tfol", "Total serum folate"),
               ("l_rbcfol", "Red-cell folate"), ("l_mefox", "MeFox"),
               ("l_umfa", "Unmetabolised folic acid"),
               ("l_pct_umfa", "UMFA fraction (% of quantified folate)")]:
    ns[lab] = int((ins & d[v].notna().to_numpy()).sum())
    print(f"  {lab:45s} n = {ns[lab]:,}")
OUT["biomarker_n"] = ns


# ---------------------------------------------- 8. by-cycle descriptive table (S2)
rows = []
for cy in CYCLES:
    mk = ins & (d.cycle == cy).to_numpy()
    r = {"Survey cycle": cy, "n": f"{int(mk.sum()):,}"}
    for v, lab in [("l_mthf", "5-MTHF"), ("l_umfa", "UMFA"), ("l_mefox", "MeFox"),
                   ("l_tfol", "Total folate"), ("l_rbcfol", "RBC folate"),
                   ("l_pct_umfa", "UMFA fraction, %")]:
        e = des.subset(mk).mean(d[v].to_numpy())
        r[lab] = f"{np.exp(e[0]):.2f} ({np.exp(e[2]):.2f}-{np.exp(e[3]):.2f})"
    for y, lab in [("dep", "PHQ-9>=10, %"), ("dep15", "PHQ-9>=15, %"),
                   ("si", "Suicidal ideation, %")]:
        pp = des.subset(mk).prop(d[y].to_numpy())
        r[lab] = f"{100*pp['est']:.1f} ({100*pp['lo']:.1f}-{100*pp['hi']:.1f})"
    rows.append(r)
OUT["by_cycle_desc"] = rows

with open("out_a6.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a6.json")

# ---------------------------------------------- 9. biomarker correlations and collinearity
print("\n" + "=" * 100)
print("BIOMARKER CORRELATION MATRIX AND COLLINEARITY OF THE PRIMARY MODEL")
cormat = d.loc[ins, ["l_mthf", "l_umfa", "l_mefox", "l_tfol", "l_rbcfol", "l_pct_umfa"]].corr()
cormat.index = cormat.columns = ["5-MTHF", "UMFA", "MeFox", "Total folate", "RBC folate",
                                 "UMFA fraction"]
print(cormat.round(2).to_string())
OUT["corr"] = cormat.round(2).reset_index().rename(columns={"index": ""}).to_dict("records")
Xc, nmc = X_model(d, 3)
Xc, nmc = add(Xc, nmc, TRIO[0], TRIO[1])
Xk = np.nan_to_num(Xc)[ins]
keep = [j for j in range(Xk.shape[1]) if np.std(Xk[:, j]) > 1e-10]
Xs = Xk[:, keep]
Xs = (Xs - Xs.mean(0)) / Xs.std(0)
# drop columns that are exact linear combinations of earlier ones before conditioning
_q, _r, _piv = __import__("scipy.linalg", fromlist=["qr"]).qr(Xs, mode="economic", pivoting=True)
_rank = int(np.sum(np.abs(np.diag(_r)) > 1e-8 * abs(_r[0, 0])))
cn = float(np.linalg.cond(Xs[:, sorted(_piv[:_rank])]))
vif = {}
for t in TRIO[1]:
    j = [nmc[k] for k in keep].index(t)
    others = [k for k in range(Xs.shape[1]) if k != j]
    b_, *_ = np.linalg.lstsq(Xs[:, others], Xs[:, j], rcond=None)
    r2 = 1 - np.sum((Xs[:, j] - Xs[:, others] @ b_) ** 2) / np.sum(Xs[:, j] ** 2)
    vif[t] = round(float(1 / max(1 - r2, 1e-12)), 2)
print(f"  design-matrix rank {_rank} of {Xs.shape[1]}; condition number (full-rank subset) {cn:.1f}")
print("  variance inflation factors:", vif)
OUT["collinearity"] = dict(condition_number=round(cn, 1), rank=int(_rank),
                           ncol=int(Xs.shape[1]), vif=vif)

# ---------------------------------------------- 10. ordinal PHQ-9 severity
print("\n" + "=" * 100)
print("ORDINAL PHQ-9 SEVERITY (proportional-odds, fitted as cumulative logits)")
cuts = [(5, "PHQ-9 >= 5"), (10, "PHQ-9 >= 10"), (15, "PHQ-9 >= 15")]
orow = []
for c_, lab in cuts:
    y = (d.phq9 >= c_).astype(float).where(d.phq9.notna())
    res, V, n, df = fit_y = None, None, None, None
    X, nm = X_model(d, 3)
    X, nm = add(X, nm, TRIO[0], TRIO[1])
    res, V, n, df = svylogit(des, y.to_numpy(float), X, mask=ins, names=nm)
    r = res[res.term == "UMFA"].iloc[0]
    orow.append(dict(cutpoint=lab, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                     _b=float(r.beta), _se=float(r.se)))
t = pd.DataFrame(orow)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
bs = np.array([r["_b"] for r in orow]); ses = np.array([r["_se"] for r in orow])
wt = 1 / ses ** 2
bbar = (wt * bs).sum() / wt.sum()
Q = (wt * (bs - bbar) ** 2).sum()
print(f"  common (proportional) odds ratio across cut-points: {np.exp(bbar):.3f}; "
      f"test of proportional odds Q = {Q:.2f} (df=2), p = {st.chi2.sf(Q, 2):.3f}")
OUT["ordinal"] = dict(rows=t.to_dict("records"), common_or=float(np.exp(bbar)),
                      Q=float(Q), p=float(st.chi2.sf(Q, 2)))

# ---------------------------------------------- 11. selection: no folate assay
print("\n" + "=" * 100)
print("PARTICIPANTS EXCLUDED BECAUSE FOLATE VITAMERS WERE NOT QUANTIFIED")
elig2 = ((d.age >= 18) & (d.preg != 1) & (d.w > 0)).to_numpy()
no_fol = elig2 & (d.mthf.isna() | d.umfa.isna()).to_numpy()
has_fol = elig2 & d.mthf.notna().to_numpy() & d.umfa.notna().to_numpy()
print(f"  n without folate assay = {no_fol.sum()}, with = {has_fol.sum()}")
rows = []
for lab, var, kind in [("Age, years", "age", "m"), ("Female, %", "female", "p"),
                       ("NH White, %", "race", "race"),
                       ("Income-to-poverty ratio", "pir", "m"),
                       ("Any supplement, %", "supp", "p"),
                       ("PHQ-9 >= 10, %", "dep", "p"),
                       ("PHQ-9 score", "phq9", "m")]:
    def v2(mk):
        if kind == "m":
            return f"{np.nanmean(d.loc[mk, var]):.2f}"
        if kind == "p":
            return f"{100*np.nanmean(d.loc[mk, var]):.1f}"
        return f"{100*(d.loc[mk, var] == 'NH White').mean():.1f}"
    rows.append({"Characteristic": lab, "Folate assayed": v2(has_fol),
                 "Not assayed": v2(no_fol)})
tsel2 = pd.DataFrame(rows)
print(tsel2.to_string(index=False))
by_cy2 = pd.crosstab(d.loc[elig2, "cycle"], (d.loc[elig2, "mthf"].isna() |
                                             d.loc[elig2, "umfa"].isna()), normalize="index")
OUT["selection_folate"] = dict(table=tsel2.to_dict("records"),
                               n_missing=int(no_fol.sum()), n_present=int(has_fol.sum()),
                               by_cycle={k: round(100 * v, 1) for k, v in by_cy2[True].items()})
print("  proportion without folate assay by cycle:")
print((100 * by_cy2[True]).round(1).to_string())

# ---------------------------------------------- 12. weighted comparison for S16
print("\n" + "=" * 100)
print("WEIGHTED COMPARISON, INCLUDED vs EXCLUDED FOR INCOMPLETE PHQ-9")
rows = []
for lab, var in [("Age, years", "age"), ("Female", "female"),
                 ("Income-to-poverty ratio", "pir"), ("Any supplement", "supp"),
                 ("UMFA, log nmol/L", "l_umfa"), ("5-MTHF, log nmol/L", "l_mthf")]:
    m1, s1 = wmean_sd(d.loc[incl, var].astype(float), d.loc[incl, "wraw"].fillna(0))
    m2, s2 = wmean_sd(d.loc[excl, var].astype(float), d.loc[excl, "wraw"].fillna(0))
    sp = np.sqrt((s1 ** 2 + s2 ** 2) / 2)
    rows.append(dict(Characteristic=lab, Included=round(float(m1), 2),
                     Excluded=round(float(m2), 2),
                     SMD=round(float((m2 - m1) / sp) if sp > 0 else np.nan, 3)))
tw = pd.DataFrame(rows)
print(tw.to_string(index=False))
OUT["selection_weighted"] = tw.to_dict("records")

# ---------------------------------------------- 13. absolute risk difference with CI
print("\n" + "=" * 100)
print("ADJUSTED ABSOLUTE RISK DIFFERENCE, 99th vs 10th PERCENTILE OF UMFA")
X, nm = X_model(d, 3)
X, nm = add(X, nm, [d.l_umfa.to_numpy(float), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()],
            ["logUMFA", "5-MTHF", "MeFox"])
res, V, n, df = svylogit(des, d.dep.to_numpy(float), X, mask=ins, names=nm)
terms = list(res.term); bb = res.beta.to_numpy()
j = terms.index("logUMFA")
Xk = np.nan_to_num(X)[:, [nm.index(t) for t in terms]]
wv = des.w * ins
q10, q99 = np.nanpercentile(d.l_umfa[ins], [10, 99])
def marg(q):
    Xs = Xk.copy(); Xs[:, j] = q
    mu = 1 / (1 + np.exp(-(Xs @ bb)))
    est = np.sum(wv * mu) / np.sum(wv)
    g = ((wv * mu * (1 - mu))[:, None] * Xs).sum(axis=0) / np.sum(wv)
    return est, g
p10, g10 = marg(q10); p99, g99 = marg(q99)
gd = g99 - g10
se_d = float(np.sqrt(gd @ V @ gd))
print(f"  p10 {100*p10:.2f}%  p99 {100*p99:.2f}%  difference {100*(p99-p10):.2f} pp "
      f"(95% CI {100*(p99-p10-1.96*se_d):.2f} to {100*(p99-p10+1.96*se_d):.2f})")
OUT["risk_diff"] = dict(p10=float(100 * p10), p99=float(100 * p99),
                        diff=float(100 * (p99 - p10)),
                        lo=float(100 * (p99 - p10 - 1.96 * se_d)),
                        hi=float(100 * (p99 - p10 + 1.96 * se_d)))

with open("out_a6.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a6.json (extended)")
