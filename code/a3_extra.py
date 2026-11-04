"""Reverse-causation checks, intake-independent (residual) UMFA, B12 interaction,
adjusted predicted prevalence, and antidepressant use as a secondary outcome."""
import json
import numpy as np
import pandas as pd
from common import get, design, X_model, add, cont, fmt, pfmt, CYCLES
from survey import svylogit, svylin, wald, contrast

d = get()
des = design(d)
ins = (d.ins == 1).to_numpy()
OUT = {}
pd.set_option("display.width", 250)


def fit(y, cols, names, level, mask=None):
    X, nm = X_model(d, level)
    X, nm = add(X, nm, cols, names)
    mk = ins if mask is None else (ins & mask)
    return svylogit(des, d[y].to_numpy(float), X, mask=mk, names=nm)


# ---------------------------------------------- 1. medication: adjust, not exclude
print("=" * 100)
print("ANTIDEPRESSANT USE (cycles with drug data: 2011-Mar 2020)")
rx = ins & d.ad.notna().to_numpy()
rows = []
for lab, extra, extranm, mk in [
        ("Drug-data cycles, model 3", [], [], rx),
        ("+ adjusted for antidepressant use", [d.ad.fillna(0).to_numpy()], ["Antidepressant"], rx),
        ("+ adjusted for any psychotropic", [d.ad.fillna(0).to_numpy(), d.anx.fillna(0).to_numpy(),
                                             d.ap.fillna(0).to_numpy()],
         ["Antidepressant", "Anxiolytic", "Antipsychotic"], rx),
        ("Antidepressant users only", [], [], rx & (d.ad == 1).to_numpy()),
        ("Non-users of antidepressants", [], [], rx & (d.ad == 0).to_numpy()),
        ("No antidepressant AND no supplement", [], [],
         rx & (d.ad == 0).to_numpy() & (d.supp == 0).to_numpy())]:
    res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                d.l_mefox_z.to_numpy()] + extra,
                        ["UMFA", "5-MTHF", "MeFox"] + extranm, 3, mask=mk)
    r = res[res.term == "UMFA"].iloc[0]
    rows.append(dict(analysis=lab, n=n, ncase=int(d.dep[mk].sum()), UMFA=fmt(r.OR, r.lo, r.hi),
                     p=pfmt(r.p), _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
t = pd.DataFrame(rows)
print(t.drop(columns=[c for c in t if c.startswith("_")]).to_string(index=False))
OUT["medication"] = t.to_dict("records")

# interaction with antidepressant use
cc = d.ad.fillna(0).to_numpy()
res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                            d.l_mefox_z.to_numpy(), cc,
                            d.l_umfa_z.to_numpy() * cc],
                    ["UMFA", "5-MTHF", "MeFox", "AD", "UMFA x AD"], 3, mask=rx)
r = res[res.term == "UMFA x AD"].iloc[0]
print(f"  interaction UMFA x antidepressant: OR {fmt(r.OR, r.lo, r.hi)}, p={pfmt(r.p)}")
OUT["ad_interaction"] = dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), p=float(r.p))

# antidepressant use as an outcome (does UMFA track treatment rather than symptoms?)
res, V, n, df = fit("ad", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                           d.l_mefox_z.to_numpy()], ["UMFA", "5-MTHF", "MeFox"],
                    3, mask=rx)
r = res[res.term == "UMFA"].iloc[0]
rm = res[res.term == "5-MTHF"].iloc[0]
print(f"  outcome = antidepressant use: UMFA {fmt(r.OR, r.lo, r.hi)} p={pfmt(r.p)}; "
      f"5-MTHF {fmt(rm.OR, rm.lo, rm.hi)} p={pfmt(rm.p)}  (n={n})")
OUT["ad_outcome"] = dict(UMFA=[float(r.OR), float(r.lo), float(r.hi), float(r.p)],
                         MTHF=[float(rm.OR), float(rm.lo), float(rm.hi), float(rm.p)], n=int(n))

# ---------------------------------------------- 2. intake-independent (residual) UMFA
print("\n" + "=" * 100)
print("INTAKE-INDEPENDENT UMFA (residual of log-UMFA on folic-acid intake and handling)")
pred = pd.DataFrame({
    "l_dfa": d.l_dfa, "l_sfa": d.l_sfa, "supp": d.supp, "l_dfol": d.l_dfol,
    "l_kcal": d.l_kcal, "fasthr": d.fasthr, "sess2": (d.session == 2).astype(float),
    "sess3": (d.session == 3).astype(float), "age": d.age, "female": d.female,
    "bmi": d.bmi, "egfr": d.egfr})
# survey-cycle fixed effects: UMFA fell markedly over the period, so cycle-related
# variation would otherwise remain in the residual
for _cy in CYCLES[1:]:
    pred["cyc_" + _cy[:7]] = (d.cycle == _cy).astype(float)
for c in pred.columns:
    pred[c] = pred[c].fillna(pred[c][ins].mean())
Xp = np.column_stack([np.ones(len(d))] + [pred[c].to_numpy() for c in pred.columns])
yv = d.l_umfa.to_numpy(float)
pnames = ["(int)"] + list(pred.columns)
res_l, Vl, nl, dfl = svylin(des, yv, Xp, mask=ins, names=pnames)
kept = [pnames.index(t) for t in res_l.term]
resid = yv - Xp[:, kept] @ res_l.beta.to_numpy()
r2 = 1 - np.nanvar(resid[ins]) / np.nanvar(yv[ins])
d["umfa_resid_z"] = (resid - np.nanmean(resid[ins])) / np.nanstd(resid[ins])
print(f"  intake/handling model R2 = {r2:.3f}")
rows = []
for lvl in [1, 2, 3]:
    res, V, n, df = fit("dep", [d.umfa_resid_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                d.l_mefox_z.to_numpy()],
                        ["Residual UMFA", "5-MTHF", "MeFox"], lvl)
    r = res[res.term == "Residual UMFA"].iloc[0]
    rows.append(dict(model=lvl, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                     _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
tr = pd.DataFrame(rows)
print(tr.drop(columns=[c for c in tr if c.startswith("_")]).to_string(index=False))
OUT["residual_umfa"] = dict(r2=float(r2), models=tr.to_dict("records"))

# ---------------------------------------------- 2b. high-UMFA thresholds
print("\n" + "=" * 100)
print("HIGH-UMFA THRESHOLDS (model 3, adjusted for 5-MTHF)")
a_ = d[d.ins == 1]
for cy in CYCLES:
    sel = d.cycle == cy
    p90 = a_.loc[a_.cycle == cy, "l_umfa"].quantile(0.90)
    p95 = a_.loc[a_.cycle == cy, "l_umfa"].quantile(0.95)
    d.loc[sel, "umfa_p90"] = (d.loc[sel, "l_umfa"] > p90).astype(float)
    d.loc[sel, "umfa_p95"] = (d.loc[sel, "l_umfa"] > p95).astype(float)
d["umfa_ge5"] = (d.umfa >= 5).astype(float).where(d.umfa.notna())
d["umfa_ge2"] = (d.umfa >= 2).astype(float).where(d.umfa.notna())
throws = []
for v, lab in [("umfa_p90", "UMFA > cycle 90th pct"), ("umfa_p95", "UMFA > cycle 95th pct"),
               ("umfa_ge2", "UMFA >= 2 nmol/L"), ("umfa_ge5", "UMFA >= 5 nmol/L")]:
    res, V, n, df = fit("dep", [d[v].to_numpy(float), d.l_mthf_z.to_numpy(),
                                d.l_mefox_z.to_numpy()], [lab, "5-MTHF", "MeFox"], 3)
    r = res[res.term == lab].iloc[0]
    mk = ins & (d[v] == 1).to_numpy()
    pr = des.subset(mk).prop(d.dep.to_numpy())
    mk0 = ins & (d[v] == 0).to_numpy()
    pr0 = des.subset(mk0).prop(d.dep.to_numpy())
    throws.append(dict(threshold=lab, n_exposed=int(mk.sum()),
                       prev_exposed=round(100 * pr["est"], 2),
                       prev_unexposed=round(100 * pr0["est"], 2),
                       OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                       _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
tt = pd.DataFrame(throws)
print(tt.drop(columns=[c for c in tt if c.startswith("_")]).to_string(index=False))
OUT["thresholds"] = tt.to_dict("records")

# ---------------------------------------------- 3. B12 interaction (folate trap)
print("\n" + "=" * 100)
print("VITAMIN B12 INTERACTION (2011-2014, n with B12)")
d["l_b12"] = np.log(d.b12.clip(lower=10))
b12m = ins & d.b12.notna().to_numpy()
print("  n =", int(b12m.sum()))
lowb12 = (d.b12 < 221).astype(float)      # <300 pg/mL ~ 221 pmol/L (marginal deficiency)
res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                            d.l_mefox_z.to_numpy(), d.l_b12.to_numpy(), lowb12.to_numpy(),
                            d.l_umfa_z.to_numpy() * lowb12.to_numpy()],
                    ["UMFA", "5-MTHF", "MeFox", "log B12", "Low B12", "UMFA x low B12"], 3,
                    mask=b12m)
r = res[res.term == "UMFA x low B12"].iloc[0]
print(f"  UMFA x low B12 (<221 pmol/L): OR {fmt(r.OR, r.lo, r.hi)}, p={pfmt(r.p)}")
b = res.beta.to_numpy()
for lab, v in [("normal B12", 0.0), ("low B12", 1.0)]:
    c = np.zeros(len(b))
    c[list(res.term).index("UMFA")] = 1
    c[list(res.term).index("UMFA x low B12")] = v
    est, se, lo, hi, p = contrast(b, V, c, df)
    print(f"    UMFA OR in {lab}: {np.exp(est):.2f} ({np.exp(lo):.2f}-{np.exp(hi):.2f}), p={pfmt(p)}")
    OUT.setdefault("b12", []).append(dict(stratum=lab, OR=float(np.exp(est)),
                                          lo=float(np.exp(lo)), hi=float(np.exp(hi)), p=float(p)))
OUT["b12_interaction"] = dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), p=float(r.p),
                              n=int(b12m.sum()))

# ---------------------------------------------- 4. adjusted predicted prevalence
print("\n" + "=" * 100)
print("ADJUSTED PREDICTED PREVALENCE of PHQ-9>=10 across UMFA percentiles (model 3)")
X, nm = X_model(d, 3)
X, nm = add(X, nm, [d.l_umfa.to_numpy(float), d.l_mthf_z.to_numpy()], ["logUMFA", "5-MTHF"])
res, V, n, df = svylogit(des, d.dep.to_numpy(float), X, mask=ins, names=nm)
terms = list(res.term)
bb = res.beta.to_numpy()
j = terms.index("logUMFA")
Xk = np.nan_to_num(X)[:, [nm.index(t) for t in terms]]
w = des.w * ins
pcts = [10, 25, 50, 75, 90, 95, 99]
qs = np.nanpercentile(d.l_umfa[ins], pcts)
pp = []
for q, pc in zip(qs, pcts):
    Xs = Xk.copy()
    Xs[:, j] = q
    eta = Xs @ bb
    mu = 1 / (1 + np.exp(-eta))
    est = np.sum(w * mu) / np.sum(w)
    # delta-method SE
    g = ((w * mu * (1 - mu))[:, None] * Xs).sum(axis=0) / np.sum(w)
    se = float(np.sqrt(g @ V @ g))
    pp.append(dict(pct=pc, umfa=float(np.exp(q)), prev=float(100 * est),
                   lo=float(100 * (est - 1.96 * se)), hi=float(100 * (est + 1.96 * se))))
    print(f"  p{pc:<3d} UMFA={np.exp(q):6.2f} nmol/L -> {100*est:5.2f}% "
          f"({100*(est-1.96*se):.2f}-{100*(est+1.96*se):.2f})")
OUT["pred_prev"] = pp
d99, d10 = pp[-1]["prev"], pp[0]["prev"]
print(f"  absolute difference p99 vs p10 = {d99-d10:.2f} percentage points")
OUT["pred_diff"] = float(d99 - d10)

# ---------------------------------------------- 5. total folate re-examined
print("\n" + "=" * 100)
print("WHY TOTAL FOLATE MISSES IT: total folate with and without UMFA in the model")
for lab, cols, nms in [("Total folate alone", [d.l_tfol_z.to_numpy()], ["Total folate"]),
                       ("Total folate + UMFA", [d.l_tfol_z.to_numpy(), d.l_umfa_z.to_numpy()],
                        ["Total folate", "UMFA"])]:
    res, V, n, df = fit("dep", cols, nms, 3)
    s = "  " + lab + ": " + "; ".join(
        f"{k} {fmt(res[res.term==k].iloc[0].OR, res[res.term==k].iloc[0].lo, res[res.term==k].iloc[0].hi)}"
        f" p={pfmt(res[res.term==k].iloc[0].p)}" for k in nms)
    print(s)
    OUT.setdefault("tfol_decomp", []).append(
        {k: [float(res[res.term == k].iloc[0].OR), float(res[res.term == k].iloc[0].lo),
             float(res[res.term == k].iloc[0].hi), float(res[res.term == k].iloc[0].p)] for k in nms}
        | {"model": lab})
# share of total folate that is UMFA
a = d[d.ins == 1]
print(f"  median UMFA as % of total folate = {(100*a.umfa/(a.tfol+a.mefox)).median():.2f}%")
print(f"  median 5-MTHF as % of total folate = {(100*a.mthf/(a.tfol+a.mefox)).median():.2f}%")
OUT["share"] = dict(umfa=float((100 * a.umfa / (a.tfol + a.mefox)).median()),
                    mthf=float((100 * a.mthf / (a.tfol + a.mefox)).median()))

with open("out_a3.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a3.json")
