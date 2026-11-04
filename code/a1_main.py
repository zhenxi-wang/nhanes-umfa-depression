"""Primary analysis: folate forms and depression."""
import json
import numpy as np
import pandas as pd
from common import (get, design, X_model, add, dum, cont, fmt, pfmt, CYCLES)
from survey import svylogit, svylin, wald, contrast

d = get()
des = design(d)
ins = (d.ins == 1).to_numpy()
OUT = {}
pd.set_option("display.width", 250)

EXPO = [("l_mthf_z", "5-MTHF (per SD log)"),
        ("l_umfa_z", "Unmetabolised folic acid (per SD log)"),
        ("l_mefox_z", "MeFox (per SD log)"),
        ("l_tfol_z", "Total serum folate (per SD log)"),
        ("l_rbcfol_z", "RBC folate (per SD log)"),
        ("l_pct_umfa_z", "UMFA as % of total folate (per SD log)")]
OUTCOMES = [("dep", "PHQ-9 >= 10"), ("dep15", "PHQ-9 >= 15"), ("si", "Suicidal ideation")]


def fit(y, extracols, extranames, level, mask=None, linear=False):
    X, nm = X_model(d, level)
    X, nm = add(X, nm, extracols, extranames)
    mk = ins if mask is None else (ins & mask)
    f = svylin if linear else svylogit
    res, V, n, df = f(des, d[y].to_numpy(float), X, mask=mk, names=nm)
    return res, V, n, df, nm


# ---------------------------------------------------------------- 1. per-SD
print("=" * 100)
print("TABLE: per-SD associations, single-exposure models")
rows = []
for ex, lab in EXPO:
    for lvl in [1, 2, 3]:
        res, V, n, df, nm = fit("dep", [d[ex].to_numpy()], [lab], lvl)
        r = res[res.term == lab].iloc[0]
        rows.append(dict(exposure=lab, model=lvl, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                         _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
t1 = pd.DataFrame(rows)
print(t1.drop(columns=[c for c in t1 if c.startswith("_")]).to_string(index=False))
OUT["persd"] = t1.to_dict("records")

# ------------------------------------------------- 2. mutually adjusted model
print("\n" + "=" * 100)
print("MUTUALLY ADJUSTED (5-MTHF + UMFA + MeFox in one model)")
mut = []
for lvl in [1, 2, 3]:
    cols = [d.l_mthf_z.to_numpy(), d.l_umfa_z.to_numpy(), d.l_mefox_z.to_numpy()]
    nms = ["5-MTHF", "UMFA", "MeFox"]
    res, V, n, df, nm = fit("dep", cols, nms, lvl)
    for k in nms:
        r = res[res.term == k].iloc[0]
        mut.append(dict(model=lvl, term=k, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                        _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
    # contrast: difference between 5-MTHF and UMFA coefficients
    i1, i2 = list(res.term).index("5-MTHF"), list(res.term).index("UMFA")
    b = res.beta.to_numpy()
    c = np.zeros(len(b)); c[i1] = 1; c[i2] = -1
    est, se, lo, hi, p = contrast(b, V, c, df)
    print(f"  model {lvl}: ratio of ORs (5-MTHF vs UMFA) = {np.exp(est):.2f} "
          f"({np.exp(lo):.2f}-{np.exp(hi):.2f}), p={pfmt(p)}")
    if lvl == 3:
        OUT["contrast_m3"] = dict(ratio=float(np.exp(est)), lo=float(np.exp(lo)),
                                  hi=float(np.exp(hi)), p=float(p))
tm = pd.DataFrame(mut)
print(tm.drop(columns=[c for c in tm if c.startswith("_")]).to_string(index=False))
OUT["mutual"] = tm.to_dict("records")

# ---------------------------------------------------------------- 3. quartiles
print("\n" + "=" * 100)
print("QUARTILES (cycle-specific), model 3, outcome PHQ-9>=10")
qrows = []
for ex, lab in [("l_mthf", "5-MTHF"), ("l_umfa", "UMFA"), ("l_tfol", "Total folate"),
                ("l_pct_umfa", "%UMFA")]:
    qv = d[ex + "_q"]
    cols = [(qv == k).astype(float).fillna(0).to_numpy() for k in [2, 3, 4]]
    nms = [f"{lab} Q{k}" for k in [2, 3, 4]]
    mk = ins & qv.notna().to_numpy()
    res, V, n, df, nm = fit("dep", cols, nms, 3, mask=mk)
    idx = [list(res.term).index(x) for x in nms]
    F, k, d2, pj = wald(res.beta.to_numpy(), V, idx, df)
    # trend on quartile median (log scale)
    trend_col = np.log(d[ex + "_qmed"].astype(float)).fillna(0).to_numpy()
    rest, Vt, nt, dft, nmt = fit("dep", [trend_col], [lab + " trend"], 3, mask=mk)
    ptr = rest[rest.term == lab + " trend"].iloc[0].p
    row = {"exposure": lab, "n": n, "Q1": "1.00 (ref)"}
    for kk, x in zip([2, 3, 4], nms):
        r = res[res.term == x].iloc[0]
        row[f"Q{kk}"] = fmt(r.OR, r.lo, r.hi)
        row[f"_q{kk}"] = [r.OR, r.lo, r.hi, r.p]
    row["p_joint"] = pfmt(pj)
    row["p_trend"] = pfmt(ptr)
    row["_p_joint"] = float(pj)
    row["_p_trend"] = float(ptr)
    qrows.append(row)
tq = pd.DataFrame(qrows)
print(tq.drop(columns=[c for c in tq if c.startswith("_")]).to_string(index=False))
OUT["quartiles"] = tq.to_dict("records")

# quartile-specific prevalence
print("\nWeighted prevalence of PHQ-9>=10 by quartile")
prev = []
for ex, lab in [("l_mthf", "5-MTHF"), ("l_umfa", "UMFA")]:
    for k in [1, 2, 3, 4]:
        mk = ins & (d[ex + "_q"] == k).to_numpy()
        p = des.subset(mk).prop(d.dep.to_numpy())
        rng = d.loc[mk, ex.replace("l_", "")]
        prev.append(dict(exposure=lab, Q=k, n=int(mk.sum()),
                         median=round(float(rng.median()), 2),
                         prev=round(100 * p["est"], 2), lo=round(100 * p["lo"], 2),
                         hi=round(100 * p["hi"], 2)))
tp = pd.DataFrame(prev)
print(tp.to_string(index=False))
OUT["prev_by_q"] = tp.to_dict("records")

# ------------------------------------------------- 4. other outcomes, model 3
print("\n" + "=" * 100)
print("OTHER OUTCOMES (model 3, mutually adjusted)")
orows = []
for y, ylab in OUTCOMES + [("dep5", "PHQ-9 >= 5")]:
    cols = [d.l_mthf_z.to_numpy(), d.l_umfa_z.to_numpy(), d.l_mefox_z.to_numpy()]
    nms = ["5-MTHF", "UMFA", "MeFox"]
    res, V, n, df, nm = fit(y, cols, nms, 3)
    row = {"outcome": ylab, "n": n}
    for k in nms:
        r = res[res.term == k].iloc[0]
        row[k] = fmt(r.OR, r.lo, r.hi)
        row["p_" + k] = pfmt(r.p)
        row["_" + k] = [r.OR, r.lo, r.hi, r.p]
    orows.append(row)
# continuous PHQ-9 and subscales
for y, ylab in [("phq9", "PHQ-9 total score"), ("som", "Somatic subscale"),
                ("cog", "Cognitive-affective subscale")]:
    cols = [d.l_mthf_z.to_numpy(), d.l_umfa_z.to_numpy(), d.l_mefox_z.to_numpy()]
    nms = ["5-MTHF", "UMFA", "MeFox"]
    res, V, n, df, nm = fit(y, cols, nms, 3, linear=True)
    row = {"outcome": ylab, "n": n}
    for k in nms:
        r = res[res.term == k].iloc[0]
        row[k] = f"{r.beta:+.3f} ({r.lo:+.3f} to {r.hi:+.3f})"
        row["p_" + k] = pfmt(r.p)
        row["_" + k] = [r.beta, r.lo, r.hi, r.p]
    orows.append(row)
to = pd.DataFrame(orows)
print(to.drop(columns=[c for c in to if c.startswith("_")]).to_string(index=False))
OUT["outcomes"] = to.to_dict("records")

# ------------------------------------------------- 5. per-cycle (consistency)
print("\n" + "=" * 100)
print("BY SURVEY CYCLE (model 2, mutually adjusted; cycle terms dropped automatically)")
crows = []
for cy in CYCLES:
    mk = ins & (d.cycle == cy).to_numpy()
    cols = [d.l_mthf_z.to_numpy(), d.l_umfa_z.to_numpy()]
    nms = ["5-MTHF", "UMFA"]
    res, V, n, df, nm = fit("dep", cols, nms, 2, mask=mk)
    row = {"cycle": cy, "n": n}
    for k in nms:
        r = res[res.term == k].iloc[0]
        row[k] = fmt(r.OR, r.lo, r.hi)
        row["_" + k] = [r.OR, r.lo, r.hi, r.p, r.beta, r.se]
    crows.append(row)
tc = pd.DataFrame(crows)
print(tc.drop(columns=[c for c in tc if c.startswith("_")]).to_string(index=False))
OUT["by_cycle"] = tc.to_dict("records")

# heterogeneity across cycles (Cochran Q on the log-OR)
for k in ["5-MTHF", "UMFA"]:
    b = np.array([r["_" + k][4] for r in tc.to_dict("records")])
    s = np.array([r["_" + k][5] for r in tc.to_dict("records")])
    wgt = 1 / s ** 2
    bbar = (wgt * b).sum() / wgt.sum()
    Q = (wgt * (b - bbar) ** 2).sum()
    from scipy import stats as st
    print(f"  {k}: pooled OR {np.exp(bbar):.2f}, Cochran Q={Q:.2f} (df=4), "
          f"p={st.chi2.sf(Q, 4):.3f}, I2={max(0, 100*(Q-4)/Q):.0f}%")
    OUT.setdefault("het", {})[k] = dict(pooled=float(np.exp(bbar)), Q=float(Q),
                                        p=float(st.chi2.sf(Q, 4)), I2=float(max(0, 100 * (Q - 4) / Q)))

with open("out_a1.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a1.json")
