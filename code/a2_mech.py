"""Mechanistic extension: dose-response, folate-conversion phenotype, effect modification,
sensitivity analyses."""
import json
import numpy as np
import pandas as pd
from scipy import stats as st
from common import get, design, X_model, add, fmt, pfmt, CYCLES
from survey import svylogit, wald, contrast

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


def rcs(x, knots):
    """Restricted cubic spline basis (Harrell), returns k-2 columns."""
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


# ------------------------------------------------------ 1. spline dose-response
print("=" * 100)
print("SPLINE: log-UMFA vs PHQ-9>=10 (model 3), 4 knots at 5/35/65/95th pct")
x = d.l_umfa.to_numpy(float)
kn = np.nanpercentile(x[ins], [5, 35, 65, 95])
sp = rcs(x, kn)
cols = ([np.nan_to_num(x - kn[0])] + [np.nan_to_num(s) for s in sp]
        + [d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()])
nms = ["umfa_lin", "umfa_s1", "umfa_s2", "5-MTHF", "MeFox"]
res, V, n, df = fit("dep", cols, nms, 3)
idx = [list(res.term).index(v) for v in ["umfa_s1", "umfa_s2"]]
F, k_, d2, p_nl = wald(res.beta.to_numpy(), V, idx, df)
print(f"  non-linearity p = {pfmt(p_nl)}  (n={n})")
b = res.beta.to_numpy()
ii = [list(res.term).index(v) for v in nms[:3]]
ref = np.nanpercentile(x[ins], 10)
grid = np.linspace(np.nanpercentile(x[ins], 2), np.nanpercentile(x[ins], 98), 60)
curve = []
for g in grid:
    cg = np.zeros(len(b))
    bas_g = [g - kn[0]] + [s[0] for s in rcs(np.array([g]), kn)]
    bas_r = [ref - kn[0]] + [s[0] for s in rcs(np.array([ref]), kn)]
    for j, i0 in enumerate(ii):
        cg[i0] = bas_g[j] - bas_r[j]
    est, se, lo, hi, p = contrast(b, V, cg, df)
    curve.append(dict(umfa=float(np.exp(g)), OR=float(np.exp(est)),
                      lo=float(np.exp(lo)), hi=float(np.exp(hi))))
OUT["spline"] = dict(p_nonlinear=float(p_nl), ref=float(np.exp(ref)), curve=curve)
for q in [0.5, 1.0, 2.0, 4.0, 8.0]:
    c = min(curve, key=lambda r: abs(r["umfa"] - q))
    print(f"   UMFA={c['umfa']:.2f} nmol/L  OR vs {np.exp(ref):.2f} = "
          f"{c['OR']:.2f} ({c['lo']:.2f}-{c['hi']:.2f})")

# --------------------------------------- 2. folate-conversion (discordance) phenotype
print("\n" + "=" * 100)
print("FOLATE-CONVERSION PHENOTYPE: cycle-specific median splits of 5-MTHF and UMFA")
a = d[d.ins == 1]
hi_m = pd.Series(np.nan, index=d.index)
hi_u = pd.Series(np.nan, index=d.index)
for cy in CYCLES:
    sel = (d.cycle == cy)
    mm = a.loc[a.cycle == cy, "l_mthf"].median()
    uu = a.loc[a.cycle == cy, "l_umfa"].median()
    hi_m[sel] = (d.loc[sel, "l_mthf"] > mm).astype(float)
    hi_u[sel] = (d.loc[sel, "l_umfa"] > uu).astype(float)
d["hi_mthf"] = hi_m
d["hi_umfa"] = hi_u
grp = np.select(
    [(d.hi_mthf == 1) & (d.hi_umfa == 0), (d.hi_mthf == 1) & (d.hi_umfa == 1),
     (d.hi_mthf == 0) & (d.hi_umfa == 0), (d.hi_mthf == 0) & (d.hi_umfa == 1)],
    ["A: high 5-MTHF / low UMFA", "B: high 5-MTHF / high UMFA",
     "C: low 5-MTHF / low UMFA", "D: low 5-MTHF / high UMFA"], default=None)
d["fgrp"] = pd.Series(grp).replace("None", np.nan)
GR = ["B: high 5-MTHF / high UMFA", "C: low 5-MTHF / low UMFA", "D: low 5-MTHF / high UMFA"]
cols = [(d.fgrp == g).astype(float).to_numpy() for g in GR]
prev_rows = []
for g in ["A: high 5-MTHF / low UMFA"] + GR:
    mk = ins & (d.fgrp == g).to_numpy()
    p = des.subset(mk).prop(d.dep.to_numpy())
    prev_rows.append(dict(group=g, n=int(mk.sum()), prev=round(100 * p["est"], 2),
                          lo=round(100 * p["lo"], 2), hi=round(100 * p["hi"], 2)))
print(pd.DataFrame(prev_rows).to_string(index=False))
grows = []
for lvl in [1, 2, 3]:
    res, V, n, df = fit("dep", cols, GR, lvl)
    row = {"model": lvl, "n": n, "A: high 5-MTHF / low UMFA": "1.00 (ref)"}
    for g in GR:
        r = res[res.term == g].iloc[0]
        row[g] = fmt(r.OR, r.lo, r.hi)
        row["_" + g] = [r.OR, r.lo, r.hi, r.p]
    grows.append(row)
tg = pd.DataFrame(grows)
print(tg.drop(columns=[c for c in tg if c.startswith("_")]).to_string(index=False))
OUT["phenotype"] = dict(prev=prev_rows, models=tg.to_dict("records"))

# multiplicative interaction (continuous)
inter = (d.l_umfa_z * d.l_mthf_z).to_numpy()
res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                            d.l_mefox_z.to_numpy(), inter],
                    ["UMFA", "5-MTHF", "MeFox", "UMFA x 5-MTHF"], 3)
r = res[res.term == "UMFA x 5-MTHF"].iloc[0]
print(f"\n  continuous interaction UMFA x 5-MTHF: OR {fmt(r.OR, r.lo, r.hi)}, p={pfmt(r.p)}")
OUT["interaction_cont"] = dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), p=float(r.p))
# UMFA effect at low vs high 5-MTHF (-1 SD / +1 SD)
b = res.beta.to_numpy()
for lab, m in [("5-MTHF -1 SD", -1.0), ("5-MTHF mean", 0.0), ("5-MTHF +1 SD", 1.0)]:
    c = np.zeros(len(b))
    c[list(res.term).index("UMFA")] = 1
    c[list(res.term).index("UMFA x 5-MTHF")] = m
    est, se, lo, hi, p = contrast(b, V, c, df)
    print(f"    UMFA OR at {lab}: {np.exp(est):.2f} ({np.exp(lo):.2f}-{np.exp(hi):.2f}), p={pfmt(p)}")
    OUT.setdefault("umfa_by_mthf", []).append(
        dict(at=lab, OR=float(np.exp(est)), lo=float(np.exp(lo)), hi=float(np.exp(hi)), p=float(p)))

# ------------------------------------------------------ 3. effect modification
print("\n" + "=" * 100)
print("SUBGROUPS: UMFA per SD (model 3, adjusted for 5-MTHF)")
subs = {
    "Age 18-39": ins & d.age.between(18, 39).to_numpy(),
    "Age 40-59": ins & d.age.between(40, 59).to_numpy(),
    "Age >=60": ins & (d.age >= 60).to_numpy(),
    "Men": ins & (d.female == 0).to_numpy(),
    "Women": ins & (d.female == 1).to_numpy(),
    "Supplement user": ins & (d.supp == 1).to_numpy(),
    "No supplements": ins & (d.supp == 0).to_numpy(),
    "Never smoker": ins & (d.smoke == "Never").to_numpy(),
    "Ever smoker": ins & d.smoke.isin(["Former", "Current"]).to_numpy(),
    "eGFR >=60": ins & (d.egfr >= 60).to_numpy(),
    "eGFR <60": ins & (d.egfr < 60).to_numpy(),
    "NH White": ins & (d.race == "NH White").to_numpy(),
    "Other race/ethnicity": ins & (d.race != "NH White").to_numpy(),
    "Fasted >=8 h": ins & (d.fasthr >= 8).to_numpy(),
    "Fasted <8 h": ins & (d.fasthr < 8).to_numpy(),
}
srows = []
for lab, mk in subs.items():
    if mk.sum() < 400:
        continue
    res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                d.l_mefox_z.to_numpy()],
                        ["UMFA", "5-MTHF", "MeFox"], 3, mask=mk)
    r = res[res.term == "UMFA"].iloc[0]
    srows.append(dict(subgroup=lab, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                      _or=r.OR, _lo=r.lo, _hi=r.hi, _b=r.beta, _se=r.se, _p=float(r.p)))
ts = pd.DataFrame(srows)
print(ts.drop(columns=[c for c in ts if c.startswith("_")]).to_string(index=False))
OUT["subgroups"] = ts.to_dict("records")

# formal interaction tests
print("\n  formal interaction tests (model 3):")
itests = {}
for lab, col in [("age>=60", (d.age >= 60).astype(float)),
                 ("female", d.female),
                 ("supplement user", d.supp.fillna(0)),
                 ("ever smoker", d.smoke.isin(["Former", "Current"]).astype(float))]:
    cc = col.to_numpy(float)
    res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                d.l_mefox_z.to_numpy(), cc,
                                d.l_umfa_z.to_numpy() * cc],
                        ["UMFA", "5-MTHF", "MeFox", lab, "UMFA x " + lab], 3)
    r = res[res.term == "UMFA x " + lab].iloc[0]
    print(f"    UMFA x {lab}: OR {fmt(r.OR, r.lo, r.hi)}, p={pfmt(r.p)}")
    itests[lab] = dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), p=float(r.p))
OUT["interaction_tests"] = itests

# ------------------------------------------------------ 4. sensitivity
print("\n" + "=" * 100)
print("SENSITIVITY (UMFA per SD, adjusted for 5-MTHF)")
sens = {
    "Main (model 3)": (ins, 3),
    "Model 4 (+hsCRP, MCV; 2015-2023)": (ins & d.crp.notna().to_numpy(), 4),
    "Excluding antidepressant users": (ins & (d.ad == 0).to_numpy(), 3),
    "Excluding all psychotropic users": (ins & (d.ad == 0).to_numpy() &
                                         (d.anx == 0).to_numpy() & (d.ap == 0).to_numpy(), 3),
    "Excluding eGFR <60": (ins & (d.egfr >= 60).to_numpy(), 3),
    "Excluding MCV >=100 fL": (ins & (d.mcv < 100).to_numpy(), 3),
    "Excluding severe obesity (BMI>=40)": (ins & (d.bmi < 40).to_numpy(), 3),
    "Complete covariate data only": (ins & d[["edu", "pir", "mar", "bmi", "smoke", "sleep",
                                              "kcal", "supp"]].notna().all(axis=1).to_numpy(), 3),
    "Excluding 2011-2012": (ins & (d.cycle != "2011-2012").to_numpy(), 3),
    "Aug 2021-Aug 2023 only": (ins & (d.cycle == "Aug 2021-Aug 2023").to_numpy(), 3),
    "2011-2020 only": (ins & (d.cycle != "Aug 2021-Aug 2023").to_numpy(), 3),
}
srows = []
for lab, (mk, lvl) in sens.items():
    res, V, n, df = fit("dep", [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                d.l_mefox_z.to_numpy()],
                        ["UMFA", "5-MTHF", "MeFox"], lvl, mask=mk)
    r = res[res.term == "UMFA"].iloc[0]
    rm = res[res.term == "5-MTHF"].iloc[0]
    srows.append(dict(analysis=lab, n=n, UMFA=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                      MTHF=fmt(rm.OR, rm.lo, rm.hi), p_mthf=pfmt(rm.p),
                      _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
tsn = pd.DataFrame(srows)
print(tsn.drop(columns=[c for c in tsn if c.startswith("_")]).to_string(index=False))
OUT["sensitivity"] = tsn.to_dict("records")

# ------------------------------------------------------ 5. E-value
def evalue(orv):
    rr = np.sqrt(orv) if orv > 1 else 1 / np.sqrt(1 / orv)
    return rr + np.sqrt(rr * (rr - 1))


m = [r for r in OUT["sensitivity"] if r["analysis"] == "Main (model 3)"][0]
print(f"\nE-value for main OR {m['_or']:.2f}: {evalue(m['_or']):.2f} "
      f"(CI limit {m['_lo']:.2f}: {evalue(max(m['_lo'],1.0001)):.2f})")
OUT["evalue"] = dict(point=float(evalue(m["_or"])), ci=float(evalue(max(m["_lo"], 1.0001))))

with open("out_a2.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a2.json")
