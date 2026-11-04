"""Table 1 (characteristics by UMFA quartile), plus remaining robustness checks."""
import json
import numpy as np
import pandas as pd
from common import get, design, X_model, add, fmt, pfmt, wmean_sd, CYCLES, RACES, EDU, MAR
from survey import svylogit, wald

d = get()
des = design(d)
ins = (d.ins == 1).to_numpy()
OUT = {}
pd.set_option("display.width", 260)

d["s_fa_z0"] = d.s_fa.fillna(0.0)
QM = {1: "Q1 (lowest)", 2: "Q2", 3: "Q3", 4: "Q4 (highest)"}
groups = [("Overall", ins)] + [(QM[k], ins & (d.l_umfa_q == k).to_numpy()) for k in [1, 2, 3, 4]]


def pct(var, val=1):
    out = {}
    y = (d[var] == val).astype(float).where(d[var].notna()).to_numpy()
    for lab, mk in groups:
        p = des.subset(mk).prop(y)
        out[lab] = f"{100*p['est']:.1f}"
    return out


def cat(var, level):
    out = {}
    y = (d[var].astype(str) == level).astype(float).where(d[var].notna()).to_numpy()
    for lab, mk in groups:
        p = des.subset(mk).prop(y)
        out[lab] = f"{100*p['est']:.1f}"
    return out


def mean(var, dec=1):
    out = {}
    for lab, mk in groups:
        mu, se, lo, hi = des.subset(mk).mean(d[var].to_numpy(float))
        out[lab] = f"{mu:.{dec}f}"
    return out


def gm(var, dec=2):
    out = {}
    for lab, mk in groups:
        mu, se, lo, hi = des.subset(mk).mean(np.log(d[var].clip(lower=1e-3)).to_numpy(float))
        out[lab] = f"{np.exp(mu):.{dec}f}"
    return out


rows = []
rows.append({"Characteristic": "Unweighted n", **{lab: str(int(mk.sum())) for lab, mk in groups}})
rows.append({"Characteristic": "UMFA, nmol/L, median (IQR)",
             **{lab: (f"{d.umfa[mk].median():.2f} "
                      f"({d.umfa[mk].quantile(.25):.2f}-{d.umfa[mk].quantile(.75):.2f})")
                for lab, mk in groups}})
rows.append({"Characteristic": "Age, years, mean", **mean("age")})
rows.append({"Characteristic": "Female, %", **pct("female")})
for r in RACES + ["NH White"]:
    rows.append({"Characteristic": f"  {r}, %", **cat("race", r)})
for e in EDU + ["College graduate+"]:
    rows.append({"Characteristic": f"  {e}, %", **cat("edu", e)})
rows.append({"Characteristic": "Income-to-poverty ratio, mean", **mean("pir", 2)})
for mm in MAR + ["Married/partnered"]:
    rows.append({"Characteristic": f"  {mm}, %", **cat("mar", mm)})
rows.append({"Characteristic": "US-born, %", **pct("usborn")})
rows.append({"Characteristic": "Insured, %", **pct("insured")})
rows.append({"Characteristic": "Food insecure, %", **pct("foodinsec")})
for s in ["Never", "Former", "Current"]:
    rows.append({"Characteristic": f"  Smoking: {s}, %", **cat("smoke", s)})
rows.append({"Characteristic": "BMI, kg/m2, mean", **mean("bmi")})
rows.append({"Characteristic": "Obesity (BMI>=30), %", **pct("obese")})
rows.append({"Characteristic": "Sleep, h/night, mean", **mean("sleep")})
rows.append({"Characteristic": "Chronic conditions, mean", **mean("ncond", 2)})
rows.append({"Characteristic": "Diabetes, %", **pct("dm")})
rows.append({"Characteristic": "eGFR, mL/min/1.73m2, mean", **mean("egfr")})
rows.append({"Characteristic": "Any dietary supplement, %", **pct("supp")})
rows.append({"Characteristic": "Folic-acid-containing supplement, %", **pct("fa_supp")})
rows.append({"Characteristic": "Supplemental folic acid, ug/d, mean", **mean("s_fa_z0", 1)})
rows.append({"Characteristic": "Dietary folic acid, ug/d, GM", **gm("d_fa", 1)})
rows.append({"Characteristic": "Energy intake, kcal/d, mean", **mean("kcal", 0)})
rows.append({"Characteristic": "Fasting, hours, mean", **mean("fasthr")})
rows.append({"Characteristic": "5-MTHF, nmol/L, GM", **gm("mthf", 1)})
rows.append({"Characteristic": "Total serum folate, nmol/L, GM", **gm("tfol", 1)})
rows.append({"Characteristic": "RBC folate, nmol/L, GM", **gm("rbcfol", 0)})
rows.append({"Characteristic": "MeFox, nmol/L, GM", **gm("mefox", 2)})
rows.append({"Characteristic": "PHQ-9 score, mean", **mean("phq9", 2)})
rows.append({"Characteristic": "PHQ-9 >=10, %", **pct("dep")})
rows.append({"Characteristic": "PHQ-9 >=15, %", **pct("dep15")})
rows.append({"Characteristic": "Suicidal ideation, %", **pct("si")})
t1 = pd.DataFrame(rows)
print(t1.to_string(index=False))
t1.to_csv("tab1.csv", index=False)
OUT["table1"] = t1.to_dict("records")

# --------------------------------------------------------- extra robustness
print("\n" + "=" * 100)
print("ADDITIONAL ROBUSTNESS (UMFA per SD, model 3, adjusted for 5-MTHF)")


def fit(y, cols, names, level, mask=None):
    X, nm = X_model(d, level)
    X, nm = add(X, nm, cols, names)
    mk = ins if mask is None else (ins & mask)
    return svylogit(des, d[y].to_numpy(float), X, mask=mk, names=nm)


p99 = d.l_umfa[ins].quantile(0.99)
extra = {
    "Main (model 3)": (ins, [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()], ["UMFA", "5-MTHF", "MeFox"]),
    "Excluding UMFA >99th percentile": (ins & (d.l_umfa <= p99).to_numpy(),
                                        [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                         d.l_mefox_z.to_numpy()],
                                        ["UMFA", "5-MTHF", "MeFox"]),
    "Additionally adjusting for RBC folate": (ins & d.rbcfol.notna().to_numpy(),
                                              [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                               d.l_mefox_z.to_numpy(), d.l_rbcfol_z.to_numpy()],
                                              ["UMFA", "5-MTHF", "MeFox", "RBC folate"]),
    "Additionally adjusting for self-rated health": (
        ins, [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy(),
              d.srh.fillna(3).to_numpy()],
        ["UMFA", "5-MTHF", "MeFox", "Self-rated health"]),
    "Age >=20 only": (ins & (d.age >= 20).to_numpy(),
                      [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()], ["UMFA", "5-MTHF", "MeFox"]),
    "Non-Hispanic White only": (ins & (d.race == "NH White").to_numpy(),
                                [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(),
                                 d.l_mefox_z.to_numpy()],
                                ["UMFA", "5-MTHF", "MeFox"]),
}
rows = []
for lab, (mk, cols, nms) in extra.items():
    res, V, n, df = fit("dep", cols, nms, 3, mask=mk)
    r = res[res.term == "UMFA"].iloc[0]
    rows.append(dict(analysis=lab, n=n, OR=fmt(r.OR, r.lo, r.hi), p=pfmt(r.p),
                     _or=r.OR, _lo=r.lo, _hi=r.hi, _p=r.p))
te = pd.DataFrame(rows)
print(te.drop(columns=[c for c in te if c.startswith("_")]).to_string(index=False))
OUT["robust_extra"] = te.to_dict("records")

# ordinal quartile score trend
print("\nQuartile trend using ordinal score (1-4):")
for ex, lab in [("l_umfa", "UMFA"), ("l_pct_umfa", "%UMFA"), ("l_mthf", "5-MTHF"),
                ("l_tfol", "Total folate")]:
    sc = d[ex + "_q"].astype(float).fillna(0).to_numpy()
    mk = ins & d[ex + "_q"].notna().to_numpy()
    res, V, n, df = fit("dep", [sc], [lab + " score"], 3, mask=mk)
    r = res[res.term == lab + " score"].iloc[0]
    print(f"  {lab:14s} OR per quartile {fmt(r.OR, r.lo, r.hi)}  p_trend={pfmt(r.p)}")
    OUT.setdefault("ordinal_trend", {})[lab] = dict(OR=float(r.OR), lo=float(r.lo),
                                                    hi=float(r.hi), p=float(r.p))

# unweighted comparison (unit weights, same design-based variance machinery)
print("\nUnweighted logistic (same covariates, model 3):")
du = d.copy()
du["w"] = ins.astype(float)
desu = design(du)
X, nm = X_model(d, 3)
X, nm = add(X, nm, [d.l_umfa_z.to_numpy(), d.l_mthf_z.to_numpy(), d.l_mefox_z.to_numpy()], ["UMFA", "5-MTHF", "MeFox"])
res, V, n, df = svylogit(desu, d.dep.to_numpy(float), X, mask=ins, names=nm)
r = res[res.term == "UMFA"].iloc[0]
print(f"  UMFA OR {fmt(r.OR, r.lo, r.hi)}, p={pfmt(r.p)}  (n={n})")
OUT["unweighted"] = dict(OR=float(r.OR), lo=float(r.lo), hi=float(r.hi), p=float(r.p))

# standardized differences Q1 vs Q4 (avoids significance testing of a baseline table)
print("\nStandardised differences, Q1 vs Q4")
sdrows = []
q1 = ins & (d.l_umfa_q == 1).to_numpy()
q4 = ins & (d.l_umfa_q == 4).to_numpy()
for lab, var in [("Age", "age"), ("Female", "female"), ("Income-to-poverty ratio", "pir"),
                 ("Food insecure", "foodinsec"), ("BMI", "bmi"), ("Sleep", "sleep"),
                 ("Chronic conditions", "ncond"), ("Diabetes", "dm"), ("eGFR", "egfr"),
                 ("Any supplement", "supp"), ("Folic-acid supplement", "fa_supp"),
                 ("Supplemental folic acid", "s_fa_z0"), ("Dietary folic acid", "d_fa"),
                 ("Energy intake", "kcal"), ("Fasting hours", "fasthr"),
                 ("5-MTHF", "mthf"), ("PHQ-9 score", "phq9")]:
    m1, s1 = wmean_sd(d.loc[q1, var].astype(float), d.loc[q1, "w"])
    m4, s4 = wmean_sd(d.loc[q4, var].astype(float), d.loc[q4, "w"])
    sp = np.sqrt((s1 ** 2 + s4 ** 2) / 2)
    smd = (m4 - m1) / sp if sp > 0 else np.nan
    sdrows.append(dict(Variable=lab, Q1=round(float(m1), 2), Q4=round(float(m4), 2),
                       SMD=round(float(smd), 3)))
tsd = pd.DataFrame(sdrows)
print(tsd.to_string(index=False))
OUT["smd"] = tsd.to_dict("records")

with open("out_a4.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a4.json")
