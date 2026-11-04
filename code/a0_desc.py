"""Descriptives: cross-cycle comparability of the folate forms and outcome prevalence."""
import numpy as np
import pandas as pd
from common import get, design, CYCLES

d = get()
a = d[d.ins == 1]
print("analytic n =", len(a))
print(a.groupby("cycle", sort=False).size())

des = design(d)
ins = (d.ins == 1).to_numpy()

print("\n--- weighted geometric means by cycle (nmol/L) ---")
rows = []
for cy in CYCLES:
    mk = ins & (d.cycle == cy).to_numpy()
    r = {"cycle": cy, "n": int(mk.sum())}
    for v, lab in [("l_mthf", "5-MTHF"), ("l_umfa", "UMFA"), ("l_mefox", "MeFox"),
                   ("l_tfol", "Total folate"), ("l_rbcfol", "RBC folate"),
                   ("l_pct_umfa", "%UMFA")]:
        s = des.subset(mk).mean(d[v].to_numpy())
        r[lab] = f"{np.exp(s[0]):.2f} ({np.exp(s[2]):.2f}-{np.exp(s[3]):.2f})"
    r["UMFA detected %"] = f"{100*des.subset(mk).mean(d.umfa_det.to_numpy())[0]:.1f}"
    r["PHQ9>=10 %"] = f"{100*des.subset(mk).mean(d.dep.to_numpy())[0]:.1f}"
    r["SI %"] = f"{100*des.subset(mk).mean(d.si.to_numpy())[0]:.1f}"
    rows.append(r)
t = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(t.to_string(index=False))

print("\n--- overall ---")
for v, lab in [("l_mthf", "5-MTHF"), ("l_umfa", "UMFA"), ("l_mefox", "MeFox"),
               ("l_tfol", "Total folate"), ("l_rbcfol", "RBC folate"), ("l_pct_umfa", "%UMFA")]:
    s = des.subset(ins).mean(d[v].to_numpy())
    q = np.nanpercentile(np.exp(d.loc[ins, v]), [5, 25, 50, 75, 95])
    print(f"{lab:14s} GM {np.exp(s[0]):8.2f} ({np.exp(s[2]):.2f}-{np.exp(s[3]):.2f})   "
          f"p5/25/50/75/95 = {q[0]:.2f}/{q[1]:.2f}/{q[2]:.2f}/{q[3]:.2f}/{q[4]:.2f}")
print("UMFA above LOD: %.1f%%" % (100 * des.subset(ins).mean(d.umfa_det.to_numpy())[0]))
print("PHQ-9>=10: %.2f%%" % (100 * des.subset(ins).mean(d.dep.to_numpy())[0]))
print("PHQ-9>=15: %.2f%%" % (100 * des.subset(ins).mean(d.dep15.to_numpy())[0]))
print("SI: %.2f%%" % (100 * des.subset(ins).mean(d.si.to_numpy())[0]))

print("\n--- correlations among forms (unweighted, log scale, analytic sample) ---")
print(a[["l_mthf", "l_umfa", "l_mefox", "l_tfol", "l_rbcfol", "l_pct_umfa"]].corr().round(2))

print("\n--- covariate missingness in analytic sample (%) ---")
cov = ["edu", "pir", "mar", "insured", "foodinsec", "usborn", "smoke", "cotinine", "drink",
       "bmi", "sed", "mvpa", "sleep", "ncond", "dm", "egfr", "kcal", "d_folate", "d_fa",
       "supp", "s_fa", "fasthr", "session", "crp", "mcv", "ad", "rbcfol"]
print((a[cov].isna().mean() * 100).round(1).to_string())
