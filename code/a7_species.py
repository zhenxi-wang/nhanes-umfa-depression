"""All measurable folate species: detection rates, discovery vs out-of-screen replication,
and a model containing every species simultaneously."""
import json
import numpy as np
import pandas as pd
from common import get, design, X_model, add, fmt, pfmt, CYCLES
from survey import svylogit

d = get()
des = design(d)
ins = (d.ins == 1).to_numpy()
OUT = {}
pd.set_option("display.width", 250)

SPECIES = [("mthf", "5-Methyltetrahydrofolate", "lc_mthf"),
           ("umfa", "Unmetabolised folic acid", "lc_umfa"),
           ("fthf", "5-Formyltetrahydrofolate", "lc_fthf"),
           ("thf", "Tetrahydrofolate", "lc_thf"),
           ("methenyl", "5,10-Methenyltetrahydrofolate", "lc_methenyl"),
           ("mefox", "MeFox", "lc_mefox")]

# ---- log-transform and standardise every species within cycle, as for the primary exposures
a = d[d.ins == 1]
for col, lab, _ in SPECIES:
    lv = "l_" + col
    if lv not in d:
        d[lv] = np.log(d[col].clip(lower=1e-3))
    mm = a.groupby("cycle")[lv].mean().to_dict() if lv in a else None
    if lv + "_z" not in d:
        aa = d[d.ins == 1]
        mu = aa.groupby("cycle")[lv].mean().to_dict()
        sd = aa.groupby("cycle")[lv].std().to_dict()
        d[lv + "_z"] = (d[lv] - d.cycle.map(mu)) / d.cycle.map(sd)


def fit(cols, names, mask=None, level=3):
    X, nm = X_model(d, level)
    X, nm = add(X, nm, cols, names)
    mk = ins if mask is None else (ins & mask)
    return svylogit(des, d.dep.to_numpy(float), X, mask=mk, names=nm)


def one(col, lab, mask):
    r, V, n, df = fit([d["l_" + col + "_z"].to_numpy()], [lab], mask=mask)
    x = r[r.term == lab].iloc[0]
    return dict(OR=fmt(x.OR, x.lo, x.hi), p=pfmt(x.p), n=int(n),
                _or=float(x.OR), _lo=float(x.lo), _hi=float(x.hi), _p=float(x.p))


# ---------------------------------------------------------------- detection rates
print("=" * 100)
print("DETECTION RATES AND DISTRIBUTION OF EACH FOLATE SPECIES")
rows = []
disc = ins & (d.cyc == "L").to_numpy()
repl = ins & (d.cyc != "L").to_numpy()
for col, lab, lc in SPECIES:
    sub = d.loc[ins, col]
    # NCHS comment flag: 0 = at/above LOD, 1 = below (released for all six species).
    # SAS writes the zero as ~5.4e-79, so compare with a tolerance.
    flag = d.loc[ins, lc]
    det = f"{100 * (flag < 0.5).mean():.1f}"
    dsc = one(col, lab, disc)
    rep = one(col, lab, repl)
    poo = one(col, lab, None)
    rows.append({"Folate species": lab, "Above LOD, %": det,
                 "Median, nmol/L": f"{sub.median():.2f}",
                 "Discovery OR (95% CI)": dsc["OR"], "p (discovery)": dsc["p"],
                 "Replication OR (95% CI)": rep["OR"], "p (replication)": rep["p"],
                 "Pooled OR (95% CI)": poo["OR"], "p (pooled)": poo["p"]})
t = pd.DataFrame(rows)
print(t.to_string(index=False))
OUT["species"] = t.to_dict("records")

# ---------------------------------------------------------------- all species together
print("\n" + "=" * 100)
print("ALL SIX FOLATE SPECIES MUTUALLY ADJUSTED (Model 3)")
cols = [d["l_" + c + "_z"].to_numpy() for c, _, _ in SPECIES]
nms = [lab for _, lab, _ in SPECIES]
for tag, mask in [("Pooled 2011-2023", None), ("Out-of-screen 2011-Mar 2020", repl),
                  ("Discovery Aug 2021-Aug 2023", disc)]:
    r, V, n, df = fit(cols, nms, mask=mask)
    print(f"  {tag} (n={n})")
    sub = []
    for lab in nms:
        x = r[r.term == lab].iloc[0]
        print(f"     {lab:32s} {fmt(x.OR, x.lo, x.hi)}  p={pfmt(x.p)}")
        sub.append(dict(species=lab, OR=fmt(x.OR, x.lo, x.hi), p=pfmt(x.p)))
    OUT.setdefault("all_species", {})[tag] = dict(n=int(n), rows=sub)

with open("out_a7.json", "w") as f:
    json.dump(OUT, f, indent=1, default=float)
print("\nsaved out_a7.json")
