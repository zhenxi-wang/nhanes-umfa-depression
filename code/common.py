"""Shared analytic sample, design and design-matrix helpers."""
import os
import numpy as np
import pandas as pd
from survey import SvyDesign, svylogit, svylin, wald, contrast

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("NHANES_DATA_DIR", os.path.join(ROOT, "data"))
FIG = os.environ.get("NHANES_LEGACY_FIGURES", os.path.join(ROOT, "outputs", "legacy", "figures"))
TAB = os.path.join(ROOT, "outputs", "legacy", "tables")

CYCLES = ["2011-2012", "2013-2014", "2015-2016", "2017-Mar 2020", "Aug 2021-Aug 2023"]
RACES = ["Mexican American", "Other Hispanic", "NH Black", "NH Asian", "Other/Multiracial"]
EDU = ["<High school", "High school/GED", "Some college"]
MAR = ["Widowed/divorced/separated", "Never married"]
SMK = ["Former", "Current"]
DRK = ["None", "Frequent"]


def get():
    d = pd.read_pickle(os.path.join(DATA, "pooled.pkl")).reset_index(drop=True)

    d["ins"] = ((d.age >= 18) & (d.preg != 1) & d.phq9.notna() & d.mthf.notna() &
                d.umfa.notna() & (d.w > 0)).astype(int)

    # exposures ------------------------------------------------------------
    d["pct_umfa"] = 100 * d.umfa / (d.tfol + d.mefox)
    d["mthf_umfa"] = d.mthf / d.umfa
    for v in ["mthf", "umfa", "mefox", "tfol", "rbcfol", "pct_umfa", "thf", "fthf"]:
        d["l_" + v] = np.log(d[v].clip(lower=1e-3))
    # z-score within cycle (removes any between-cycle assay/level shift)
    a = d[d.ins == 1]
    for v in ["l_mthf", "l_umfa", "l_mefox", "l_tfol", "l_rbcfol", "l_pct_umfa"]:
        mu = a.groupby("cycle")[v].transform("mean")
        sd = a.groupby("cycle")[v].transform("std")
        mm = a.groupby("cycle")[v].mean().to_dict()
        ss = a.groupby("cycle")[v].std().to_dict()
        d[v + "_z"] = (d[v] - d.cycle.map(mm)) / d.cycle.map(ss)
    # cycle-specific quartiles for the two primary exposures
    for v in ["l_mthf", "l_umfa", "l_tfol", "l_pct_umfa"]:
        q = pd.Series(np.nan, index=d.index)
        for cy in CYCLES:
            s = a.loc[a.cycle == cy, v]
            if s.notna().sum() < 50:
                continue
            cut = np.quantile(s.dropna(), [.25, .5, .75])
            sel = d.cycle == cy
            q[sel] = np.select([d.loc[sel, v] <= cut[0], d.loc[sel, v] <= cut[1],
                                d.loc[sel, v] <= cut[2], d.loc[sel, v] > cut[2]],
                               [1, 2, 3, 4], default=np.nan)
        d[v + "_q"] = q
    # median of each quartile on the natural scale (for trend tests)
    for v in ["l_mthf", "l_umfa", "l_tfol", "l_pct_umfa"]:
        med = a.assign(_q=d.loc[a.index, v + "_q"]).groupby("_q")[v.replace("l_", "")].median()
        d[v + "_qmed"] = d[v + "_q"].map(med)

    d["umfa_det"] = (d.umfa_bloq < 0.5).astype(float).where(d.umfa_bloq.notna())
    d["egfr"] = ckdepi(d)
    d["l_crp"] = np.log(d.crp.clip(lower=0.01))
    d["l_cot"] = np.log(d.cotinine.clip(lower=0.011))
    d["l_kcal"] = np.log(d.kcal.clip(lower=100))
    d["l_dfa"] = np.log1p(d.d_fa.clip(lower=0))
    d["l_sfa"] = np.log1p(d.s_fa.fillna(0).clip(lower=0))
    d["l_dfol"] = np.log1p(d.d_folate.clip(lower=0))
    d["l_mvpa"] = np.log1p(d.mvpa_min.clip(lower=0))
    d["obese"] = (d.bmi >= 30).astype(float).where(d.bmi.notna())
    # folic-acid-containing supplement use (DSQTFA is missing for non-users)
    d["fa_supp"] = ((d.s_fa.fillna(0) > 0) | (d.l_sfa > 0)).astype(float)
    d.loc[d.supp.isna() & d.s_fa.isna(), "fa_supp"] = np.nan
    return d


def wquantile(values, weights, q):
    """Survey-weighted quantile (linear interpolation on the weighted CDF)."""
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    ok = ~np.isnan(v) & (w > 0)
    v, w = v[ok], w[ok]
    o = np.argsort(v)
    v, w = v[o], w[o]
    cw = np.cumsum(w) - 0.5 * w
    cw /= w.sum()
    return np.interp(q, cw, v)


def wmean_sd(values, weights):
    v = np.asarray(values, float)
    w = np.asarray(weights, float)
    ok = ~np.isnan(v) & (w > 0)
    v, w = v[ok], w[ok]
    mu = np.sum(w * v) / np.sum(w)
    var = np.sum(w * (v - mu) ** 2) / np.sum(w)
    return mu, np.sqrt(var)


def ckdepi(d):
    """2021 CKD-EPI creatinine equation (race-free)."""
    scr = d.creat
    k = np.where(d.female == 1, 0.7, 0.9)
    a = np.where(d.female == 1, -0.241, -0.302)
    e = 142 * np.minimum(scr / k, 1) ** a * np.maximum(scr / k, 1) ** -1.200 \
        * 0.9938 ** d.age * np.where(d.female == 1, 1.012, 1.0)
    return pd.Series(e, index=d.index).where(scr.notna())


def design(d):
    return SvyDesign(d, "w", "strata", "psu")


def dum(s, levels, prefix, ins=None):
    """Indicator coding with an explicit 'missing' category.

    Without the extra column, participants with a missing value would be silently
    absorbed into the reference level, which is not a defensible default.
    """
    cols = [np.where(s.isna(), 0.0, (s.astype(str) == lv).astype(float)) for lv in levels]
    names = [f"{prefix}: {lv}" for lv in levels]
    miss = s.isna().astype(float)
    n_miss = miss[ins].sum() if ins is not None else miss.sum()
    if n_miss > 0:
        cols.append(miss.to_numpy())
        names.append(f"{prefix}: missing")
    return cols, names


def cont(d, col, name, scale=1.0, ref=None):
    """Continuous covariate, mean-imputed with a missing indicator."""
    s = d[col].astype(float)
    mu = s[d.ins == 1].mean() if ref is None else ref
    v = (s.fillna(mu) - mu) / scale
    miss = s.isna().astype(float)
    if miss[d.ins == 1].sum() == 0:
        return [v.to_numpy()], [name]
    return [v.to_numpy(), miss.to_numpy()], [name, name + " (miss)"]


def X_model(d, level):
    """Covariate blocks.

    1 = demographic
    2 = + socioeconomic
    3 = + lifestyle / clinical / folate intake  (extended, primary)
    4 = model 3 + inflammation and MCV
    5 = DAG-guided reduced adjustment set: keeps only variables judged to be
        upstream of both circulating folate species and depressive symptoms, and drops
        those plausibly downstream of depression (BMI, sleep, physical activity, alcohol,
        chronic conditions, diabetes, marital status, insurance).
    """
    ins = (d.ins == 1)
    n = len(d)
    cols = [np.ones(n)]
    names = ["(Intercept)"]
    c, nm = cont(d, "age", "Age, per 10 y", 10.0)
    cols += c; names += nm
    cols.append(((d.age - d[ins].age.mean()) / 10.0).pow(2).to_numpy())
    names.append("Age^2")
    cols.append(d.female.to_numpy(float)); names.append("Female")
    c, nm = dum(d.race, RACES, "Race/ethnicity", ins); cols += c; names += nm
    c, nm = dum(d.cycle, CYCLES[1:], "Cycle", ins); cols += c; names += nm
    if level in (2, 3, 4, 5):
        c, nm = dum(d.edu, EDU, "Education", ins); cols += c; names += nm
        c, nm = cont(d, "pir", "Income-to-poverty ratio"); cols += c; names += nm
        c, nm = cont(d, "foodinsec", "Food insecure"); cols += c; names += nm
        c, nm = cont(d, "usborn", "US-born"); cols += c; names += nm
    if level in (2, 3, 4):
        c, nm = dum(d.mar, MAR, "Marital status", ins); cols += c; names += nm
        c, nm = cont(d, "insured", "Insured"); cols += c; names += nm
    if level in (3, 4, 5):
        c, nm = dum(d.smoke, SMK, "Smoking", ins); cols += c; names += nm
        c, nm = cont(d, "l_cot", "log cotinine"); cols += c; names += nm
        c, nm = cont(d, "egfr", "eGFR, per 10", 10.0); cols += c; names += nm
        c, nm = cont(d, "l_kcal", "log energy intake"); cols += c; names += nm
        c, nm = cont(d, "l_dfol", "log dietary folate"); cols += c; names += nm
        c, nm = cont(d, "l_dfa", "log dietary folic acid"); cols += c; names += nm
        c, nm = cont(d, "supp", "Any supplement use"); cols += c; names += nm
        c, nm = cont(d, "l_sfa", "log supplemental folic acid"); cols += c; names += nm
        c, nm = cont(d, "fasthr", "Fasting hours"); cols += c; names += nm
        c, nm = dum(d.session, ["2.0", "3.0"], "MEC session", ins); cols += c; names += nm
    if level in (3, 4):
        c, nm = dum(d.drink, DRK, "Alcohol", ins); cols += c; names += nm
        c, nm = cont(d, "bmi", "BMI, per 5", 5.0); cols += c; names += nm
        c, nm = cont(d, "sed", "Sedentary h/day", 60.0); cols += c; names += nm
        c, nm = cont(d, "l_mvpa", "log recreational MVPA (MET-min/week)"); cols += c; names += nm
        c, nm = cont(d, "sleep", "Sleep h/night"); cols += c; names += nm
        c, nm = cont(d, "ncond", "Chronic conditions"); cols += c; names += nm
        c, nm = cont(d, "dm", "Diabetes"); cols += c; names += nm
    if level == 4:
        c, nm = cont(d, "l_crp", "log hsCRP"); cols += c; names += nm
        c, nm = cont(d, "mcv", "MCV"); cols += c; names += nm
    return _dedupe(cols, names, ins.to_numpy())


def _dedupe(cols, names, ins):
    """Return a full-rank design matrix.

    Two kinds of redundancy arise here. Energy intake, dietary folate and dietary folic
    acid all come from the same 24-hour recall and share identical missingness, so their
    missingness indicators are the same column three times over. And some indicator levels
    (the evening examination session) are empty within the analytic sample. Our
    pseudo-inverse tolerates both, but an external solver need not, so we remove them at
    construction rather than relying on the estimator to cope.
    """
    keep, seen, out_names = [], {}, []
    for j, (c, n) in enumerate(zip(cols, names)):
        col = np.asarray(c, float)
        if j > 0 and np.nanstd(np.nan_to_num(col)[ins]) < 1e-12:
            continue                                   # empty within the analytic sample
        key = np.nan_to_num(col)[ins].tobytes()
        if key in seen:
            prev = out_names[seen[key]]
            if n.endswith("(miss)") and prev.endswith("(miss)"):
                out_names[seen[key]] = "24-h dietary recall (miss)"
            continue
        seen[key] = len(keep)
        keep.append(col)
        out_names.append(n)
    return np.column_stack(keep), out_names


def add(X, names, newcols, newnames):
    return np.column_stack([X] + list(newcols)), names + list(newnames)


def fmt(o, lo, hi):
    return f"{o:.2f} ({lo:.2f}-{hi:.2f})"


def pfmt(p):
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"
