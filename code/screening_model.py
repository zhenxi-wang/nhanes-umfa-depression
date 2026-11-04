"""Shared demographic screening model; extracted without numerical changes."""
import numpy as np
from survey import SvyDesign, svylogit

def base_X(dd, extra_cols, extra_names, adjust_sex=True):
    n = len(dd)
    cols = [np.ones(n), dd.age.to_numpy(float) / 10.0]
    names = ["(int)", "age/10"]
    if adjust_sex:
        cols.append(dd.female.to_numpy(float)); names.append("female")
    for lv in ["Mexican American", "Other Hispanic", "NH Black", "NH Asian", "Other/Multiracial"]:
        cols.append((dd.race == lv).astype(float).to_numpy()); names.append("race:" + lv)
    for lv in ["<High school", "High school/GED", "Some college"]:
        cols.append((dd.edu == lv).astype(float).to_numpy()); names.append("edu:" + lv)
    for lv in ["<1.3", "1.3-3.5"]:
        cols.append((dd.pircat.astype(str) == lv).astype(float).to_numpy()); names.append("pir:" + lv)
    pm = dd.pir.isna().astype(float).to_numpy()
    cols.append(pm); names.append("pir:miss")
    cols += extra_cols
    names += extra_names
    return np.column_stack(cols), names

def run(dd, wcol, expo, label, mask=None, adjust_sex=True, outcome="dep"):
    dd = dd.reset_index(drop=True)
    dg = dd.copy()
    dg["w"] = dg[wcol].fillna(0.0)
    des = SvyDesign(dg, "w", "SDMVSTRA", "SDMVPSU")
    x = dg[expo].to_numpy(float)
    sd = np.nanstd(x[(dg.w > 0) & (mask if mask is not None else True)])
    z = (x - np.nanmean(x)) / sd
    X, names = base_X(dg, [z], [label], adjust_sex)
    try:
        res, V, n, df = svylogit(des, dg[outcome].to_numpy(float), X, mask=mask, names=names)
    except Exception as e:
        return None
    r = res[res.term == label]
    if len(r) == 0:
        return None
    r = r.iloc[0]
    return dict(exposure=label, n=n, OR=r.OR, lo=r.lo, hi=r.hi, p=r.p)
