"""Systematic screen of candidate exposures against PHQ-9 depression."""
import os
import numpy as np
import pandas as pd
from screening_model import run

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("NHANES_DATA_DIR", os.path.join(ROOT, "data"))
d = pd.read_pickle(os.path.join(DATA, "wide.pkl")).reset_index(drop=True)

d["logamh"] = np.log(d.amh.clip(lower=0.015))
for v in ["fsh", "lh", "prog", "e2", "tst", "shbg", "dheas", "andro", "ohp17", "e1", "e1s",
          "crp", "fol_5mthf", "fol_fa", "fol_mefox", "fol_tot", "rbc_folate", "lsm", "cap",
          "bche", "vitd"]:
    d["log_" + v] = np.log(d[v].clip(lower=d[v][d[v] > 0].min() if (d[v] > 0).any() else 0.01))






adult = (d.age >= 18) & d.phq9.notna()
out = []

# ---------- hormones, by sex and life stage --------------------------------
horm = ["log_fsh", "logamh", "log_lh", "log_e2", "log_prog", "log_tst", "log_shbg",
        "log_dheas", "log_andro", "log_ohp17", "log_e1", "log_e1s"]
groups = {
    "women 18+": adult & (d.female == 1),
    "women 40-60": adult & (d.female == 1) & d.age.between(40, 60),
    "women 45-64": adult & (d.female == 1) & d.age.between(45, 64),
    "women 18-44": adult & (d.female == 1) & d.age.between(18, 44),
    "men 18+": adult & (d.female == 0),
    "men 40-70": adult & (d.female == 0) & d.age.between(40, 70),
}
for gname, gmask in groups.items():
    for h in horm:
        m = (gmask & d[h].notna()).to_numpy()
        if m.sum() < 250:
            continue
        r = run(d, "wph", h, h, mask=m, adjust_sex=False)
        if r:
            r["group"] = gname
            out.append(r)

# ---------- other biomarker / exam exposures --------------------------------
others = {
    "log_fol_5mthf": adult, "log_fol_fa": adult, "log_fol_mefox": adult,
    "log_rbc_folate": adult, "log_bche": adult, "log_lsm": adult, "log_cap": adult,
    "log_vitd": adult,
}
for h, gm in others.items():
    m = (gm & d[h].notna()).to_numpy()
    r = run(d, "wph", h, h, mask=m)
    if r:
        r["group"] = "adults 18+"
        out.append(r)

# binary/ordinal exposures with MEC weight
bin_exp = ["romberg_fail4", "romberg_fail3", "dizzy_any", "vertigo", "imbalance",
           "fall12", "wg_see_any", "wg_hear_any", "wg_walk_any", "wg_comm_any",
           "wg_cog_any", "wg_self_any", "wg_lift_any", "wg_hand_any", "foodinsec",
           "tinnitus"]
for b in bin_exp:
    if b not in d:
        continue
    m = (adult & d[b].notna()).to_numpy()
    if m.sum() < 250:
        continue
    r = run(d, "wmec", b, b, mask=m)
    if r:
        r["group"] = "adults 18+ (binary, per SD)"
        out.append(r)

res = pd.DataFrame(out).sort_values("p")
pd.set_option("display.width", 200)
print(res.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
res.to_csv("screen.csv", index=False)
