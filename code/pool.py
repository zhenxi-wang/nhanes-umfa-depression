"""Pool NHANES 2011-2012 .. August 2021-August 2023 for the folate-forms / depression analysis.

Cycles and folate-form weights
  G  2011-2012        2.0 y   WTMEC2YR   (folate forms on full MEC sample)
  H  2013-2014        2.0 y   WTMEC2YR
  I  2015-2016        2.0 y   WTMEC2YR
  P  2017-Mar 2020    3.2 y   WTFOLPRP   (pre-pandemic folate subsample weight)
  L  Aug2021-Aug2023  2.0 y   WTPH2YR    (phlebotomy weight)

Pooled weight = cycle weight x (cycle years / 11.2).  SDMVSTRA/SDMVPSU are already unique
across cycles (90-103, 104-118, 119-133, 149-172, 173-187), so no re-coding is needed.
"""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("NHANES_DATA_DIR", os.path.join(ROOT, "data"))

CYC = [("G", "2011-2012", 2.0, "WTMEC2YR"),
       ("H", "2013-2014", 2.0, "WTMEC2YR"),
       ("I", "2015-2016", 2.0, "WTMEC2YR"),
       ("P", "2017-Mar 2020", 3.2, "WTFOLPRP"),
       ("L", "Aug 2021-Aug 2023", 2.0, "WTPH2YR")]
TOTY = sum(c[2] for c in CYC)
MISS = (7, 9, 77, 99, 777, 999, 7777, 9999)


def fn(base, suf):
    return ("P_" + base) if suf == "P" else (base + "_" + suf)


def rd(base, suf, cols=None):
    p = os.path.join(DATA, fn(base, suf) + ".xpt")
    if not os.path.exists(p):
        return None
    d = pd.read_sas(p)
    if cols is not None:
        d = d[["SEQN"] + [c for c in cols if c in d.columns]]
    return d


def m(s, codes=MISS):
    return s.where(~s.isin(codes))


def bstr(s):
    return s.astype(str).str.replace("^b'", "", regex=True).str.replace("'$", "", regex=True).str.strip()


# --- antidepressant drug-id lookup ------------------------------------------
_dg = pd.read_sas(os.path.join(DATA, "RXQ_DRUG.xpt"))
_dg["name"] = bstr(_dg.RXDDRUG)
_cls = [c for c in _dg.columns if c.startswith("RXDDCN") or c.startswith("RXDICN")]
_txt = _dg[_cls].apply(lambda c: bstr(c)).agg(" | ".join, axis=1).str.upper()
_dg["is_ad"] = _txt.str.contains("ANTIDEPRESSANT", na=False)
_dg["is_anx"] = _txt.str.contains("ANXIOLYTIC|BENZODIAZEPINE", na=False)
_dg["is_ap"] = _txt.str.contains("ANTIPSYCHOTIC", na=False)
_dg["is_ac"] = _txt.str.contains("ANTICONVULSANT", na=False)
AD_ID = set(bstr(_dg.RXDDRGID[_dg.is_ad]))
ANX_ID = set(bstr(_dg.RXDDRGID[_dg.is_anx]))
AP_ID = set(bstr(_dg.RXDDRGID[_dg.is_ap]))


def one(suf, label, yrs, wvar):
    de = rd("DEMO", suf)
    d = pd.DataFrame({"SEQN": de.SEQN})
    d["cycle"] = label
    d["cyc"] = suf
    d["age"] = de.RIDAGEYR
    d["female"] = (de.RIAGENDR == 2).astype(float)
    rmap = {1: "Mexican American", 2: "Other Hispanic", 3: "NH White", 4: "NH Black",
            6: "NH Asian", 7: "Other/Multiracial"}
    d["race"] = de.RIDRETH3.map(rmap) if "RIDRETH3" in de else de.RIDRETH1.map(rmap)
    d["edu"] = m(de.DMDEDUC2).map({1: "<High school", 2: "<High school", 3: "High school/GED",
                                   4: "Some college", 5: "College graduate+"})
    if "DMDMARTZ" in de:
        d["mar"] = m(de.DMDMARTZ).map({1: "Married/partnered", 2: "Widowed/divorced/separated",
                                       3: "Never married"})
    else:
        d["mar"] = m(de.DMDMARTL).map({1: "Married/partnered", 6: "Married/partnered",
                                       2: "Widowed/divorced/separated", 3: "Widowed/divorced/separated",
                                       4: "Widowed/divorced/separated", 5: "Never married"})
    d["pir"] = de.INDFMPIR
    d["usborn"] = (de.DMDBORN4 == 1).astype(float).where(de.DMDBORN4.isin([1, 2]))
    d["preg"] = (de.RIDEXPRG == 1).astype(float) if "RIDEXPRG" in de else 0.0
    d["strata"] = de.SDMVSTRA
    d["psu"] = de.SDMVPSU
    d["wraw"] = de[wvar] if wvar in de else np.nan
    if wvar not in de:                      # folate weight lives in the lab file
        wf = rd("FOLFMS", suf)
        d = d.merge(wf[["SEQN", wvar]].rename(columns={wvar: "wraw2"}), on="SEQN", how="left")
        d["wraw"] = d.wraw2
        d = d.drop(columns=["wraw2"])
    d["yrs"] = yrs

    # PHQ-9
    dp = rd("DPQ", suf)
    it = [f"DPQ0{i}0" for i in range(1, 10)]
    for c in it:
        dp[c] = m(dp[c])
    dp["nm"] = dp[it].isna().sum(axis=1)
    dp["phq9"] = dp[it].sum(axis=1).where(dp.nm == 0)
    dp["dep"] = (dp.phq9 >= 10).astype(float).where(dp.phq9.notna())
    dp["dep15"] = (dp.phq9 >= 15).astype(float).where(dp.phq9.notna())
    dp["dep5"] = (dp.phq9 >= 5).astype(float).where(dp.phq9.notna())
    dp["si"] = (dp.DPQ090 >= 1).astype(float).where(dp.DPQ090.notna())
    dp["som"] = dp[["DPQ030", "DPQ040", "DPQ050", "DPQ080"]].sum(axis=1).where(dp.nm == 0)
    dp["cog"] = dp[["DPQ010", "DPQ020", "DPQ060", "DPQ070", "DPQ090"]].sum(axis=1).where(dp.nm == 0)
    d = d.merge(dp[["SEQN", "phq9", "dep", "dep15", "dep5", "si", "som", "cog"] + it],
                on="SEQN", how="left")

    # folate forms
    ff = rd("FOLFMS", suf)
    # NCHS releases a below-LOD comment flag for every species (0 = at/above LOD, 1 = below),
    # so detection is read from the flag rather than inferred from the imputed fill value,
    # whose magnitude differs between cycles.
    d = d.merge(ff[["SEQN", "LBDFOTSI", "LBXSF1SI", "LBXSF2SI", "LBXSF3SI", "LBXSF4SI",
                    "LBXSF5SI", "LBXSF6SI",
                    "LBDSF1LC", "LBDSF2LC", "LBDSF3LC", "LBDSF4LC", "LBDSF5LC", "LBDSF6LC"]],
                on="SEQN", how="left")
    d = d.rename(columns={"LBDFOTSI": "tfol", "LBXSF1SI": "mthf", "LBXSF2SI": "umfa",
                          "LBXSF3SI": "fthf", "LBXSF4SI": "thf", "LBXSF5SI": "methenyl",
                          "LBXSF6SI": "mefox", "LBDSF2LC": "umfa_bloq",
                          "LBDSF1LC": "lc_mthf", "LBDSF3LC": "lc_fthf", "LBDSF4LC": "lc_thf",
                          "LBDSF5LC": "lc_methenyl", "LBDSF6LC": "lc_mefox"})
    d["lc_umfa"] = d["umfa_bloq"]
    fo = rd("FOLATE", suf)
    if fo is not None and "LBDRFOSI" in fo:
        d = d.merge(fo[["SEQN", "LBDRFOSI"]].rename(columns={"LBDRFOSI": "rbcfol"}),
                    on="SEQN", how="left")
    else:
        d["rbcfol"] = np.nan

    # fasting / session
    fq = rd("FASTQX", suf)
    if fq is not None:
        sess = "PHDSESN" if "PHDSESN" in fq else "PHDSESNZ"
        fq2 = fq[["SEQN", "PHAFSTHR", sess]].rename(columns={sess: "session"})
        fq2["PHAFSTHR"] = fq2.PHAFSTHR.where(fq2.PHAFSTHR < 40)
        d = d.merge(fq2.rename(columns={"PHAFSTHR": "fasthr"}), on="SEQN", how="left")
    else:
        d["fasthr"] = np.nan
        d["session"] = np.nan

    # covariates
    bm = rd("BMX", suf, ["BMXBMI", "BMXWAIST"])
    d = d.merge(bm.rename(columns={"BMXBMI": "bmi", "BMXWAIST": "waist"}), on="SEQN", how="left")
    sm = rd("SMQ", suf, ["SMQ020", "SMQ040"])
    sm["smoke"] = np.where(sm.SMQ020 == 2, "Never",
                           np.where(sm.SMQ040.isin([1, 2]), "Current",
                                    np.where(sm.SMQ040 == 3, "Former", None)))
    d = d.merge(sm[["SEQN", "smoke"]], on="SEQN", how="left")
    co = rd("COT", suf, ["LBXCOT"])
    if co is None:
        co = rd("COTNAL", suf, ["LBXCOT"])
    if co is not None:
        d = d.merge(co.rename(columns={"LBXCOT": "cotinine"}), on="SEQN", how="left")
    else:
        d["cotinine"] = np.nan
    al = rd("ALQ", suf)
    if "ALQ121" in al:
        a121 = m(al.ALQ121)
        never = al.ALQ111 == 2 if "ALQ111" in al else pd.Series(False, index=al.index)
        al["drink"] = np.where(never | (a121 == 0), "None",
                               np.where(a121 <= 3, "Frequent", "Occasional"))
        al["drinks_day"] = al.ALQ130.where(al.ALQ130 < 30) if "ALQ130" in al else np.nan
    else:
        a101 = al.ALQ101 if "ALQ101" in al else pd.Series(np.nan, index=al.index)
        a120 = m(al.ALQ120Q, (777, 999)) if "ALQ120Q" in al else pd.Series(np.nan, index=al.index)
        al["drink"] = np.where(a101 == 2, "None", np.where(a120 == 0, "None", "Occasional"))
        al["drinks_day"] = al.ALQ130.where(al.ALQ130 < 30) if "ALQ130" in al else np.nan
    d = d.merge(al[["SEQN", "drink", "drinks_day"]], on="SEQN", how="left")
    pa = rd("PAQ", suf, ["PAD680", "PAQ650", "PAQ655", "PAD660", "PAQ665", "PAQ670", "PAD675",
                         "PAD790Q", "PAD790U", "PAD800", "PAD810Q", "PAD810U", "PAD820"])
    pa["sed"] = pa.PAD680.where(pa.PAD680 < 1300)
    # Recreational moderate-to-vigorous physical activity, MET-weighted minutes/week.
    # 2011-2020 used the GPAQ-style items (PAQ650/655/660 vigorous, PAQ665/670/675 moderate);
    # Aug 2021-Aug 2023 replaced them with PAD790/800 (moderate) and PAD810/820 (vigorous),
    # which record a frequency plus a unit (D/W/M/Y) rather than days per week.
    def _z(s):
        """SAS writes an exact zero as ~5.4e-79, so compare with a tolerance."""
        return s.mask(s.abs() < 1e-10, 0.0)

    def _minwk_old(yn, days, mins):
        y = m(pa[yn]) if yn in pa else pd.Series(np.nan, index=pa.index)
        dd = _z(m(pa[days], (77, 99))) if days in pa else pd.Series(np.nan, index=pa.index)
        mm = pa[mins].where(pa[mins] < 1000) if mins in pa else pd.Series(np.nan, index=pa.index)
        out = dd * mm
        return out.where(y != 2, 0.0)

    def _minwk_new(q, u, mins):
        qq = _z(pa[q].where(pa[q] < 1000))
        uu = pa[u].astype(str).str.replace("^b'", "", regex=True).str.replace("'$", "", regex=True)
        per_wk = pd.Series(np.select([uu == "D", uu == "W", uu == "M", uu == "Y"],
                                     [qq * 7.0, qq, qq / 4.345, qq / 52.18], default=np.nan),
                           index=pa.index)
        mm = pa[mins].where(pa[mins] < 1000)
        out = per_wk * mm
        return out.where(~(qq == 0), 0.0)

    if "PAQ650" in pa:
        vig = _minwk_old("PAQ650", "PAQ655", "PAD660")
        mod = _minwk_old("PAQ665", "PAQ670", "PAD675")
    else:
        vig = _minwk_new("PAD810Q", "PAD810U", "PAD820")
        mod = _minwk_new("PAD790Q", "PAD790U", "PAD800")
    pa["mvpa_min"] = (mod.fillna(0) + 2 * vig.fillna(0)).where(
        ~(mod.isna() & vig.isna())).clip(upper=5040)
    pa["mvpa"] = (pa.mvpa_min > 0).astype(float).where(pa.mvpa_min.notna())
    d = d.merge(pa[["SEQN", "sed", "mvpa", "mvpa_min"]], on="SEQN", how="left")
    sl = rd("SLQ", suf)
    scol = "SLD012" if "SLD012" in sl else "SLD010H"
    sl["sleep"] = sl[scol].where(sl[scol] < 25)
    d = d.merge(sl[["SEQN", "sleep"]], on="SEQN", how="left")
    fs = rd("FSQ", suf, ["FSDAD"])
    fs["foodinsec"] = (fs.FSDAD >= 3).astype(float).where(fs.FSDAD.notna())
    d = d.merge(fs[["SEQN", "foodinsec"]], on="SEQN", how="left")
    hi = rd("HIQ", suf, ["HIQ011"])
    hi["insured"] = m(hi.HIQ011).map({1: 1, 2: 0})
    d = d.merge(hi[["SEQN", "insured"]], on="SEQN", how="left")
    hu = rd("HUQ", suf)
    c010 = "HUQ010" if "HUQ010" in hu else "HUAQUEX"
    hu["srh"] = m(hu[c010])
    d = d.merge(hu[["SEQN", "srh"]], on="SEQN", how="left")
    mc = rd("MCQ", suf)
    for k, v in {"MCQ160B": "chf", "MCQ160C": "chd", "MCQ160F": "stroke", "MCQ160P": "copd",
                 "MCQ160L": "liver", "MCQ160A": "arthritis", "MCQ220": "cancer",
                 "MCQ160M": "thyroid"}.items():
        mc[v] = m(mc[k]).map({1: 1, 2: 0}) if k in mc else np.nan
    mc["ncond"] = mc[["chf", "chd", "stroke", "copd", "liver", "arthritis", "cancer"]].sum(axis=1)
    d = d.merge(mc[["SEQN", "chf", "chd", "stroke", "copd", "liver", "arthritis", "cancer",
                    "thyroid", "ncond"]], on="SEQN", how="left")
    di = rd("DIQ", suf, ["DIQ010"])
    di["dm"] = m(di.DIQ010).map({1: 1, 2: 0, 3: 0})
    d = d.merge(di[["SEQN", "dm"]], on="SEQN", how="left")
    bp = rd("BPQ", suf, ["BPQ020"])
    bp["htn"] = m(bp.BPQ020).map({1: 1, 2: 0})
    d = d.merge(bp[["SEQN", "htn"]], on="SEQN", how="left")
    cr = rd("HSCRP", suf, ["LBXHSCRP"])
    if cr is not None:
        d = d.merge(cr.rename(columns={"LBXHSCRP": "crp"}), on="SEQN", how="left")
    else:
        d["crp"] = np.nan
    bi = rd("BIOPRO", suf, ["LBXSCR", "LBXSAL", "LBXSGTSI"])
    d = d.merge(bi.rename(columns={"LBXSCR": "creat", "LBXSAL": "alb", "LBXSGTSI": "ggt"}),
                on="SEQN", how="left")
    cb = rd("CBC", suf, ["LBXMCVSI", "LBXHGB", "LBXRDW"])
    d = d.merge(cb.rename(columns={"LBXMCVSI": "mcv", "LBXHGB": "hgb", "LBXRDW": "rdw"}),
                on="SEQN", how="left")
    b12 = rd("VITB12", suf, ["LBDB12SI"])
    if b12 is not None and "LBDB12SI" in b12:
        d = d.merge(b12.rename(columns={"LBDB12SI": "b12"}), on="SEQN", how="left")
    else:
        d["b12"] = np.nan
    dr = rd("DR1TOT", suf, ["DR1DRSTZ", "DR1TKCAL", "DR1TFOLA", "DR1TFA", "DR1TFDFE",
                            "DR1TVB12", "DR1TALCO"])
    ok = dr.DR1DRSTZ == 1
    for c in ["DR1TKCAL", "DR1TFOLA", "DR1TFA", "DR1TFDFE", "DR1TVB12", "DR1TALCO"]:
        dr[c] = dr[c].where(ok)
    d = d.merge(dr.drop(columns=["DR1DRSTZ"]).rename(columns={
        "DR1TKCAL": "kcal", "DR1TFOLA": "d_folate", "DR1TFA": "d_fa", "DR1TFDFE": "d_dfe",
        "DR1TVB12": "d_b12", "DR1TALCO": "d_alc"}), on="SEQN", how="left")
    ds = rd("DSQTOT", suf, ["DSD010", "DSQTFA", "DSQTVB12", "DSDCOUNT"])
    ds["supp"] = m(ds.DSD010).map({1: 1, 2: 0})
    d = d.merge(ds[["SEQN", "supp", "DSQTFA", "DSQTVB12"]].rename(
        columns={"DSQTFA": "s_fa", "DSQTVB12": "s_b12"}), on="SEQN", how="left")

    # medications
    rx = rd("RXQ_RX", suf)
    if rx is not None and "RXDDRGID" in rx:
        rx["gid"] = bstr(rx.RXDDRGID)
        g = rx.groupby("SEQN").agg(
            ad=("gid", lambda s: float(any(x in AD_ID for x in s))),
            anx=("gid", lambda s: float(any(x in ANX_ID for x in s))),
            ap=("gid", lambda s: float(any(x in AP_ID for x in s))),
            nrx=("gid", lambda s: float(sum(x not in ("", "nan") for x in s))))
        d = d.merge(g.reset_index(), on="SEQN", how="left")
        use = rx.groupby("SEQN").RXDUSE.min()
        d = d.merge(use.rename("rxuse").reset_index(), on="SEQN", how="left")
        for c in ["ad", "anx", "ap"]:
            d[c] = np.where(d.rxuse == 2, 0.0, d[c])
    else:
        for c in ["ad", "anx", "ap", "nrx"]:
            d[c] = np.nan
    return d


def build():
    frames = [one(s, lab, y, w) for s, lab, y, w in CYC]
    d = pd.concat(frames, ignore_index=True)
    d["w"] = d.wraw.fillna(0.0) * d.yrs / TOTY
    return d


if __name__ == "__main__":
    d = build()
    d.to_pickle(os.path.join(DATA, "pooled.pkl"))
    print(d.shape)
    print(d.groupby("cycle", sort=False).apply(
        lambda g: pd.Series({
            "n": len(g),
            "adults18": int((g.age >= 18).sum()),
            "phq9": int(((g.age >= 18) & g.phq9.notna()).sum()),
            "folate": int(((g.age >= 18) & g.mthf.notna()).sum()),
            "both": int(((g.age >= 18) & g.phq9.notna() & g.mthf.notna() & (g.w > 0)).sum()),
            "ad": int(((g.age >= 18) & g.ad.notna()).sum()),
        }), include_groups=False))
