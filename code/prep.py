"""Build a wide NHANES August 2021-August 2023 analytic frame for screening."""
import os
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.environ.get("NHANES_DATA_DIR", os.path.join(ROOT, "data"))


def rd(name, cols=None):
    d = pd.read_sas(os.path.join(DATA, name + ".xpt"))
    if cols:
        cols = ["SEQN"] + [c for c in cols if c in d.columns]
        d = d[cols]
    return d


def miss(s, codes=(7, 9, 77, 99, 777, 999, 7777, 9999)):
    return s.where(~s.isin(codes))


def build():
    d = rd("DEMO_L")
    d["age"] = d.RIDAGEYR
    d["female"] = (d.RIAGENDR == 2).astype(float)
    d["race"] = d.RIDRETH3.map({1: "Mexican American", 2: "Other Hispanic", 3: "NH White",
                                4: "NH Black", 6: "NH Asian", 7: "Other/Multiracial"})
    d["edu"] = miss(d.DMDEDUC2).map({1: "<High school", 2: "<High school", 3: "High school/GED",
                                     4: "Some college", 5: "College graduate+"})
    d["mar"] = miss(d.DMDMARTZ).map({1: "Married/partnered", 2: "Widowed/divorced/separated",
                                     3: "Never married"})
    d["pir"] = d.INDFMPIR
    d["pircat"] = pd.cut(d.pir, [-.01, 1.3, 3.5, 99], labels=["<1.3", "1.3-3.5", ">3.5"])
    d["usborn"] = (d.DMDBORN4 == 1).astype(float).where(d.DMDBORN4.isin([1, 2]))

    # ---- PHQ-9 -------------------------------------------------------------
    dp = rd("DPQ_L")
    items = [f"DPQ0{i}0" for i in range(1, 10)]
    for c in items:
        dp[c] = miss(dp[c])
    dp["nmiss"] = dp[items].isna().sum(axis=1)
    dp["phq9"] = dp[items].sum(axis=1).where(dp.nmiss == 0)
    dp["dep"] = (dp.phq9 >= 10).astype(float).where(dp.phq9.notna())
    dp["dep_modsev"] = (dp.phq9 >= 15).astype(float).where(dp.phq9.notna())
    dp["dep_mild"] = (dp.phq9 >= 5).astype(float).where(dp.phq9.notna())
    dp["si"] = (dp.DPQ090 >= 1).astype(float).where(dp.DPQ090.notna())
    dp["phq_impair"] = miss(dp.DPQ100).map({0: 0, 1: 1, 2: 1, 3: 1})
    # symptom dimensions (Krause/Elhai): somatic = sleep, energy, appetite, psychomotor
    dp["phq_som"] = dp[["DPQ030", "DPQ040", "DPQ050", "DPQ080"]].sum(axis=1).where(dp.nmiss == 0)
    dp["phq_cog"] = dp[["DPQ010", "DPQ020", "DPQ060", "DPQ070", "DPQ090"]].sum(axis=1).where(dp.nmiss == 0)
    d = d.merge(dp[["SEQN", "phq9", "dep", "dep_modsev", "dep_mild", "si", "phq_impair",
                    "phq_som", "phq_cog"] + items], on="SEQN", how="left")

    # ---- Washington Group functioning -------------------------------------
    fn = rd("FNQ_L")
    wg = {"FNQ410": "wg_see", "FNQ430": "wg_hear", "FNQ440": "wg_walk", "FNQ450": "wg_comm",
          "FNQ460": "wg_cog", "FNQ470": "wg_self", "FNQ480": "wg_lift", "FNQ490": "wg_hand"}
    for k, v in wg.items():
        s = miss(fn[k])
        fn[v] = s
        fn[v + "_any"] = (s >= 2).astype(float).where(s.notna())
        fn[v + "_lot"] = (s >= 3).astype(float).where(s.notna())
    fn["wg_ndom"] = fn[[v + "_any" for v in wg.values()]].sum(axis=1)
    fn["wg_dis"] = (fn.FNDADI == 1).astype(float).where(fn.FNDADI.notna())
    fn["wg_anxfreq"] = miss(fn.FNQ510)
    fn["wg_depfreq"] = miss(fn.FNQ530)
    d = d.merge(fn[["SEQN"] + [v for v in wg.values()] + [v + "_any" for v in wg.values()] +
                   [v + "_lot" for v in wg.values()] +
                   ["wg_ndom", "wg_dis", "wg_anxfreq", "wg_depfreq"]], on="SEQN", how="left")

    # ---- sex steroid hormones ---------------------------------------------
    ts = rd("TST_L")
    hn = {"LBXTST": "tst", "LBXEST": "e2", "LBXSHBG": "shbg", "LBXAMH": "amh", "LBXFSH": "fsh",
          "LBXLUH": "lh", "LBXPG4": "prog", "LBXDHE": "dheas", "LBXAND": "andro",
          "LBX17H": "ohp17", "LBXESO": "e1", "LBXES1": "e1s"}
    for k, v in hn.items():
        ts[v] = ts[k]
    # free testosterone (Vermeulen) needs albumin - merged later
    d = d.merge(ts[["SEQN", "WTPH2YR"] + list(hn.values())], on="SEQN", how="left")

    # ---- reproductive history ---------------------------------------------
    rh = rd("RHQ_L")
    rh["menarche"] = rh.RHQ010.where(rh.RHQ010 < 30)
    rh["period12"] = miss(rh.RHQ031).map({1: 1, 2: 0})
    rh["noperiod_why"] = miss(rh.RHD043, (77, 99))
    rh["age_lastperiod"] = rh.RHQ060.where(rh.RHQ060 < 70)
    rh["parity"] = miss(rh.RHD167, (77, 99))
    rh["pregnow"] = miss(rh.RHD143).map({1: 1, 2: 0})
    rh["hyster"] = miss(rh.RHD280).map({1: 1, 2: 0})
    rh["ooph"] = miss(rh.RHQ305).map({1: 1, 2: 0})
    rh["age_ooph"] = rh.RHQ332.where(rh.RHQ332 < 70)
    d = d.merge(rh[["SEQN", "menarche", "period12", "noperiod_why", "age_lastperiod",
                    "parity", "pregnow", "hyster", "ooph", "age_ooph"]], on="SEQN", how="left")

    # ---- balance / vestibular ---------------------------------------------
    bq = rd("BAQ_L")
    for c in ["BAQ321A", "BAQ321B", "BAQ321C", "BAQ321D"]:
        bq[c] = miss(bq[c])
    bq["dizzy_any"] = (bq[["BAQ321A", "BAQ321B", "BAQ321C", "BAQ321D"]] == 1).any(axis=1).astype(float)
    bq["vertigo"] = (bq.BAQ321A == 1).astype(float).where(bq.BAQ321A.notna())
    bq["lighthead"] = (bq.BAQ321B == 1).astype(float).where(bq.BAQ321B.notna())
    bq["floating"] = (bq.BAQ321C == 1).astype(float).where(bq.BAQ321C.notna())
    bq["imbalance"] = (bq.BAQ321D == 1).astype(float).where(bq.BAQ321D.notna())
    bq["dizzy_freq"] = miss(bq.BAQ341)
    bq["fall12"] = miss(bq.BAQ491).map({1: 1, 2: 0})
    d = d.merge(bq[["SEQN", "dizzy_any", "vertigo", "lighthead", "floating", "imbalance",
                    "dizzy_freq", "fall12"]], on="SEQN", how="left")

    bx = rd("BAX_L")
    # modified Romberg: condition 4 = eyes closed, compliant surface (vestibular)
    for cond in [1, 2, 3, 4]:
        t1 = miss(bx[f"BAXPF{cond}1"], (7, 9))
        t2 = miss(bx[f"BAXPF{cond}2"], (7, 9))
        # 1 = pass, 2 = fail ; fail overall if failed both trials (NHANES protocol: trial2 only if trial1 fails)
        fail = np.where(t1 == 1, 0.0, np.where(t2.notna(), (t2 == 2).astype(float), np.nan))
        bx[f"romberg_fail{cond}"] = fail
    bx["romberg_ok"] = (bx.BAXMSTAT == 1).astype(float)
    d = d.merge(bx[["SEQN", "romberg_ok"] + [f"romberg_fail{c}" for c in [1, 2, 3, 4]]],
                on="SEQN", how="left")

    # ---- folate forms ------------------------------------------------------
    ff = rd("FOLFMS_L")
    ff["fol_tot"] = ff.LBDFOTSI
    ff["fol_5mthf"] = ff.LBXSF1SI
    ff["fol_fa"] = ff.LBXSF2SI          # unmetabolized folic acid
    ff["fol_mefox"] = ff.LBXSF6SI
    d = d.merge(ff[["SEQN", "fol_tot", "fol_5mthf", "fol_fa", "fol_mefox"]], on="SEQN", how="left")
    d = d.merge(rd("FOLATE_L", ["LBDRFOSI"]).rename(columns={"LBDRFOSI": "rbc_folate"}),
                on="SEQN", how="left")

    # ---- other labs --------------------------------------------------------
    d = d.merge(rd("HSCRP_L", ["LBXHSCRP"]).rename(columns={"LBXHSCRP": "crp"}), on="SEQN", how="left")
    d = d.merge(rd("BCHE_L", ["LBXBA", "LBXCH"]).rename(columns={"LBXBA": "bche", "LBXCH": "ache"}),
                on="SEQN", how="left")
    d = d.merge(rd("GHB_L", ["LBXGH"]).rename(columns={"LBXGH": "hba1c"}), on="SEQN", how="left")
    d = d.merge(rd("VID_L", ["LBXVIDMS"]).rename(columns={"LBXVIDMS": "vitd"}), on="SEQN", how="left")
    d = d.merge(rd("BIOPRO_L", ["LBXSAL", "LBXSCR", "LBXSGTSI", "LBXSATSI", "LBXSASSI",
                                "LBXSTP", "LBXSGL"]), on="SEQN", how="left")
    d = d.merge(rd("CBC_L", ["LBXWBCSI", "LBDLYMNO", "LBDNENO", "LBXPLTSI", "LBXRDW", "LBXHGB",
                             "LBXMCVSI"]), on="SEQN", how="left")
    d = d.merge(rd("LUX_L", ["LUXSMED", "LUXCAPM"]).rename(
        columns={"LUXSMED": "lsm", "LUXCAPM": "cap"}), on="SEQN", how="left")

    # ---- covariates --------------------------------------------------------
    d = d.merge(rd("BMX_L", ["BMXBMI", "BMXWAIST"]), on="SEQN", how="left")
    d["bmi"] = d.BMXBMI
    sm = rd("SMQ_L")
    sm["smoke"] = np.where(sm.SMQ020 == 2, "Never",
                           np.where(sm.SMQ040.isin([1, 2]), "Current",
                                    np.where(sm.SMQ040 == 3, "Former", None)))
    d = d.merge(sm[["SEQN", "smoke"]], on="SEQN", how="left")
    d = d.merge(rd("COT_L", ["LBXCOT"]).rename(columns={"LBXCOT": "cotinine"}), on="SEQN", how="left")
    al = rd("ALQ_L")
    al["alq121"] = miss(al.ALQ121)
    al["drinker"] = np.where(al.ALQ111 == 2, "None",
                             np.where(al.alq121 == 0, "None", "Current"))
    d = d.merge(al[["SEQN", "drinker", "alq121"]], on="SEQN", how="left")
    pa = rd("PAQ_L")
    pa["sed_min"] = pa.PAD680.where(pa.PAD680 < 1300)
    for q, u, m in [("PAD790Q", "PAD790U", "PAD800"), ("PAD810Q", "PAD810U", "PAD820")]:
        pass
    d = d.merge(pa[["SEQN", "sed_min"]], on="SEQN", how="left")
    sl = rd("SLQ_L")
    sl["sleep_wk"] = sl.SLD012.where(sl.SLD012 < 25)
    sl["sleep_we"] = sl.SLD013.where(sl.SLD013 < 25)
    d = d.merge(sl[["SEQN", "sleep_wk", "sleep_we"]], on="SEQN", how="left")
    fs = rd("FSQ_L", ["FSDAD"])
    fs["foodinsec"] = (fs.FSDAD >= 3).astype(float).where(fs.FSDAD.notna())
    d = d.merge(fs[["SEQN", "foodinsec"]], on="SEQN", how="left")
    hi = rd("HIQ_L")
    hi["insured"] = miss(hi.HIQ011).map({1: 1, 2: 0})
    d = d.merge(hi[["SEQN", "insured"]], on="SEQN", how="left")
    hu = rd("HUQ_L")
    hu["srh"] = miss(hu.HUQ010)
    hu["mh_visit"] = miss(hu.HUQ090).map({1: 1, 2: 0})
    hu["hc_visits"] = miss(hu.HUQ042, (77, 99))
    d = d.merge(hu[["SEQN", "srh", "mh_visit", "hc_visits"]], on="SEQN", how="left")
    mc = rd("MCQ_L")
    for k, v in {"MCQ160B": "chf", "MCQ160C": "chd", "MCQ160E": "mi", "MCQ160F": "stroke",
                 "MCQ160M": "thyroid", "MCQ160P": "copd", "MCQ160L": "liver", "MCQ160A": "arthritis",
                 "MCQ220": "cancer", "MCQ010": "asthma"}.items():
        mc[v] = miss(mc[k]).map({1: 1, 2: 0})
    mc["ncond"] = mc[["chf", "chd", "mi", "stroke", "copd", "liver", "arthritis", "cancer"]].sum(axis=1)
    d = d.merge(mc[["SEQN", "chf", "chd", "mi", "stroke", "thyroid", "copd", "liver",
                    "arthritis", "cancer", "asthma", "ncond"]], on="SEQN", how="left")
    di = rd("DIQ_L")
    di["dm"] = miss(di.DIQ010).map({1: 1, 2: 0, 3: 0})
    d = d.merge(di[["SEQN", "dm"]], on="SEQN", how="left")
    bp = rd("BPQ_L")
    bp["htn"] = miss(bp.BPQ020).map({1: 1, 2: 0})
    d = d.merge(bp[["SEQN", "htn"]], on="SEQN", how="left")
    ki = rd("KIQ_U_L")
    ki["ckd_sr"] = miss(ki.KIQ022).map({1: 1, 2: 0})
    d = d.merge(ki[["SEQN", "ckd_sr"]], on="SEQN", how="left")
    au = rd("AUQ_L")
    au["hear_diff"] = miss(au.AUQ054)
    au["tinnitus"] = miss(au.AUQ144)
    d = d.merge(au[["SEQN", "hear_diff", "tinnitus"]], on="SEQN", how="left")
    oh = rd("OHQ_L")
    oh["oral_health"] = miss(oh.OHQ845)
    d = d.merge(oh[["SEQN", "oral_health"]], on="SEQN", how="left")
    oc = rd("OCQ_L")
    oc["work"] = miss(oc.OCD150).map({1: "Working", 2: "Working", 3: "Not working", 4: "Not working"})
    oc["hours"] = oc.OCQ180.where(oc.OCQ180 <= 120)
    d = d.merge(oc[["SEQN", "work", "hours"]], on="SEQN", how="left")
    ds = rd("DSQTOT_L", ["DSQTFA", "DSQTVB12", "DSQTVD", "DSD010"])
    ds["supp_any"] = miss(ds.DSD010).map({1: 1, 2: 0})
    d = d.merge(ds[["SEQN", "DSQTFA", "DSQTVB12", "DSQTVD", "supp_any"]], on="SEQN", how="left")
    dr = rd("DR1TOT_L", ["WTDRD1", "DR1TKCAL", "DR1TFOLA", "DR1TFA", "DR1TVB12", "DR1TALCO"])
    d = d.merge(dr, on="SEQN", how="left")
    rx = rd("RXQ_RX_L")
    rx["rx_any"] = miss(rx.RXQ033).map({1: 1, 2: 0})
    rx["rx_n"] = rx.RXQ050.where(rx.RXQ050 < 50)
    d = d.merge(rx[["SEQN", "rx_any", "rx_n"]], on="SEQN", how="left")

    d["wmec"] = d.WTMEC2YR.fillna(0.0)
    d["wph"] = d.WTPH2YR.fillna(0.0)
    return d


if __name__ == "__main__":
    d = build()
    d.to_pickle(os.path.join(DATA, "wide.pkl"))
    print(d.shape)
    ad = d[(d.age >= 18) & (d.phq9.notna())]
    print("adults with PHQ-9:", len(ad))
    for c in ["amh", "fsh", "lh", "prog", "e2", "tst", "shbg", "dheas", "romberg_fail4",
              "dizzy_any", "fol_5mthf", "fol_fa", "wg_ndom", "lsm", "bche"]:
        print(f"  {c:16s} n={ad[c].notna().sum():5d}")
