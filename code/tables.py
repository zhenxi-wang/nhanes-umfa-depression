"""Emit Table 1-4 and all supplementary tables as markdown."""
import json
import os
import numpy as np
import pandas as pd

A1 = json.load(open("out_a1.json"))
A2 = json.load(open("out_a2.json"))
A3 = json.load(open("out_a3.json"))
A4 = json.load(open("out_a4.json"))
A5 = json.load(open("out_a5.json"))
A6 = json.load(open("out_a6.json"))
A7 = json.load(open("out_a7.json"))
FL = json.load(open("flow.json"))
out, sup = [], []

LAB = {"5-MTHF (per SD log)": "5-Methyltetrahydrofolate",
       "Total serum folate (per SD log)": "Total serum folate",
       "RBC folate (per SD log)": "Red-cell folate",
       "MeFox (per SD log)": "MeFox (oxidation product)",
       "Unmetabolised folic acid (per SD log)": "Unmetabolised folic acid",
       "UMFA as % of total folate (per SD log)": "UMFA fraction (% of quantified folate)"}
QLAB = {"5-MTHF": "5-Methyltetrahydrofolate", "UMFA": "Unmetabolised folic acid",
        "Total folate": "Total serum folate",
        "%UMFA": "UMFA fraction (% of quantified folate)"}


def md(df):
    hdr = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "| " + " | ".join("---" for _ in df.columns) + " |"
    rows = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
            for r in df.itertuples(index=False)]
    return "\n".join([hdr, sep] + rows)


def clean(recs, cols=None, rename=None):
    t = pd.DataFrame(recs)
    t = t[[c for c in t.columns if not c.startswith("_")]]
    if cols:
        t = t[cols]
    if rename:
        t = t.rename(columns=rename)
    return t


# ============================================================== Table 1
t1 = pd.DataFrame(A4["table1"]).rename(columns={"Q1 (lowest)": "Q1 (lowest UMFA)",
                                                "Q4 (highest)": "Q4 (highest UMFA)"})
out.append("## Table 1. Characteristics of 24,326 US adults by cycle-specific quartile of serum "
           "unmetabolised folic acid, NHANES 2011-2023\n")
out.append(md(t1))
out.append("\nValues are survey-weighted percentages or means unless stated; unweighted n is given "
           "for sample size only. Medians and interquartile ranges are unweighted sample quantiles, "
           "as are the quartile cut-points; the survey-weighted counterparts are given in "
           "Supplementary Table S14B. No significance tests are presented for baseline characteristics; "
           "standardised mean differences between the extreme quartiles are given in Supplementary "
           "Table S12. GM, geometric mean. Supplemental folic acid is set to zero in non-users. "
           "eGFR, estimated glomerular filtration rate (2021 CKD-EPI creatinine equation). PHQ-9, "
           "nine-item Patient Health Questionnaire; suicidal ideation, any endorsement of item 9.\n")

# ============================================================== Table 2
persd = pd.DataFrame(A1["persd"])
rows = []
for ex in LAB:
    _n = A6["biomarker_n"].get(LAB[ex], A6["biomarker_n"].get(LAB[ex].split(" (")[0]))
    r = {"Biomarker": LAB[ex], "n": f"{_n:,}"}
    for mo in [1, 2, 3]:
        x = persd[(persd.exposure == ex) & (persd.model == mo)].iloc[0]
        r[f"Model {mo}"] = x.OR
        r[f"p (M{mo})"] = x.p
    rows.append(r)
out.append("\n\n## Table 2. Association of each serum folate species with clinically significant "
           "depressive symptoms (PHQ-9 >= 10), per standard deviation of the log-transformed "
           "concentration\n")
out.append("**A. Single-biomarker models**\n")
out.append(md(pd.DataFrame(rows)))
mut = pd.DataFrame(A1["mutual"])
rows = []
for term, nm in [("5-MTHF", "5-Methyltetrahydrofolate"), ("UMFA", "Unmetabolised folic acid"),
                 ("MeFox", "MeFox")]:
    r = {"Biomarker": nm}
    for mo in [1, 2, 3]:
        x = mut[(mut.term == term) & (mut.model == mo)].iloc[0]
        r[f"Model {mo}"] = x.OR
        r[f"p (M{mo})"] = x.p
    rows.append(r)
out.append("\n**B. Mutually adjusted model containing all three species simultaneously**\n")
out.append(md(pd.DataFrame(rows)))
c = A1["contrast_m3"]
out.append(f"\nModel 1: age, age-squared, sex, race and Hispanic origin, survey cycle. Model 2: "
           f"+ education, income-to-poverty ratio, marital status, insurance, food security, "
           f"nativity. Model 3 (extended, primary): + smoking status and cotinine, alcohol, body "
           f"mass index, sedentary time, physical activity, sleep duration, chronic conditions, "
           f"diabetes, eGFR, energy intake, dietary folate, dietary folic acid, supplement use, "
           f"supplemental folic acid, fasting hours and examination session. In panel B, model 3, "
           f"the contrast between the 5-MTHF and UMFA coefficients gave a ratio of odds ratios of "
           f"{c['ratio']:.2f} (95% CI {c['lo']:.2f}-{c['hi']:.2f}; p = {c['p']:.3f}). Alternative "
           f"model specifications, including a directed-acyclic-graph-guided minimal adjustment "
           f"set, are given in Supplementary Table S4A.\n")

# ============================================================== Table 3
q = pd.DataFrame(A1["quartiles"])
ordt = A4["ordinal_trend"]
rows = []
for ex in ["5-MTHF", "UMFA", "Total folate", "%UMFA"]:
    r = q[q.exposure == ex].iloc[0]
    key = {"5-MTHF": "5-MTHF", "UMFA": "UMFA", "Total folate": "Total folate",
           "%UMFA": "%UMFA"}[ex]
    rows.append({"Biomarker": QLAB[ex], "Q1": "1.00 (ref)", "Q2": r.Q2, "Q3": r.Q3, "Q4": r.Q4,
                 "p (joint)": r.p_joint, "p (trend, quartile median)": r.p_trend,
                 "p (trend, quartile rank)": f"{ordt[key]['p']:.3f}"})
out.append("\n\n## Table 3. Depressive symptoms by cycle-specific quartile of folate species, "
           "extended model (Model 3)\n")
out.append(md(pd.DataFrame(rows)))
out.append("\nOdds ratios (95% confidence intervals) for PHQ-9 >= 10. Each biomarker is modelled "
           "separately: the quartile indicators for one species are entered into Model 3 without "
           "the other species, so these are not the mutually adjusted estimates of Table 2B. "
           "Quartiles were defined within survey cycle. p (joint) is a 3-degree-of-freedom Wald "
           "test. Both trend tests are "
           "reported because they can differ: one fits the median of each quartile, the other the "
           "quartile rank. Neither the joint test nor either trend test for unmetabolised folic acid "
           "reaches conventional significance, and the quartile pattern is not strictly monotonic; "
           "the continuous and spline analyses (Table 2, Figure 1) are the more informative "
           "parameterisations.\n")
thr = clean(A3["thresholds"])
thr.columns = ["Threshold", "n exposed", "Prevalence if exposed (%)",
               "Prevalence if unexposed (%)", "Adjusted OR (95% CI)", "p"]
out.append("\n**Post-hoc threshold analyses (Model 3, additionally adjusted for 5-MTHF and MeFox)**\n")
out.append(md(thr))
out.append("\nThese cut-points were chosen after inspecting the distribution, for interpretability, "
           "and were not specified in advance. Lower cut-points show no crude gradient, so these "
           "results should not be read as a clean threshold effect.\n")

# ============================================================== Table 4
o = pd.DataFrame(A1["outcomes"])
o2 = o[["outcome", "n", "UMFA", "p_UMFA", "5-MTHF", "p_5-MTHF", "MeFox", "p_MeFox"]]
o2.columns = ["Outcome", "n", "UMFA", "p", "5-MTHF", "p ", "MeFox", "p  "]
out.append("\n\n## Table 4. Secondary outcomes, mutually adjusted extended model\n")
out.append(md(o2))
out.append("\nBinary outcomes are odds ratios; the last three rows are linear-regression "
           "coefficients per standard deviation of the log biomarker. Somatic subscale = PHQ-9 "
           "items 3, 4, 5, 8; cognitive-affective subscale = items 1, 2, 6, 7, 9. All outcomes "
           "other than PHQ-9 >= 10 are exploratory; see Supplementary Table S11 for "
           "false-discovery-rate-adjusted q values.\n")

# ============================================================== Supplement
sup.append("# Supplementary Material\n")
sup.append("## Supplementary Table S1. Participant flow and exclusions\n")
sup.append(md(pd.DataFrame([{"Step": a.replace("\n", " "), "n remaining": f"{b:,}",
                             "Excluded at this step": ("-" if c is None
                                                       else f"{c[0]} (n = {c[1]:,})")}
                            for a, b, c in FL["steps"]])))
sup.append("\n**Analytic sample by cycle:** " +
           "; ".join(f"{k} n = {v:,}" for k, v in FL["by_cycle"].items()) + ".\n")

sup.append("\n\n## Supplementary Table S2. Weighted distribution of folate species by survey cycle\n")
sup.append(md(pd.DataFrame(A6["by_cycle_desc"])))
sup.append("\nGeometric means (95% confidence intervals) in nmol/L for the folate species, and "
           "survey-weighted outcome prevalences with 95% confidence intervals. RBC, red blood cell. "
           "See Supplementary Figure S3 for the corresponding trend plot.\n")

sup.append("\n\n## Supplementary Table S3. Subgroup analyses: UMFA per SD of the log "
           "concentration, mutually adjusted three-species Model 3\n")
sg = clean(A2["subgroups"], ["subgroup", "n", "OR", "p"],
           {"subgroup": "Subgroup stratum", "OR": "Adjusted OR (95% CI)"})
sup.append(md(sg))
it = A2["interaction_tests"]
sup.append("\nParticipants with a missing value for the variable defining a stratum are omitted "
           "from that stratified analysis, so strata do not sum to 24,326. These are strata, not "
           "independent comparisons. Formal multiplicative interaction "
           "tests: " + "; ".join(
               f"UMFA x {k}: OR {v['OR']:.2f} ({v['lo']:.2f}-{v['hi']:.2f}), p = {v['p']:.2f}"
               for k, v in it.items()) +
           f"; UMFA x folic-acid-containing supplement use: OR "
           f"{A5['fa_supp_interaction']['OR']:.2f} "
           f"({A5['fa_supp_interaction']['lo']:.2f}-{A5['fa_supp_interaction']['hi']:.2f}), "
           f"p = {A5['fa_supp_interaction']['p']:.2f}. There is no evidence of effect modification "
           f"in any stratum examined.\n")

sup.append("\n\n## Supplementary Table S4. Sensitivity analyses\n")
sup.append("**A. Model specification**\n")
sup.append(md(clean(A5["specifications"], ["analysis", "n", "OR", "p"],
                    {"analysis": "Specification", "OR": "UMFA OR (95% CI)"})))
sup.append("\nModel 5 retains only covariates judged upstream of both circulating folate species "
           "and depressive symptoms (age, sex, race and Hispanic origin, cycle, education, "
           "income-to-poverty ratio, food security, nativity, smoking and cotinine, eGFR, energy "
           "intake, dietary folate and folic acid, supplement use, supplemental folic acid, fasting "
           "hours, examination session) and drops those plausibly downstream of depressive symptoms "
           "(body mass index, sleep, physical activity, alcohol, chronic conditions, diabetes, "
           "marital status, insurance).\n")
sup.append("\n**B. Survey cycle and measurement mode**\n")
sup.append(md(clean(A5["cycle_mode"], ["analysis", "n", "OR", "p"],
                    {"analysis": "Analysis", "OR": "UMFA OR (95% CI)"})))
sup.append("\nThe dietary-supplement and prescription-medication instruments changed in August "
           "2021-August 2023 (supplement interviewing moved to a post-examination telephone "
           "contact, the 24-hour supplement recall was discontinued and drug-level prescription "
           "medication information was no longer collected, although participants were still asked "
           "whether they took prescription medicines and how many), so the analysis excluding that "
           "cycle is the relevant check on measurement-mode effects.\n")
sup.append("\n**C. Folic-acid-specific supplement use, fasting and examination timing**\n")
sup.append(md(clean(A5["supp_fast"], ["analysis", "n", "OR", "p"],
                    {"analysis": "Stratum", "OR": "UMFA OR (95% CI)"})))
sup.append("\n**D. Medication, exclusions and alternative adjustments**\n")
sn = pd.concat([pd.DataFrame(A2["sensitivity"])[["analysis", "n", "_or", "_lo", "_hi", "_p"]],
                pd.DataFrame(A4["robust_extra"])[["analysis", "n", "_or", "_lo", "_hi", "_p"]],
                pd.DataFrame(A3["medication"])[["analysis", "n", "_or", "_lo", "_hi", "_p"]]],
               ignore_index=True).drop_duplicates("analysis")
sn["UMFA OR (95% CI)"] = [f"{a:.2f} ({b:.2f}-{c:.2f})" for a, b, c in zip(sn._or, sn._lo, sn._hi)]
sn["p"] = ["<0.001" if p < 0.001 else f"{p:.3f}" for p in sn._p]
sn.loc[sn.analysis == "Main (model 3)", "analysis"] = \
    "Reference specification (mutually adjusted three-species Model 3)"
sn = pd.concat([sn[["analysis", "n", "UMFA OR (95% CI)", "p"]],
                pd.DataFrame([{"analysis": "Antidepressant use as the outcome (rather than "
                                           "depressive symptoms)",
                               "n": A5["ad_as_outcome"]["n"],
                               "UMFA OR (95% CI)": A5["ad_as_outcome"]["OR"],
                               "p": A5["ad_as_outcome"]["p"]}])], ignore_index=True)
sup.append(md(sn[["analysis", "n", "UMFA OR (95% CI)", "p"]].rename(
    columns={"analysis": "Analysis"})))
ev = A2["evalue"]
mi_ = A5["mice"]
sup.append(f"\nEvery row uses the mutually adjusted three-species Model 3 specification. The "
           f"UMFA x antidepressant-use interaction was not significant "
           f"(p = {A3['ad_interaction']['p']:.2f}), so the user and non-user strata should not be "
           f"read as differing. "
           f"Multiple imputation by chained equations (m = {mi_['m']}, Rubin rules): OR "
           f"{mi_['OR']:.2f} ({mi_['lo']:.2f}-{mi_['hi']:.2f}), p = {mi_['p']:.3f}. Unweighted "
           f"estimate: OR {A4['unweighted']['OR']:.2f} ({A4['unweighted']['lo']:.2f}-"
           f"{A4['unweighted']['hi']:.2f}). E-value for the primary estimate: {ev['point']:.2f} "
           f"(confidence limit {ev['ci']:.2f}); because the outcome is not rare, the odds ratio was "
           f"first converted to an approximate risk ratio by the conservative square-root "
           f"transformation before applying the VanderWeele-Ding formula. Applying the formula "
           f"directly to the odds ratio would give a larger E-value.\n")

sup.append("\n\n## Supplementary Table S5. Residual UMFA after accounting for measured intake "
           "and specimen-related predictors\n")
sup.append(md(clean(A3["residual_umfa"]["models"], ["model", "n", "OR", "p"],
                    {"model": "Model", "OR": "Adjusted OR (95% CI)"})))
sup.append(f"\nResidual UMFA is the standardised residual of log-UMFA regressed on dietary folic "
           f"acid, supplemental folic acid, dietary folate, supplement use, energy intake, fasting "
           f"hours, examination session, age, sex, body mass index, eGFR and survey cycle "
           f"(model R-squared = {A3['residual_umfa']['r2']:.3f}). Because these predictors explain less than a fifth "
           f"of the variance, the residual also contains dietary and supplement recall error, "
           f"unmeasured dose and formulation, time since last intake, absorption, hepatic handling, "
           f"renal clearance and assay variability; it should not be interpreted as a measure of "
           f"dihydrofolate reductase activity. The outcome models additionally contained "
           f"5-MTHF and MeFox, as elsewhere.\n")

sup.append("\n\n## Supplementary Table S6. Exploratory folate-conversion phenotype (median splits "
           "within cycle)\n")
pv = pd.DataFrame(A2["phenotype"]["prev"])
pv["Unadjusted prevalence (95% CI), %"] = [f"{a:.2f} ({b:.2f}-{c:.2f})"
                                           for a, b, c in zip(pv.prev, pv.lo, pv.hi)]
mo1 = [r for r in A2["phenotype"]["models"] if r["model"] == 1][0]
mo3 = [r for r in A2["phenotype"]["models"] if r["model"] == 3][0]
pv["Model 1 OR"] = [mo1.get(g, "1.00 (ref)") for g in pv.group]
pv["Model 3 OR"] = [mo3.get(g, "1.00 (ref)") for g in pv.group]
sup.append(md(pv[["group", "n", "Unadjusted prevalence (95% CI), %", "Model 1 OR", "Model 3 OR"]]
              .rename(columns={"group": "Group", "n": "Unweighted n"})))
ic = A2["interaction_cont"]
sup.append(f"\nA phenotype combining relatively low 5-MTHF with relatively high UMFA - the pattern "
           f"expected if impaired folic-acid conversion drove the association - showed a strong "
           f"unadjusted association that disappeared after multivariable adjustment. The continuous "
           f"UMFA-by-5-MTHF interaction was also null (OR {ic['OR']:.2f}, {ic['lo']:.2f}-"
           f"{ic['hi']:.2f}, p = {ic['p']:.2f}). These data therefore do not establish impaired "
           f"conversion as the mechanism underlying the UMFA association.\n")

b12 = A3["b12_interaction"]
sup.append(f"\n\n## Supplementary Table S7. Vitamin B12 interaction (2011-2014, n = {b12['n']:,})\n")
sup.append(md(pd.DataFrame([{"Stratum": r["stratum"],
                             "UMFA OR (95% CI)": f"{r['OR']:.2f} ({r['lo']:.2f}-{r['hi']:.2f})",
                             "p": f"{r['p']:.3f}"} for r in A3["b12"]])))
sup.append(f"\nInteraction UMFA x low vitamin B12 (<221 pmol/L): OR {b12['OR']:.2f} "
           f"({b12['lo']:.2f}-{b12['hi']:.2f}), p = {b12['p']:.2f}.\n")

# S8 exploratory screen
sup.append("\n\n## Supplementary Table S8. Exploratory screen that led to the choice of exposure\n")
if os.path.exists("screen.csv"):
    sc = pd.read_csv("screen.csv")
    NAME = {"log_fol_5mthf": "5-Methyltetrahydrofolate", "log_fol_fa": "Unmetabolised folic acid",
            "log_fol_mefox": "MeFox", "log_rbc_folate": "Red-cell folate",
            "logamh": "Anti-Mullerian hormone", "log_fsh": "Follicle-stimulating hormone",
            "log_lh": "Luteinising hormone", "log_e2": "Oestradiol", "log_prog": "Progesterone",
            "log_tst": "Total testosterone", "log_shbg": "Sex hormone-binding globulin",
            "log_dheas": "DHEAS", "log_andro": "Androstenedione",
            "log_ohp17": "17-hydroxyprogesterone", "log_e1": "Oestrone",
            "log_e1s": "Oestrone sulfate",
            "romberg_fail4": "Modified Romberg, condition 4",
            "romberg_fail3": "Modified Romberg, condition 3",
            "log_lsm": "Liver stiffness", "log_cap": "Controlled attenuation parameter",
            "log_bche": "Butyrylcholinesterase", "log_vitd": "25-hydroxyvitamin D"}
    DOM = {"5-Methyltetrahydrofolate": "Serum folate vitamers",
           "Unmetabolised folic acid": "Serum folate vitamers", "MeFox": "Serum folate vitamers",
           "Red-cell folate": "Serum folate vitamers",
           "Anti-Mullerian hormone": "Sex steroid hormones",
           "Follicle-stimulating hormone": "Sex steroid hormones",
           "Luteinising hormone": "Sex steroid hormones", "Oestradiol": "Sex steroid hormones",
           "Progesterone": "Sex steroid hormones", "Total testosterone": "Sex steroid hormones",
           "Sex hormone-binding globulin": "Sex steroid hormones", "DHEAS": "Sex steroid hormones",
           "Androstenedione": "Sex steroid hormones",
           "17-hydroxyprogesterone": "Sex steroid hormones", "Oestrone": "Sex steroid hormones",
           "Oestrone sulfate": "Sex steroid hormones",
           "Modified Romberg, condition 4": "Objective balance",
           "Modified Romberg, condition 3": "Objective balance",
           "Liver stiffness": "Liver elastography",
           "Controlled attenuation parameter": "Liver elastography",
           "Butyrylcholinesterase": "Butyrylcholinesterase", "25-hydroxyvitamin D": "Vitamin D"}
    sc["Biomarker"] = sc.exposure.map(NAME)
    sc = sc[sc.Biomarker.notna()].copy()
    sc["Domain"] = sc.Biomarker.map(DOM)
    sc["Group"] = sc["group"].str.replace(" (binary, per SD)", "", regex=False)
    sc = sc.sort_values(["Domain", "p"])
    sc["OR (95% CI)"] = [f"{o:.2f} ({l:.2f}-{h:.2f})" for o, l, h in zip(sc.OR, sc.lo, sc.hi)]
    sc["p"] = [f"{v:.4f}" for v in sc.p]
    sup.append(md(sc[["Domain", "Biomarker", "Group", "n", "OR (95% CI)", "p"]].rename(
        columns={"Group": "Sample screened"})))
if os.path.exists("screen_folate_final.json"):
    sfin = pd.DataFrame(json.load(open("screen_folate_final.json")))
    sup.append("\n**Folate vitamers re-screened on the final discovery eligibility "
               "(n = 5,028)**\n")
    sup.append(md(sfin.rename(columns={"OR": "OR (95% CI)"})))
    sup.append("\nThe screen as executed used all August 2021-August 2023 participants with a "
               "folate result and a complete PHQ-9 (n = 5,060), before the pregnancy and "
               "analytic-weight exclusions that define the discovery sample. Repeating it on the "
               "final discovery eligibility changes nothing material.\n")
sup.append("\nCandidate biomarker domains available in the August 2021-August 2023 cycle were "
           "screened against PHQ-9-defined depressive symptoms using an age-, sex-, race-, "
           "education- and income-adjusted logistic model before the folate vitamers were selected "
           "for full analysis. Folate was not the domain with the smallest screening p value; it "
           "was chosen for biological specificity, novelty and availability across earlier cycles. "
           "Because the exposure was selected after inspecting the outcome, the August "
           "2021-August 2023 cycle is treated as a discovery sample and the four earlier cycles as "
           "an out-of-screen replication sample (Supplementary Table S13).\n")

sup.append("\n\n## Supplementary Table S9. Covariate missingness in the analytic sample "
           "(n = 24,326)\n")
sup.append(md(pd.DataFrame(A5["missingness"])))
sup.append("\nContinuous and binary covariates were mean-imputed with a missingness indicator "
           "(algebraically equivalent to an explicit missing category for a binary variable); "
           "multi-level categorical covariates were assigned an explicit missing category rather "
           "than collapsed into the reference level. Two variables are largely structurally absent "
           "rather than missing at random: high-sensitivity C-reactive protein was not measured in "
           "2011-2012 or 2013-2014 (9,649 participants), leaving 153 additional specimen-level "
           "missing values in the cycles where it was measured; drug-level medication data were not "
           "released for August 2021-August 2023 (5,028 participants), with no additional "
           "missingness. Both are used only in cycle-restricted analyses. Complete-case and "
           "multiple-imputation analyses are in Supplementary Table S4D.\n")

sup.append("\n\n## Supplementary Table S10. Validation of the design-based estimator\n")
sup.append(md(pd.DataFrame(A5["validation"]).rename(
    columns={"estimate": "Estimate", "ours": "This analysis", "published": "Published figure"})))
sup.append("\nThe two figures are not expected to coincide exactly because the published estimate "
           "covers ages 12 years and older using the full examined sample, whereas this analysis is "
           "restricted to adults with folate vitamers measured and uses the corresponding subsample "
           "weight.\n")

sup.append("\n\n## Supplementary Table S11. Benjamini-Hochberg q values for secondary outcomes, "
           "all reported quartile inference and threshold analyses\n")
sup.append(md(pd.DataFrame(A5["fdr"])))
sup.append("\nThe family comprises the secondary outcomes, every quartile inference reported in "
           "Table 3 (the Q4-versus-Q1 contrast, the 3-degree-of-freedom joint test and both trend "
           "tests) and the threshold analyses - hypotheses about the same exposure examined in "
           "different ways. The nominally significant trend tests are included rather than set "
           "aside, since excluding them would make the correction self-serving. Subgroup "
           "strata are deliberately excluded because subgroup inference belongs to the interaction "
           "tests in Supplementary Table S3, and robustness analyses are excluded because they "
           "re-examine a single hypothesis rather than test new ones. **No exploratory analysis "
           "remains below a false-discovery rate of 0.05.** The post-hoc 5 nmol/L cut-point in "
           "addition cannot be corrected for the cut-points that were inspected but not reported.\n")

sup.append("\n\n## Supplementary Table S12. Standardised mean differences between extreme UMFA "
           "quartiles\n")
sup.append(md(pd.DataFrame(A4["smd"])))
sup.append("\nStandardised mean differences are reported in preference to significance testing of "
           "baseline characteristics. Q1 and Q4 values are survey-weighted arithmetic means or "
           "proportions, as required to compute standardised differences; for right-skewed "
           "biomarkers these therefore differ from the geometric means reported in Table 1.\n")

sup.append("\n\n## Supplementary Table S13. Discovery, out-of-screen replication and pooled "
           "estimate\n")
sup.append(md(clean(A6["discovery_replication"], ["sample", "n", "OR", "p"],
                    {"sample": "Sample", "OR": "UMFA OR (95% CI)"})))
sup.append("\nThe exposure screen (Supplementary Table S8) used only the August 2021-August 2023 "
           "cycle, so the four earlier cycles provide an estimate independent of the data in which "
           "the exposure was selected.\n")

sup.append("\n\n## Supplementary Table S14. Survey-weighted standardisation and cut-points\n")
sup.append("**A. Standardisation of the exposure**\n")
sup.append(md(clean(A6["weighted_sd"]["models"], ["analysis", "n", "OR", "p"],
                    {"analysis": "Standardisation", "OR": "UMFA OR (95% CI)"})))
sup.append("\nWithin-cycle standard deviation of log-UMFA, unweighted vs survey-weighted: " +
           "; ".join(f"{r['cycle']} {r['sd_unweighted']:.3f} vs {r['sd_weighted']:.3f}"
                     for r in A6["weighted_sd"]["sds"]) + ".\n")
sup.append("\n**B. Quartiles and upper-tail cut-points on the survey-weighted distribution**\n")
sup.append(md(pd.DataFrame([A6["weighted_quantiles"]["quartiles"]])))
sup.append("\n" + md(clean(A6["weighted_quantiles"]["thresholds"],
                           ["threshold", "n_exposed", "OR", "p"],
                           {"threshold": "Cut-point", "n_exposed": "n exposed",
                            "OR": "Adjusted OR (95% CI)"})))
sup.append("\nQuartiles and cut-points in the main text are unweighted sample quantiles computed "
           "within cycle. Defining them on the survey-weighted distribution leaves the top-quartile "
           "estimate essentially unchanged but widens its confidence interval; the continuous and "
           "spline analyses do not depend on this choice.\n")

sup.append("\n\n## Supplementary Table S15. Non-linear adjustment for key confounders\n")
sup.append(md(pd.DataFrame([{"Specification": "Model 3 with restricted cubic splines for age, "
                                              "eGFR, dietary folic acid, supplemental folic acid "
                                              "and fasting time",
                             "n": f"{A6['spline_confounders']['n']:,}",
                             "UMFA OR (95% CI)": A6["spline_confounders"]["OR"],
                             "p": A6["spline_confounders"]["p"]}])))
sup.append("\nAge, renal function, folic acid intake and fasting time show the largest imbalance "
           "across UMFA quartiles (Supplementary Table S12), so their functional form was relaxed "
           "as a robustness check.\n")

sup.append("\n\n## Supplementary Table S16. Participants excluded for an incomplete PHQ-9\n")
sup.append("**A. Survey-weighted comparison with standardised mean differences**\n")
sup.append(md(clean(A6["selection_weighted"])))
sup.append("\n**B. Unweighted descriptive comparison**\n")
sup.append(md(pd.DataFrame(A6["selection"]["table"])))
sup.append("\nProportion excluded by cycle: " +
           "; ".join(f"{k} {v}%" for k, v in A6["selection"]["by_cycle"].items()) +
           ". Panel A is survey-weighted and reports standardised mean differences in preference "
           "to significance tests; panel B gives the corresponding unweighted descriptive values. "
           "Excluded participants had somewhat lower income and were less often non-Hispanic "
           "White, but the exposure distribution was very similar (standardised mean difference in "
           "log-UMFA -0.04).\n")

sup.append("\n\n## Supplementary Table S17. Participants excluded because folate vitamers were "
           "not quantified\n")
sup.append(md(pd.DataFrame(A6["selection_folate"]["table"])))
sup.append("\nUnweighted descriptive values. Proportion without a folate assay by cycle: " +
           "; ".join(f"{k} {v}%" for k, v in A6["selection_folate"]["by_cycle"].items()) +
           ". Those without an assay had lower income, were less often non-Hispanic White and had "
           "slightly more depressive symptoms (11.0% vs 9.9% with a PHQ-9 of 10 or more), so this "
           "exclusion is not completely at random; the direction implies that, if anything, "
           "including them would raise rather than lower the outcome prevalence in the analytic "
           "sample.\n")

sup.append("\n\n## Supplementary Table S18. Correlations among folate biomarkers and "
           "collinearity of the primary model\n")
sup.append(md(pd.DataFrame(A6["corr"])))
_cl = A6["collinearity"]
sup.append("\nPearson correlations on the log scale in the analytic sample. In the mutually "
           "adjusted three-species Model 3 the variance inflation factors were " +
           "; ".join(f"{k} {v}" for k, v in _cl["vif"].items()) +
           f". The {_cl['ncol']} predictor columns are linearly independent; including the intercept the "
           f"design matrix has rank 61 of 61 (Supplementary Table S10). The condition number of the "
           f"standardised predictor matrix was "
           f"{_cl['condition_number']}. Collinearity is therefore not a material concern despite "
           "the correlation between 5-MTHF and UMFA. Redundant columns are removed at "
           "construction: energy intake, dietary folate and dietary folic acid come from the same "
           "24-hour recall and so share a single missingness indicator, and the evening "
           "examination session is empty in the analytic sample.\n")

sup.append("\n\n## Supplementary Table S19. Severity gradient across PHQ-9 cut-points\n")
sup.append(md(clean(A6["ordinal"]["rows"], ["cutpoint", "n", "OR", "p"],
                    {"cutpoint": "Cut-point", "OR": "UMFA OR (95% CI)"})))
_o = A6["ordinal"]
sup.append("\nEach row is a separate survey-weighted logistic model for the cumulative outcome at "
           "that cut-point, using the mutually adjusted three-species Model 3. No formal test of "
           "equality across cut-points is reported: the three outcomes are nested, so the "
           "coefficients are strongly dependent and a test treating them as independent would be "
           "anti-conservative. The estimates are presented descriptively. The pattern is "
           "compatible with an association concentrated at clinically significant symptom levels, "
           "but these exploratory estimates should not be read as evidence of statistically "
           "different effects across thresholds.\n")

sup.append("\n\n## Supplementary Table S20. Covariate definitions and NHANES source variables\n")
DICT = [
 ("Age, sex, race and Hispanic origin", "RIDAGEYR, RIAGENDR, RIDRETH3", "same",
  "Age modelled as linear and quadratic terms; race and Hispanic origin as indicators with "
  "non-Hispanic White as reference."),
 ("Survey cycle", "SDDSRVYR / file suffix", "same", "Indicator per cycle."),
 ("Education", "DMDEDUC2", "same",
  "Less than high school / high school or GED / some college / college graduate."),
 ("Income", "INDFMPIR", "same", "Family income-to-poverty ratio, continuous."),
 ("Marital status", "DMDMARTL (earlier), DMDMARTZ", "DMDMARTZ",
  "Married or partnered / widowed, divorced or separated / never married."),
 ("Health insurance", "HIQ011", "same", "Any coverage, yes/no."),
 ("Food security", "FSDAD", "same", "Adult food-security category, low or very low vs higher."),
 ("Nativity", "DMDBORN4", "same", "Born in the USA, yes/no."),
 ("Depressive symptoms", "DPQ010-DPQ090", "same", "Sum of nine items, complete responses only; "
  "codes 7/9 set to missing. Interviewer-administered to March 2020, ACASI thereafter."),
 ("Serum folate vitamers", "LBXSF1SI (5-MTHF), LBXSF2SI (folic acid), LBXSF3SI "
  "(5-formyl-THF), LBXSF4SI (THF), LBXSF5SI (5,10-methenyl-THF), LBXSF6SI (MeFox), LBDFOTSI "
  "(total)", "same", "As released; below-LOD values already replaced by NCHS with "
  "LOD/sqrt(2). UMFA fraction = 100 x UMFA / (total folate + MeFox)."),
 ("Folate below-LOD comment flags", "LBDSF1LC to LBDSF6LC", "same",
  "0 = at or above the limit of detection, 1 = below; used for the detection rates in "
  "Supplementary Table S21."),
 ("Red-cell folate", "LBDRFOSI", "same", "Microbiologic assay."),
 ("Dietary folate", "DR1TFOLA (total food folate, ug/day)", "same",
  "Day 1 recall, reliable recalls only (DR1DRSTZ = 1); log(1+x)."),
 ("Dietary folic acid", "DR1TFA (synthetic folic acid, ug/day)", "same",
  "Day 1 recall, reliable recalls only; log(1+x)."),
 ("Total energy intake", "DR1TKCAL", "same", "Day 1 recall, reliable recalls only; log."),
 ("Supplemental folic acid", "DSQTFA (30-day supplement total, ug/day)", "same",
  "Missing set to zero for non-users; log(1+x). Folic-acid-containing supplement use defined "
  "as DSQTFA > 0."),
 ("Any supplement use", "DSD010", "same", "Yes/no in the past 30 days."),
 ("Recreational physical activity", "PAQ650/PAQ655/PAD660 (vigorous), "
  "PAQ665/PAQ670/PAD675 (moderate)", "PAD810Q/PAD810U/PAD820 (vigorous), "
  "PAD790Q/PAD790U/PAD800 (moderate)",
  "Both converted to moderate-equivalent minutes/week = moderate minutes + 2 x vigorous "
  "minutes. In the newer instrument the reported frequency is converted to times per week using "
  "the accompanying unit (day/week/month/year). Analysed as log(1+x)."),
 ("Sedentary time", "PAD680", "same", "Minutes/day; values >1300 set to missing."),
 ("Smoking", "SMQ020, SMQ040", "same", "Never / former / current."),
 ("Serum cotinine", "LBXCOT (COTNAL_G in 2011-2012)", "LBXCOT", "log-transformed."),
 ("Alcohol", "ALQ101/ALQ120Q (earlier), ALQ111/ALQ121", "ALQ111/ALQ121",
  "None / occasional / frequent."),
 ("Body mass index", "BMXBMI", "same", "kg/m2, continuous, per 5 units."),
 ("Sleep duration", "SLD010H (earlier), SLD012", "SLD012", "Hours per night; >24 set to missing."),
 ("Chronic conditions", "MCQ160B, MCQ160C, MCQ160F, MCQ160P, MCQ160L, MCQ160A, MCQ220",
  "same", "Count of heart failure, coronary disease, stroke, chronic lung disease, liver "
  "disease, arthritis and cancer."),
 ("Diabetes", "DIQ010", "same", "Doctor-diagnosed, yes/no (borderline coded as no)."),
 ("Renal function", "LBXSCR", "same", "2021 CKD-EPI creatinine equation, race-free."),
 ("C-reactive protein", "LBXHSCRP (2015-2016 onwards)", "LBXHSCRP",
  "High-sensitivity assay; not measured in 2011-2014. log-transformed."),
 ("Mean corpuscular volume", "LBXMCVSI", "same", "fL, continuous."),
 ("Vitamin B12", "LBDB12SI (2011-2014 only)", "not available",
  "Used only in the folate-trap interaction analysis."),
 ("Fasting time and session", "PHAFSTHR, PHDSESN", "PHAFSTHR, PHDSESNZ",
  "Hours since last food; morning / afternoon examination (no evening sessions occur in the "
  "analytic sample)."),
 ("Medication", "RXQ_RX drug identifiers matched to the NCHS therapeutic classification",
  "not available at drug level",
  "Antidepressant, anxiolytic and antipsychotic use; earlier cycles only."),
 ("Survey design", "WTMEC2YR (2011-2016), WTFOLPRP (2017-Mar 2020); SDMVSTRA, SDMVPSU",
  "WTPH2YR; SDMVSTRA, SDMVPSU",
  "Pooled weight = cycle weight x cycle years / 11.2."),
]
sup.append(md(pd.DataFrame(DICT, columns=["Construct", "2011-March 2020 variables",
                                          "August 2021-August 2023 variables", "Derivation"])))
sup.append("\n\n## Supplementary Table S21. All measurable serum folate species: detection, "
           "discovery and out-of-screen replication\n")
sup.append("**A. One species at a time, fully adjusted (Model 3)**\n")
sup.append(md(pd.DataFrame(A7["species"])))
sup.append("\nDetection is the percentage of analytic-sample specimens at or above the limit of "
           "detection, taken from the NCHS comment flag that accompanies each of the six species "
           "(LBDSF1LC to LBDSF6LC). Discovery is August 2021-August 2023, the cycle used for the "
           "exposure screen; replication is the four earlier cycles.\n")
sup.append("\n**B. All six species entered simultaneously (Model 3)**\n")
_rows = []
for _tag, _blk in A7["all_species"].items():
    for _r in _blk["rows"]:
        _rows.append({"Sample": _tag, "n": f"{_blk['n']:,}", "Folate species": _r["species"],
                      "OR (95% CI)": _r["OR"], "p": _r["p"]})
sup.append(md(pd.DataFrame(_rows)))
sup.append("\nTwo minor vitamers also showed a discovery-cycle association: 5,10-methenyl-THF "
           "(1.08, 1.02-1.15) and 5-formyl-THF (1.05, 1.01-1.09). Neither showed clear evidence of "
           "out-of-screen replication: 5,10-methenyl-THF was attenuated and imprecise in the "
           "earlier cycles (1.05, 0.99-1.12) and 5-formyl-THF reversed direction (0.94, "
           "0.88-1.00). Neither retained an association once all six species were entered "
           "together, whereas unmetabolised folic acid did. Both are at or above the limit of "
           "detection in only a small minority of specimens (5.3% and 4.2% respectively, against "
           "99.4% for unmetabolised folic acid and 100% for 5-MTHF), so roughly 95% of their "
           "values are the imputed below-detection fill and they are poor candidates for a stable "
           "exposure. In the six-species replication model the UMFA point estimate was unchanged "
           "at 1.11 but less precise (0.98-1.25), consistent with the loss of precision expected "
           "on adding heavily censored covariates rather than with attenuation of the estimate "
           "itself.\n")

sup.append("\n\n"
           "## Supplementary Figures\n")

sup.append("**Figure S1.** Participant flow and exclusions.\n")
sup.append("**Figure S2.** Subgroup (a) and sensitivity (b) analyses for unmetabolised folic acid "
           "and depressive symptoms.\n")
sup.append("**Figure S3.** Weighted geometric mean 5-methyltetrahydrofolate, unmetabolised folic "
           "acid and prevalence of depressive symptoms by survey cycle, 2011-2023. The dotted "
           "segment and shaded divider mark the change in dietary-supplement and PHQ-9 instruments "
           "before the final cycle.\n")
sup.append("**Figure S4.** Hypothesised biological pathways linking unmetabolised folic acid with "
           "depressive symptoms. Solid arrows and filled boxes mark steps established in previous "
           "biochemical and clinical work; dashed arrows and unfilled boxes mark hypothesised "
           "pathways, none of which was measured or tested in the present study. The schematic "
           "presents no result from this analysis.\n")
sup.append("**Figure S5.** Directed acyclic graph encoding the causal assumptions used to derive "
           "the reduced adjustment set (Model 5). Every edge is directed and the graph is acyclic "
           "by construction: where the real-world direction is arguable - body mass index, sleep, "
           "physical activity and depressive symptoms plausibly influence one another - the graph "
           "commits to the single direction assumed for identification rather than drawing a "
           "reciprocal pair. It is therefore not intended to represent every biologically "
           "plausible reciprocal relationship. Determinants of the exposure itself are separated "
           "from specimen factors that act on the measured concentration; dashed edges denote "
           "unmeasured variables.\n")

open("tables.md", "w", encoding="utf-8", newline="").write("\n".join(out))
open("supplement.md", "w", encoding="utf-8", newline="").write("\n".join(sup))
print("wrote tables.md and supplement.md")
