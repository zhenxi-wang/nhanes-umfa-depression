"""Main figures 1-2 and supplementary figures."""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from common import get, FIG

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "savefig.dpi": 400,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
})
NAVY, RED, GREY, TEAL, AMBER = "#1f3864", "#c0392b", "#7f8c8d", "#0e8a8a", "#d68910"

A1 = json.load(open("out_a1.json"))
A2 = json.load(open("out_a2.json"))
A3 = json.load(open("out_a3.json"))
A4 = json.load(open("out_a4.json"))
A6 = json.load(open("out_a6.json"))
d = get()
a = d[d.ins == 1]


def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print("wrote", name)


# =========================================================== FIGURE 1
fig = plt.figure(figsize=(7.1, 3.5))
gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 0.85], width_ratios=[1.25, 1.0],
                      hspace=0.12, wspace=0.32)
ax = fig.add_subplot(gs[0, 0])
axh = fig.add_subplot(gs[1, 0], sharex=ax)
c = pd.DataFrame(A2["spline"]["curve"])
ax.fill_between(c.umfa, c.lo, c.hi, color=NAVY, alpha=0.15, lw=0)
ax.plot(c.umfa, c.OR, color=NAVY, lw=1.9)
ax.axhline(1, color=GREY, lw=0.8, ls="--")
ax.set_xscale("log")
ax.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20])
ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
ax.set_ylabel("Adjusted odds ratio for depressive\nsymptoms (PHQ-9 $\\geq$ 10)")
ax.set_ylim(0.75, 2.3)
ax.tick_params(labelbottom=False)
ax.text(0.03, 0.94, "a", transform=ax.transAxes, fontweight="bold", fontsize=11, va="top")
ax.text(0.97, 0.06, f"$P$ for non-linearity = {A2['spline']['p_nonlinear']:.2f}",
        transform=ax.transAxes, ha="right", fontsize=7.5, color=GREY)
ax.annotate(f"reference:\n{A2['spline']['ref']:.2f} nmol/L", xy=(A2["spline"]["ref"], 1.0),
            xytext=(0.24, 1.42), fontsize=7, color=GREY,
            arrowprops=dict(arrowstyle="-", color=GREY, lw=0.7))
axh.hist(a.umfa.clip(upper=30), bins=np.logspace(np.log10(0.1), np.log10(30), 70),
         color=GREY, alpha=0.55, lw=0)
axh.set_xscale("log")
axh.set_yticks([])
axh.spines["left"].set_visible(False)
axh.set_xlim(0.15, 30)
axh.set_xticks([0.2, 0.5, 1, 2, 5, 10, 20])
axh.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10", "20"])
axh.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axh.set_xlabel("Serum unmetabolised folic acid (nmol/L, log scale)")

ax2 = fig.add_subplot(gs[:, 1])
pp = pd.DataFrame(A3["pred_prev"])
ax2.fill_between(pp.pct, pp.lo, pp.hi, color=NAVY, alpha=0.15, lw=0)
ax2.plot(pp.pct, pp.prev, "o-", color=NAVY, lw=1.8, ms=4)
ax2.set_ylabel("Model-adjusted prevalence of\ndepressive symptoms (%)")
ax2.set_xlabel("Percentile of serum unmetabolised folic acid\n(analytic-sample distribution)")
ax2.set_xticks([10, 25, 50, 75, 90, 99])
ax2.set_xlim(4, 105)
ax2.text(0.03, 0.94, "b", transform=ax2.transAxes, fontweight="bold", fontsize=11, va="top")
for _, r in pp.iterrows():
    if r.pct in (10, 99):
        ax2.annotate(f"{r.prev:.1f}%", (r.pct, r.prev), textcoords="offset points",
                     xytext=(0, 9 if r.pct == 99 else -14), ha="center", fontsize=7.5,
                     color=NAVY, fontweight="bold")
ax2.set_ylim(6.9, 16.6)
_rd = A6["risk_diff"]
ax2.text(0.5, -0.21,
         f"10th to 99th percentile difference:\n{_rd['diff']:.1f} percentage points "
         f"(95% CI, {_rd['lo']:.1f}-{_rd['hi']:.1f})",
         transform=ax2.transAxes, ha="center", va="top", fontsize=7.3, color=GREY,
         linespacing=1.35)
save(fig, "fig1_doseresponse")

# =========================================================== FIGURE 2
persd = pd.DataFrame(A1["persd"])
order = ["5-MTHF (per SD log)", "Total serum folate (per SD log)", "RBC folate (per SD log)",
         "MeFox (per SD log)", "Unmetabolised folic acid (per SD log)",
         "UMFA as % of total folate (per SD log)"]
short = {"5-MTHF (per SD log)": "5-Methyl-THF\n(bioactive folate)",
         "Total serum folate (per SD log)": "Total serum folate",
         "RBC folate (per SD log)": "Red-cell folate",
         "MeFox (per SD log)": "MeFox\n(oxidation product)",
         "Unmetabolised folic acid (per SD log)": "Unmetabolised\nfolic acid",
         "UMFA as % of total folate (per SD log)": "UMFA fraction\n(% of quantified folate)"}
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.6), gridspec_kw={"width_ratios": [1.35, 1]})
ax = axes[0]
mods = [1, 2, 3]
cols = {1: "#b6c4de", 2: "#5f7cb5", 3: NAVY}
ypos = np.arange(len(order))[::-1]
for k, mo in enumerate(mods):
    off = (k - 1) * 0.24
    for i, ex in enumerate(order):
        r = persd[(persd.exposure == ex) & (persd.model == mo)].iloc[0]
        y = ypos[i] + off
        hl = "highlight" in ex
        ax.plot([r._lo, r._hi], [y, y], color=cols[mo], lw=1.5, solid_capstyle="butt")
        ax.plot(r._or, y, "o", color=cols[mo], ms=4.2,
                mec="white", mew=0.6)
ax.axvline(1, color=GREY, lw=0.8, ls="--")
ax.set_yticks(ypos)
ax.set_yticklabels([short[o] for o in order], fontsize=8)
ax.set_xscale("log")
ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
ax.set_xticklabels(["0.7", "0.8", "0.9", "1.0", "1.1", "1.2", "1.3"])
ax.set_xlabel("Odds ratio for depressive symptoms per SD of log biomarker")
ax.axhspan(ypos[5] - 0.5, ypos[4] + 0.5, color=RED, alpha=0.06, zorder=0)
ax.set_ylim(ypos[5] - 0.55, ypos[0] + 0.55)
ax.legend(handles=[Line2D([], [], color=cols[m], lw=2,
                          label=["Model 1: demographic", "Model 2: + socioeconomic",
                                 "Model 3: + lifestyle, clinical, folate intake"][m - 1])
                   for m in mods], loc="upper center", frameon=False, fontsize=7,
          bbox_to_anchor=(0.42, -0.17), ncol=1)
ax.text(-0.62, 1.02, "a", transform=ax.transAxes, fontweight="bold", fontsize=11)

ax = axes[1]
q = pd.DataFrame(A1["prev_by_q"])
for ex, col, mk in [("5-MTHF", TEAL, "s"), ("UMFA", RED, "o")]:
    s = q[q.exposure == ex]
    ax.errorbar(s.Q, s.prev, yerr=[s.prev - s.lo, s.hi - s.prev], fmt=mk + "-", color=col,
                lw=1.7, ms=5, capsize=2.5, label=("5-Methyl-THF" if ex == "5-MTHF"
                                                  else "Unmetabolised folic acid"))
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"])
ax.set_xlabel("Quartile of biomarker (cycle-specific)")
ax.set_ylabel("Prevalence of depressive symptoms (%)")
ax.legend(frameon=False, fontsize=7.5, loc="upper center")
ax.set_ylim(5.5, 15.2)
ax.text(-0.28, 1.02, "b", transform=ax.transAxes, fontweight="bold", fontsize=11)
save(fig, "fig2_forest")

# =========================================================== SUPP: subgroups + sensitivity
sub = pd.DataFrame(A2["subgroups"])
sen = pd.DataFrame(A2["sensitivity"])
med = pd.DataFrame(A3["medication"])
rob = pd.DataFrame(A4["robust_extra"])
fig, axes = plt.subplots(2, 1, figsize=(6.9, 8.8),
                         gridspec_kw={"hspace": 0.30, "height_ratios": [15, 20]})


def forest(ax, labs, ors, los, his, ns, title, letter):
    y = np.arange(len(labs))[::-1]
    ax.axvline(1, color=GREY, lw=0.8, ls="--")
    for i in range(len(labs)):
        ax.plot([los[i], his[i]], [y[i], y[i]], color=NAVY, lw=1.4)
        ax.plot(ors[i], y[i], "o", color=NAVY, ms=4, mec="white", mew=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=7.3)
    ax.set_xscale("log")
    ax.set_xticks([0.8, 1.0, 1.25, 1.6])
    ax.set_xticklabels(["0.8", "1.0", "1.25", "1.6"])
    ax.set_xlim(0.7, 1.9)
    ax.set_xlabel("OR per SD of log UMFA")
    ax.set_title(title, fontsize=8.5, loc="left", fontweight="bold")
    tr = ax.get_yaxis_transform()
    for i in range(len(labs)):
        ax.text(1.04, y[i], f"{ors[i]:.2f} ({los[i]:.2f}-{his[i]:.2f})", fontsize=6.5,
                va="center", ha="left", transform=tr, clip_on=False)
    ax.text(-0.40, 1.05, letter, transform=ax.transAxes, fontweight="bold", fontsize=11)


PINT = {"Age 18-39": None, "Age 40-59": None, "Age >=60": A2["interaction_tests"]["age>=60"]["p"],
        "Men": None, "Women": A2["interaction_tests"]["female"]["p"],
        "Supplement user": None, "No supplements": A2["interaction_tests"]["supplement user"]["p"],
        "Never smoker": None, "Ever smoker": A2["interaction_tests"]["ever smoker"]["p"]}
forest(axes[0], list(sub.subgroup), list(sub._or), list(sub._lo), list(sub._hi), list(sub.n),
       "Subgroups", "a")
_yy = np.arange(len(sub))[::-1]
_tr = axes[0].get_yaxis_transform()
for _i, _lab in enumerate(sub.subgroup):
    if PINT.get(_lab) is not None:
        axes[0].text(1.30, _yy[_i], f"$P_{{interaction}}$ = {PINT[_lab]:.2f}", fontsize=6.3,
                     va="center", ha="left", color=GREY, transform=_tr, clip_on=False)
sn = pd.concat([sen, rob[~rob.analysis.isin(sen.analysis)], med.rename(
    columns={"analysis": "analysis"})[["analysis", "n", "_or", "_lo", "_hi", "_p"]]],
    ignore_index=True)
sn = sn.drop_duplicates("analysis")
forest(axes[1], list(sn.analysis), list(sn._or), list(sn._lo), list(sn._hi), list(sn.n),
       "Sensitivity analyses", "b")
save(fig, "figS2_sensitivity")

# =========================================================== SUPP: temporal trends
fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.0), gridspec_kw={"wspace": 0.45})
cy = ["2011-2012", "2013-2014", "2015-2016", "2017-Mar 2020", "Aug 2021-Aug 2023"]
xl = ["2011-\n2012", "2013-\n2014", "2015-\n2016", "2017-\nMar 2020", "Aug 2021-\nAug 2023"]
from common import design
des = design(d)
ins = (d.ins == 1).to_numpy()
gm_u, gm_m, dep = [], [], []
for c_ in cy:
    mk = ins & (d.cycle == c_).to_numpy()
    gm_u.append([np.exp(v) for v in des.subset(mk).mean(d.l_umfa.to_numpy())[:1] +
                 tuple(des.subset(mk).mean(d.l_umfa.to_numpy())[2:])])
    gm_m.append([np.exp(v) for v in des.subset(mk).mean(d.l_mthf.to_numpy())[:1] +
                 tuple(des.subset(mk).mean(d.l_mthf.to_numpy())[2:])])
    p = des.subset(mk).prop(d.dep.to_numpy())
    dep.append([100 * p["est"], 100 * p["lo"], 100 * p["hi"]])
for ax, vals, col, lab in [(axes[0], gm_m, TEAL, "5-Methyl-THF (nmol/L)"),
                           (axes[1], gm_u, RED, "Unmetabolised folic acid (nmol/L)"),
                           (axes[2], dep, NAVY, "Depressive symptoms, %")]:
    v = np.array(vals)
    ax.errorbar(range(4), v[:4, 0], yerr=[v[:4, 0] - v[:4, 1], v[:4, 2] - v[:4, 0]], fmt="o-",
                color=col, lw=1.6, ms=4.5, capsize=2.5)
    ax.errorbar([4], v[4:, 0], yerr=[v[4:, 0] - v[4:, 1], v[4:, 2] - v[4:, 0]], fmt="o",
                color=col, lw=1.6, ms=4.5, capsize=2.5)
    ax.plot([3, 4], v[3:5, 0], color=col, lw=1.2, ls=":")
    ax.axvspan(3.5, 3.6, color="#c9d1d9", alpha=0.7, lw=0)
    ax.set_xticks(range(5))
    ax.set_xticklabels(xl, fontsize=6.2, rotation=35, ha="right")
    ax.set_ylabel(lab, fontsize=7.5)
    ax.set_xlim(-0.4, 4.4)
save(fig, "figS3_trends")
