"""Figure 3: discovery, out-of-screen replication, pooled estimate and cycle-specific forest,
all on the primary three-species Model 3 specification."""
import json
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from common import FIG

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.linewidth": 0.8,
    "savefig.dpi": 400, "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42})
NAVY, RED, GREY = "#1f3864", "#c0392b", "#7f8c8d"

A6 = json.load(open("out_a6.json"))
dr = {r["sample"].split(":")[0]: r for r in A6["discovery_replication"]}
het = A6["heterogeneity_m3"]

rows = []
rows.append(("Discovery: August 2021-August 2023", dr["Discovery"], True))
rows.append(("Out-of-screen replication: 2011-March 2020", dr["Out-of-screen replication"], True))
rows.append(("Pooled: 2011-2023", dr["Pooled"], True))
rows.append((None, None, None))
for r in het["per_cycle"]:
    b, se = r["b"], r["se"]
    rows.append(("   " + r["cycle"], dict(_or=np.exp(b), _lo=np.exp(b - 1.96 * se),
                                          _hi=np.exp(b + 1.96 * se), n=""), False))

fig, ax = plt.subplots(figsize=(6.8, 4.1))
y = np.arange(len(rows))[::-1]
ax.axvline(1, color=GREY, lw=0.8, ls="--")
for i, (lab, r, is_main) in enumerate(rows):
    if lab is None:
        continue
    col = NAVY if is_main else GREY
    ms = 6.5 if is_main else 4.5
    mk = "D" if is_main else "o"
    ax.plot([r["_lo"], r["_hi"]], [y[i], y[i]], color=col, lw=2.0 if is_main else 1.3)
    ax.plot(r["_or"], y[i], mk, color=col, ms=ms, mec="white", mew=0.7)
ax.set_yticks([y[i] for i, (l, _, _) in enumerate(rows) if l is not None])
ax.set_yticklabels([l for l, _, _ in rows if l is not None], fontsize=8)
for t, (_, _, is_main) in zip(ax.get_yticklabels(),
                              [r for r in rows if r[0] is not None]):
    if is_main:
        t.set_fontweight("bold")
ax.set_xscale("log")
ax.set_xticks([0.9, 1.0, 1.1, 1.25, 1.5])
ax.set_xticklabels(["0.9", "1.0", "1.1", "1.25", "1.5"])
ax.set_xlim(0.82, 1.95)
ax.set_xlabel("Adjusted odds ratio for depressive symptoms (PHQ-9 $\\geq$ 10)\n"
              "per SD of log unmetabolised folic acid")
tr = ax.get_yaxis_transform()
for i, (lab, r, _) in enumerate(rows):
    if lab is None:
        continue
    ax.text(1.03, y[i], f"{r['_or']:.2f} ({r['_lo']:.2f}-{r['_hi']:.2f})", fontsize=7.2,
            va="center", ha="left", transform=tr, clip_on=False)
ax.text(0.01, y[4] + 0.55, "Individual survey cycles", fontsize=7.5, color=GREY,
        style="italic", transform=ax.get_yaxis_transform(), va="bottom")
yy = [y[i] for i, (l, _, _) in enumerate(rows) if l is not None and l.startswith("   ")]
ax.text(-0.005, (max(yy) + 0.6 - min(y)) / (max(y) - min(y) + 1) - 0.02, "",
        transform=ax.transAxes)
ax.axhline(y[3] + 0.0, color=GREY, lw=0.6, ls=":")
ax.text(0.5, -0.30,
        f"Between-cycle heterogeneity: Cochran $Q$ = {het['Q']:.2f} (df = 4), "
        f"$P$ = {het['p']:.2f}, $I^2$ = {het['I2']:.0f}%",
        transform=ax.transAxes, ha="center", fontsize=7.2, color=GREY)
for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(FIG, f"fig3_replication.{ext}"), bbox_inches="tight",
                facecolor="white")
plt.close(fig)
print("wrote fig3_replication")
