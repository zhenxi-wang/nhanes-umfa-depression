"""Figure S1: participant flow, with exclusion counts computed from the data."""
import json
import os
import numpy as np
import textwrap
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from common import get, FIG

d = get()
n_all = len(d)
s = d.copy()
step = []
m = pd.Series(True, index=s.index)
step.append(("Participants in NHANES 2011-2012, 2013-2014, 2015-2016,\n"
             "2017-March 2020 and August 2021-August 2023", int(m.sum()), None))
ex1 = m & (s.age < 18)
m = m & ~ex1
step.append(("Aged 18 years or older", int(m.sum()), ("Aged under 18 years", int(ex1.sum()))))
ex2 = m & (s.preg == 1)
m = m & ~ex2
step.append(("Not pregnant at examination", int(m.sum()), ("Pregnant", int(ex2.sum()))))
ex3 = m & (s.wraw.fillna(0) <= 0)
m = m & ~ex3
step.append(("Non-zero folate-subsample examination weight", int(m.sum()),
             ("No examination / zero analytic weight", int(ex3.sum()))))
ex4 = m & (s.mthf.isna() | s.umfa.isna())
m = m & ~ex4
step.append(("Serum folate forms measured by LC-MS/MS", int(m.sum()),
             ("5-MTHF or folic acid not quantified", int(ex4.sum()))))
ex5 = m & s.phq9.isna()
m = m & ~ex5
step.append(("Complete PHQ-9 (all nine items)", int(m.sum()),
             ("Incomplete PHQ-9", int(ex5.sum()))))
assert int(m.sum()) == int((d.ins == 1).sum()), (int(m.sum()), int((d.ins == 1).sum()))

by_cycle = d[d.ins == 1].groupby("cycle", sort=False).size()
print(by_cycle)
json.dump({"steps": [(a, b, c) for a, b, c in step],
           "by_cycle": {k: int(v) for k, v in by_cycle.items()}},
          open("flow.json", "w"), indent=1)

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
INK, MUT, EDGE, BLUEL, BLUE = "#1a2433", "#5c6b7a", "#c9d1d9", "#dceaf7", "#2f6fb0"
fig, ax = plt.subplots(figsize=(7.0, 8.2))
ax.set_xlim(0, 100)
ax.set_ylim(-15, 100)
ax.axis("off")
ytop, gap, bh = 93.0, 14.0, 8.6
for i, (lab, n, ex) in enumerate(step):
    y = ytop - i * gap
    fc = BLUEL if i in (0, len(step) - 1) else "white"
    ax.add_patch(FancyBboxPatch((14, y - bh), 52, bh, boxstyle="round,pad=0.5,rounding_size=1.2",
                                fc=fc, ec=BLUE if i in (0, len(step) - 1) else EDGE, lw=1.0))
    ax.text(40, y - bh / 2 + 1.3, lab, ha="center", va="center", fontsize=6.8, color=INK,
            linespacing=1.3)
    ax.text(40, y - bh + 2.0, f"n = {n:,}", ha="center", va="center", fontsize=7.6,
            fontweight="bold", color=BLUE)
    if i < len(step) - 1:
        ax.add_patch(FancyArrowPatch((40, y - bh), (40, y - gap), arrowstyle="-|>",
                                     color=MUT, lw=1.1, mutation_scale=9))
    if ex is not None:
        # The exclusion reasons vary in length, so wrap them to the box width and grow the
        # box to fit rather than letting long text run outside its border.
        lines = textwrap.wrap(f"Excluded: {ex[0]}", width=30) + [f"(n = {ex[1]:,})"]
        ebh = max(7.6, 2.55 * len(lines) + 2.5)
        ymid = y - bh + 3.5
        ax.plot([66, 68], [ymid, ymid], color=MUT, lw=0.8)
        ax.add_patch(FancyArrowPatch((68, ymid), (70, ymid),
                                     arrowstyle="-|>", color=MUT, lw=0.8, mutation_scale=8))
        ax.add_patch(FancyBboxPatch((70, ymid - ebh / 2), 29, ebh,
                                    boxstyle="round,pad=0.4,rounding_size=1.0",
                                    fc="#f7f8fa", ec=EDGE, lw=0.8))
        ax.text(84.5, ymid, "\n".join(lines), ha="center",
                va="center", fontsize=6.0, color=MUT, linespacing=1.3)
yb = ytop - (len(step) - 1) * gap - bh - 7
ax.add_patch(FancyBboxPatch((14, yb - 16.6), 52, 16.6,
                            boxstyle="round,pad=0.5,rounding_size=1.2", fc="white",
                            ec=EDGE, lw=1.0))
ax.text(40, yb - 1.8, "Analytic sample by survey cycle", ha="center", va="top",
        fontsize=6.8, fontweight="bold", color=INK)
txt = "\n".join(f"{k}:  n = {v:,}" for k, v in by_cycle.items())
ax.text(40, yb - 9.7, txt, ha="center", va="center", fontsize=6.3, color=MUT,
        linespacing=1.55)
ax.add_patch(FancyArrowPatch((40, yb + 6), (40, yb), arrowstyle="-|>", color=MUT, lw=1.1,
                             mutation_scale=9))
for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(FIG, f"figS1_flow.{ext}"), bbox_inches="tight",
                facecolor="white", dpi=400)
print("wrote figS1_flow")
