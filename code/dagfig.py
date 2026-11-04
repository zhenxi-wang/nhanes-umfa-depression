"""Supplementary Figure S5: the directed acyclic graph used to derive the reduced adjustment set.

Every edge is directed and the graph is acyclic by construction. Where the real-world direction is
arguable - body mass index, sleep, physical activity and depressive symptoms plausibly influence
one another - the graph commits to the single direction assumed for identification rather than
drawing a reciprocal pair, which would not be a DAG. Exposure determinants and specimen factors
are separated because they enter at different points: the former cause circulating UMFA, the
latter only the measured concentration.

One accent colour (the exposure); everything else is the neutral scale. Labels are fitted to their
boxes by figstyle.node.
"""
import os
from figstyle import setup, node, label, arrow, INK, MUT, LINE, FILL, ACC, ACC_F, DASH
from common import FIG

fig, ax = setup((7.2, 4.6), ylim=(0, 100))

# ---------------------------------------------------------------- exposure and outcome
node(ax, 18, 58, 30, 9.5, "Unmetabolised folic acid\n(measured serum concentration)",
     fc=ACC_F, ec=ACC, tc=ACC, bold=True, fs=7.0)
node(ax, 82, 58, 26, 9.5, "Depressive symptoms", fc=FILL, ec=INK, bold=True, fs=7.6, lw=1.3)
arrow(ax, 33.5, 58, 68.5, 58, color=ACC, lw=2.0, ms=11)
label(ax, 51, 61.6, "target causal effect", color=ACC, italic=True, fs=6.0, box=False)

# ---------------------------------------------------------------- common causes (adjusted)
node(ax, 50, 88, 78, 9.0,
     "Age  ·  sex  ·  race and Hispanic origin  ·  survey cycle  ·  education  ·  income\n"
     "food security  ·  nativity  ·  smoking  ·  renal function", fc=FILL, fs=6.4)
label(ax, 50, 95.6, "Common causes of both — retained in the reduced adjustment set (Model 5)",
      color=INK, fs=6.8, box=False)
arrow(ax, 13, 83.5, 15, 63.0, rad=0.16)
arrow(ax, 87, 83.5, 85, 63.0, rad=-0.16)

# ---------------------------------------------------------------- unmeasured
node(ax, 51, 74, 34, 7.5, "Unmeasured metabolic and genetic factors\n"
                          "(e.g. DHFR and MTHFR variation)", ls=DASH, fs=6.0)
arrow(ax, 34.5, 71.5, 24.0, 63.0, rad=0.12, ls=DASH)
arrow(ax, 67.5, 71.5, 78.0, 63.0, rad=-0.12, ls=DASH)

# ---------------------------------------------------------------- causes of the exposure
node(ax, 15, 31, 28, 12.0,
     "Dietary folic acid\nSupplemental folic acid\nSupplement use  ·  energy intake",
     fc=FILL, fs=6.2)
arrow(ax, 15, 37.0, 16.5, 53.0)
label(ax, 15, 19.4, "Dietary determinants of the exposure —\nretained in Model 5 to isolate the\n"
                    "component not explained by intake", fs=5.9, box=False)

node(ax, 47, 31, 22, 12.0, "Fasting time\nExamination session", fc=FILL, fs=6.2)
arrow(ax, 44, 37.0, 24.5, 53.0, rad=0.20)
label(ax, 47, 19.4, "Specimen factors — act on the\nmeasured concentration;\nretained in Model 5",
      fs=5.9, box=False)

# ---------------------------------------------------------------- downstream of the outcome
node(ax, 82, 31, 32, 12.0,
     "Body mass index  ·  sleep  ·  physical activity\nalcohol  ·  chronic conditions  ·  diabetes\n"
     "marital status  ·  health insurance", fc="white", fs=6.0)
arrow(ax, 82, 53.0, 82, 37.0, color=INK)
label(ax, 82, 19.4, "Assumed downstream of depressive symptoms —\nexcluded from Model 5, retained "
                    "in the extended Model 3;\nconditioning on them can open a collider path",
      fs=5.9, box=False)

# ---------------------------------------------------------------- key
ax.plot([2, 98], [12.5, 12.5], color=LINE, lw=0.8)
label(ax, 50, 7.0,
      "The graph encodes the causal assumptions used to identify the reduced adjustment set "
      "(Model 5); it is not intended to represent every biologically plausible reciprocal\n"
      "relationship. Where the direction is arguable the graph commits to one direction, so it "
      "remains acyclic. Dashed edges denote unmeasured variables.",
      fs=6.0, box=False)

for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(FIG, f"figS5_dag.{ext}"), bbox_inches="tight",
                facecolor="white", dpi=400)
print("wrote figS5_dag")
