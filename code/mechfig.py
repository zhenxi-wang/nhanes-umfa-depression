"""Supplementary Figure S4: hypothesised biological pathways linking unmetabolised folic acid
with depressive symptoms.

Shows no study result: the observational estimates are in Table 2 and Figures 1-3, and repeating
them under this schematic invited the reading that they test the pathway. Steps established in
previous work are solid; everything this study could not observe is dashed.

Three columns on a fixed grid with wide gutters, so annotations sit in empty space rather than on
top of a box or an arrow. Every label is fitted to its box by figstyle.node.
"""
import os
from figstyle import (setup, node, label, arrow, heading,
                      INK, MUT, LINE, FILL, ACC, ACC_F, ACC2, ACC2_F, DASH)
from common import FIG

fig, ax = setup((7.2, 4.0), ylim=(15, 100))

C1, C2, C3 = 14.0, 50.0, 86.0     # column centres
BW = 24.0                          # standard box width
HEAD = 94.0

heading(ax, 2.0, 26.0, HEAD, "1", "Gut and liver")
heading(ax, 38.0, 62.0, HEAD, "2", "Circulation")
heading(ax, 74.0, 98.0, HEAD, "3", "Brain")

# ------------------------------------------------------------------ 1. gut and liver
node(ax, C1, 82.0, BW, 7.5, "Fortified grain foods and\nfolic acid supplements", fc=FILL)
arrow(ax, C1, 78.25, C1, 74.75)
node(ax, C1, 71.0, BW, 6.0, "Synthetic folic acid", fc=ACC_F, ec=ACC, tc=ACC, bold=True, fs=7.2)
arrow(ax, C1, 68.0, C1, 56.5)
label(ax, C1, 62.0, "DHFR - limited and variable\nactivity in humans", italic=True, fs=5.6)
node(ax, C1, 53.0, BW, 5.5, "Tetrahydrofolate")
arrow(ax, C1, 50.25, C1, 42.5)
label(ax, C1, 46.4, "MTHFR", italic=True, fs=5.6)
node(ax, C1, 39.0, BW, 5.5, "5-methyl-THF", fc=ACC2_F, ec=ACC2, tc=ACC2, bold=True, fs=7.0)

# ------------------------------------------------------------------ 2. circulation
node(ax, C2, 84.0, BW, 7.0, "Unmetabolised folic acid", fc=ACC_F, ec=ACC, tc=ACC, bold=True,
     fs=7.4)
label(ax, C2, 75.5, "median 0.64 nmol/L\nabout 2% of quantified folate\n"
                    "detected in more than 99% of adults", fs=5.5, box=False)
node(ax, C2, 58.0, BW, 7.0, "5-methyl-THF\n(bioactive folate)", fc=ACC2_F, ec=ACC2, tc=ACC2,
     bold=True, fs=7.0)
label(ax, C2, 51.6, "about 90% of total serum folate", fs=5.7, box=False)
node(ax, C2, 42.0, BW + 2, 9.0,
     "Total serum folate is the sum of the two,\nbut is dominated by 5-methyl-THF and\n"
     "does not resolve them separately", fc=FILL, fs=5.9)

# spill-over from column 1 into column 2
arrow(ax, 26.5, 71.0, 37.5, 82.0, color=ACC, lw=1.8, rad=-0.18, ms=10)
label(ax, 32.2, 88.8, "intake exceeds\nreduction capacity", color=ACC, italic=True, fs=5.4)
arrow(ax, 26.5, 39.0, 37.5, 55.5, color=ACC2, lw=1.5, rad=0.18)

# ------------------------------------------------------------------ 3. brain (hypothesised)
node(ax, C3, 82.0, BW, 7.5, "Reduced delivery of bioactive\nfolate to the central nervous system",
     ec=ACC, ls=DASH, fs=6.2)
node(ax, 80.0, 65.0, 11.0, 8.5, "Impaired\none-carbon\nmethylation", ec=MUT, ls=DASH, fs=5.8)
node(ax, 92.0, 65.0, 11.0, 8.5, "Impaired BH$_4$\nregeneration", ec=MUT, ls=DASH, fs=5.8)
arrow(ax, 83.0, 78.25, 80.0, 69.5, ls=DASH)
arrow(ax, 89.0, 78.25, 92.0, 69.5, ls=DASH)
node(ax, C3, 50.0, BW, 7.0, "Altered monoamine synthesis\nand availability", ec=MUT, ls=DASH,
     fs=6.2)
arrow(ax, 80.0, 60.75, 83.0, 53.5, ls=DASH)
arrow(ax, 92.0, 60.75, 89.0, 53.5, ls=DASH)
label(ax, C3, 44.6, "parallel routes, not a sequential chain", fs=5.7, box=False)

# column 2 into column 3
arrow(ax, 62.5, 84.0, 73.5, 83.0, color=ACC, lw=1.4, ls=DASH)
arrow(ax, 62.5, 58.0, 73.5, 78.5, color=ACC2, lw=1.4, rad=-0.16, ls=DASH)
label(ax, C2, 66.4, "UMFA competes with 5-methyl-THF for folate\n"
                    "receptor alpha at the choroid plexus", color=ACC, fs=5.6, box=False)

# ------------------------------------------------------------------ key
ax.plot([2, 98], [27.0, 27.0], color=LINE, lw=0.8)
label(ax, 2.0, 21.0,
      "Solid lines and filled boxes: steps established in previous biochemical and clinical work.   "
      "Dashed lines and open boxes: hypothesised pathways, none of\nwhich was measured or tested "
      "in the present study.   The schematic reports no result from this analysis; the "
      "observational estimates are in Table 2 and Figures 1-3.",
      fs=6.0, ha="left", box=False)

for ext in ("png", "pdf", "svg"):
    fig.savefig(os.path.join(FIG, f"figS4_mechanism.{ext}"), bbox_inches="tight",
                facecolor="white", dpi=400)
print("wrote figS4_mechanism")
