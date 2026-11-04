"""Shared drawing helpers for the two schematic supplementary figures.

The important one is `node`: it measures the rendered text and shrinks the font until the label
actually fits inside its box. Hand-tuned font sizes drifted out of their borders every time a
label was reworded, which is what made the earlier versions look untidy.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# One neutral scale plus two accents. Everything that is not the exposure or the outcome is grey.
INK = "#252d38"
MUT = "#69757f"
LINE = "#b9c2cb"
FILL = "#f4f6f8"
ACC = "#9e4038"     # unmetabolised folic acid
ACC_F = "#f7eae8"
ACC2 = "#37716c"    # 5-methyltetrahydrofolate
ACC2_F = "#e6f0ef"
DASH = (0, (3.0, 2.2))


def measure(ax, t):
    """Text extent in data units."""
    fig = ax.figure
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    p0 = ax.transData.inverted().transform((bb.x0, bb.y0))
    p1 = ax.transData.inverted().transform((bb.x1, bb.y1))
    return abs(p1[0] - p0[0]), abs(p1[1] - p0[1])


def node(ax, x, y, w, h, text, *, fc="white", ec=None, fs=6.6, bold=False, tc=None,
         ls="-", lw=1.0, zorder=4, pad=2.2):
    """Rounded box with text guaranteed to fit inside it."""
    ec = MUT if ec is None else ec
    tc = INK if tc is None else tc
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0,rounding_size=0.9",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=zorder))
    t = ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc,
                zorder=zorder + 1, fontweight="bold" if bold else "normal", linespacing=1.32)
    for _ in range(24):
        tw, th = measure(ax, t)
        if tw <= w - pad and th <= h - 1.2:
            break
        fs -= 0.25
        t.set_fontsize(fs)
    return t


def label(ax, x, y, text, *, fs=5.7, color=None, italic=False, box=True, ha="center",
          va="center", zorder=6):
    """Free-standing annotation, optionally masked so it never sits on top of a line."""
    return ax.text(x, y, text, ha=ha, va=va, fontsize=fs,
                   color=MUT if color is None else color, zorder=zorder, linespacing=1.35,
                   style="italic" if italic else "normal",
                   bbox=dict(fc="white", ec="none", pad=1.6) if box else None)


def arrow(ax, x1, y1, x2, y2, *, color=None, lw=1.1, rad=0.0, ls="-", ms=8, zorder=3):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 color=MUT if color is None else color, lw=lw, linestyle=ls,
                                 mutation_scale=ms, zorder=zorder, shrinkA=1, shrinkB=1,
                                 connectionstyle=f"arc3,rad={rad}"))


def blunt(ax, x, y, color=ACC, half=1.9, lw=1.6, zorder=8):
    """Inhibition symbol: a short bar, drawn on an opaque white pad so it reads cleanly."""
    ax.plot([x - half, x + half], [y, y], color=color, lw=lw, zorder=zorder,
            solid_capstyle="butt")


def heading(ax, x0, x1, y, num, title):
    ax.text(x0, y, f"{num}", fontsize=7.2, fontweight="bold", color="white", ha="center",
            va="center", zorder=6,
            bbox=dict(boxstyle="circle,pad=0.34", fc=INK, ec="none"))
    ax.text(x0 + 3.4, y, title, fontsize=8.4, fontweight="bold", color=INK, ha="left",
            va="center", zorder=6)
    ax.plot([x0 - 2.4, x1], [y - 3.6, y - 3.6], color=LINE, lw=0.9, zorder=2)


def setup(figsize, xlim=(0, 100), ylim=(0, 100)):
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42,
                         "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.axis("off")
    return fig, ax
