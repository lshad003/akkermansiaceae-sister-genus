#!/usr/bin/env python3
# Placement figure rendered from the collapsed tree
# Source: ch3-chitin-evolution/scripts/fig1_render.py
# Output: results/figures/fig1_collapsed/Figure1_placement.pdf
import os, re, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
NWK  = f"{BASE}/results/figures/fig1_collapsed/fig1_iqtree.nwk"
COL  = f"{BASE}/results/figures/fig1_collapsed/fig1_iqtree_colors.txt"
OUTD = f"{BASE}/results/figures/fig1_collapsed"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from ete3 import Tree

for p in (NWK, COL):
    if not os.path.exists(p):
        print("MISSING:", p); sys.exit(1)

# ---- colours, verified mapping
C_AKK_OTHER = "#b8b6ac"
C_AKK_AMPH  = "#d9a441"
C_NOV_UHM   = "#4a7ebb"
C_NOV_EHI   = "#c96a3f"
C_NOV_MIX   = "#7a5ea8"
C_FREE      = "#5c5c5c"
C_OUT       = "#000000"

LEGEND = [
    (C_NOV_EHI,   "candidate genus, EHI species clusters"),
    (C_NOV_UHM,   "candidate genus, UHM species clusters"),
    (C_NOV_MIX,   "candidate genus, mixed species cluster"),
    (C_AKK_AMPH,  "Akkermansia, amphibian-derived clades"),
    (C_AKK_OTHER, "Akkermansia, other clades"),
    (C_FREE,      "free-living genera"),
    (C_OUT,       "outgroup families"),
]

colour = {}
for line in open(COL):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 3 and f[1] == "range":
        colour[f[0]] = f[2]

def tipcolour(n):
    if n in colour:
        return colour[n]
    return C_OUT if n.startswith("OUTGROUP") else C_FREE

t = Tree(NWK, format=1)
t.ladderize()
leaves = t.get_leaves()
print("tips:", len(leaves))
missing = [l.name for l in leaves if l.name not in colour
           and not l.name.startswith("OUTGROUP")]
print("uncoloured (drawn as free-living):", len(missing))

for i, l in enumerate(leaves):
    l.y = i
for node in t.traverse("postorder"):
    if not node.is_leaf():
        node.y = sum(c.y for c in node.children) / float(len(node.children))

t.x = 0.0
for node in t.traverse("preorder"):
    if node.up is not None:
        node.x = node.up.x + (node.dist or 0.0)

for node in t.traverse():
    node.x_draw = node.x
    node.broken = False
truncated = []
xmax = max(n.x for n in t.traverse())

def support_pair(node):
    m = re.match(r"^([0-9.]+)/([0-9.]+)$", (node.name or "").strip())
    return (float(m.group(1)), float(m.group(2))) if m else None

nov = [l for l in leaves if l.name.startswith("NOVEL")]
akk = [l for l in leaves if l.name.startswith("Akkermansia")]
key = {"candidate genus": t.get_common_ancestor(nov),
       "Akkermansia": t.get_common_ancestor(akk),
       "shared ancestor": t.get_common_ancestor(nov + akk)}
print("key nodes:", ", ".join("%s=%s" % (k, v.name) for k, v in key.items()))

fig, ax = plt.subplots(figsize=(12, 13))

for node in t.traverse():
    if node.up is not None:
        ax.plot([node.up.x_draw, node.x_draw], [node.y, node.y],
                color="#333333", lw=0.9, zorder=1)
    if not node.is_leaf():
        ys = [c.y for c in node.children]
        ax.plot([node.x_draw, node.x_draw], [min(ys), max(ys)],
                color="#333333", lw=0.9, zorder=1)

# break marker on truncated branches
for l in leaves:
    if getattr(l, "broken", False):
        bx = l.x_draw - xmax * 0.012
        for off in (-xmax * 0.006, xmax * 0.006):
            ax.plot([bx + off - xmax * 0.004, bx + off + xmax * 0.004],
                    [l.y - 0.42, l.y + 0.42], color="#ffffff", lw=2.4, zorder=2)
            ax.plot([bx + off - xmax * 0.004, bx + off + xmax * 0.004],
                    [l.y - 0.42, l.y + 0.42], color="#333333", lw=0.8, zorder=3)

keyset = {id(v) for v in key.values()}
for node in t.traverse():
    if node.is_leaf() or node.up is None or id(node) in keyset:
        continue
    sp = support_pair(node)
    if sp and sp[0] >= 80 and sp[1] >= 95:
        ax.plot([node.x_draw], [node.y], marker="o", ms=3.4,
                color="#333333", mec="none", zorder=4)

for lab, node in key.items():
    sp = support_pair(node)
    ax.plot([node.x_draw], [node.y], marker="o", ms=9, color="#ffffff",
            mec="#b03030", mew=1.9, zorder=5)
    ax.annotate("%s\n%s" % (lab, "%g/%g" % sp if sp else "?"),
                (node.x_draw, node.y), xytext=(-13, 0), textcoords="offset points",
                ha="right", va="center", fontsize=8.5, color="#b03030",
                linespacing=1.3, zorder=6)

for l in leaves:
    ax.plot([l.x_draw, xmax * 1.015], [l.y, l.y], color="#d5d5d5",
            lw=0.5, ls=":", zorder=0)
    ax.text(xmax * 1.03, l.y, l.name.replace("_", " "), va="center", ha="left",
            fontsize=7.5, color=tipcolour(l.name))

ybase = -2.0
sc = xmax * 0.2
ax.plot([0, sc], [ybase, ybase], color="#333333", lw=1.5)
ax.text(sc / 2.0, ybase - 0.7, "%.2f substitutions per site" % sc,
        ha="center", va="top", fontsize=8, color="#333333")

handles = [Line2D([], [], color=c, lw=3.2, label=lab) for c, lab in LEGEND]
handles += [
    Line2D([], [], marker="o", ls="none", ms=3.4, color="#333333",
           label="SH-aLRT at least 80 and UFBoot at least 95"),
    Line2D([], [], marker="o", ls="none", ms=9, color="#ffffff", mec="#b03030",
           mew=1.9, label="nodes central to the placement"),
]
ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=False,
          bbox_to_anchor=(0.0, -0.035), ncol=2, handlelength=1.6,
          labelspacing=0.55, columnspacing=2.2)

if truncated:
    ax.text(0.0, ybase - 2.6,
            "Outgroup branches are truncated at the break marks; all other "
            "branch lengths are drawn to scale.",
            ha="left", va="top", fontsize=7.5, color="#555555")

ax.set_xlim(-xmax * 0.03, xmax * 1.60)
ax.set_ylim(ybase - 1.5, len(leaves) + 0.5)
ax.axis("off")
fig.subplots_adjust(bottom=0.16)

for ext in ("pdf", "png"):
    out = f"{OUTD}/Figure1_placement.{ext}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("wrote", out)
