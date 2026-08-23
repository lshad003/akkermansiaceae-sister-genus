#!/usr/bin/env python3
# Orthogroup occupancy figure
# Source: ch3-chitin-evolution/scripts/fig_orthogroup_polarity.py
# Output: results/figures/fig1_collapsed/Figure_orthogroup_polarity.pdf
import csv, collections, os, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
POL  = f"{BASE}/results/pangenome/orthogroup_polarity.tsv"
OUTD = f"{BASE}/results/figures/fig1_collapsed"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

rows = list(csv.DictReader(open(POL), delimiter="\t"))
cnt = collections.Counter(r["category"] for r in rows)
print("rows:", len(rows))
for k, v in cnt.most_common():
    print("  %-26s %d" % (k, v))
nf = [int(float(r["n_free_genera"])) for r in rows]
print("n_free_genera: min %d max %d" % (min(nf), max(nf)))

B = cnt["B_lost_at_gut_ancestor"]; A = cnt["A_lost_in_Akkermansia"]
S = cnt["shared_gut"]; C = cnt["C_Akkermansia_enriched"]; D = cnt["D_novel_specific"]

COL = {"B_lost_at_gut_ancestor": "#b03030",
       "shared_gut":             "#9a9a9a",
       "A_lost_in_Akkermansia":  "#3f7f9f",
       "D_novel_specific":       "#c96a3f",
       "C_Akkermansia_enriched": "#7a5ea8"}
ORDER = ["shared_gut", "B_lost_at_gut_ancestor", "A_lost_in_Akkermansia",
         "D_novel_specific", "C_Akkermansia_enriched"]

fig, ax = plt.subplots(figsize=(8.4, 7.4))
rng = np.random.default_rng(20260822)

for cat in ORDER:
    sel = [r for r in rows if r["category"] == cat]
    if not sel:
        continue
    x = np.array([float(r["pct_novel"]) for r in sel])
    y = np.array([float(r["pct_akk"]) for r in sel])
    # jitter: values are quantised by panel size (17 and 24 representatives)
    x = x + rng.normal(0, 1.05, len(x))
    y = y + rng.normal(0, 0.85, len(y))
    ax.scatter(x, y, s=13, c=COL[cat], alpha=0.55, linewidths=0, zorder=2)

for v in (10, 70):
    ax.axvline(v, color="#cccccc", lw=0.8, ls="--", zorder=1)
    ax.axhline(v, color="#cccccc", lw=0.8, ls="--", zorder=1)

ax.text(3, 103, "absent from the\ncandidate genus", fontsize=8.5, color="#777777",
        ha="left", va="top")
ax.text(85, -4, "present in the\ncandidate genus", fontsize=8.5, color="#777777",
        ha="center", va="top")

ann = [(88, 60, "retained across\nboth gut genera\n%s" % "{:,}".format(S), "#6a6a6a"),
       (30, 4,  "lost at the shared\ngut ancestor\n%s" % "{:,}".format(B), "#b03030"),
       (88, 22, "present in the candidate genus,\nabsent from Akkermansia\n%d branch losses, %d restricted" % (A, D), "#3f7f9f"),
       (22, 92, "Akkermansia\nenriched\n%d" % C, "#7a5ea8")]
for x, y, t, c in ann:
    ax.text(x, y, t, fontsize=9, color=c, ha="center", va="center",
            linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.35", fc="#ffffff", ec="none", alpha=0.82),
            zorder=4)

ax.set_xlabel("percentage of candidate-genus representatives carrying the orthogroup (n = 17)",
              fontsize=9.5)
ax.set_ylabel("percentage of Akkermansia representatives carrying the orthogroup (n = 24)",
              fontsize=9.5)
ax.set_xlim(-6, 106); ax.set_ylim(-6, 106)
ax.set_xticks([0, 10, 25, 50, 70, 100]); ax.set_yticks([0, 10, 25, 50, 70, 100])
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(labelsize=8.5)
ax.set_title("Orthogroup occupancy across %s polarized orthogroups"
             % "{:,}".format(len(rows)), fontsize=11, loc="left", pad=12)

handles = [Line2D([], [], marker="o", ls="none", ms=6, color=COL[c],
                  label={"shared_gut": "retained across both gut genera",
                         "B_lost_at_gut_ancestor": "lost at the shared gut ancestor",
                         "A_lost_in_Akkermansia": "lost on the Akkermansia branch",
                         "D_novel_specific": "candidate-genus restricted",
                         "C_Akkermansia_enriched": "Akkermansia enriched"}[c])
           for c in ORDER]
ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.13),
          ncol=2, fontsize=8.5, frameon=False, handletextpad=0.4)

# ---- inset: which branch the two loss events sit on
ins = ax.inset_axes([0.16, 0.30, 0.42, 0.34])
CG = "#333333"; LW = 1.6
ins.plot([0.02, 0.24], [0.50, 0.50], color=CG, lw=LW)
ins.plot([0.24, 0.24], [0.14, 0.78], color=CG, lw=LW)
ins.plot([0.24, 0.60], [0.14, 0.14], color=CG, lw=LW)
ins.plot([0.24, 0.46], [0.78, 0.78], color=CG, lw=LW)
ins.plot([0.46, 0.46], [0.62, 0.92], color=CG, lw=LW)
ins.plot([0.46, 0.74], [0.92, 0.92], color=CG, lw=LW)
ins.plot([0.46, 0.74], [0.62, 0.62], color=CG, lw=LW)
ins.text(0.77, 0.92, "candidate\ngenus", fontsize=7.5, va="center", color="#c96a3f",
         linespacing=1.2)
ins.text(0.77, 0.62, "Akkermansia", fontsize=7.5, va="center", color="#3f7f9f",
         style="italic")
ins.text(0.62, 0.14, "free-living", fontsize=7.5, va="center", color="#6a6a6a")
ins.plot([0.35], [0.78], marker="o", ms=7, color="#ffffff", mec="#b03030", mew=1.7)
ins.text(0.35, 0.99, "%s lost" % "{:,}".format(B), fontsize=7.5, color="#b03030",
         ha="center", va="bottom")
ins.plot([0.60], [0.62], marker="o", ms=7, color="#ffffff", mec="#3f7f9f", mew=1.7)
ins.text(0.60, 0.47, "%d lost" % A, fontsize=7.5, color="#3f7f9f",
         ha="center", va="top")
ins.set_xlim(0, 1.18); ins.set_ylim(0, 1.14); ins.axis("off")
ins.patch.set_facecolor("#ffffff"); ins.patch.set_alpha(0.90)

fig.text(0.02, -0.10,
         "Points are jittered; occupancy is quantised by panel size. Dashed lines mark the "
         "10 and 70 percent thresholds. Categories A and B additionally require occurrence "
         "in free-living genera and are read as branch events; C and D carry no such "
         "requirement and are distributional only.",
         fontsize=7.5, color="#555555", wrap=True)

for ext in ("pdf", "png"):
    out = f"{OUTD}/Figure_orthogroup_polarity.{ext}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("wrote", out)
