#!/usr/bin/env python3
# Orthogroup occupancy figure
# Source: ch3-chitin-evolution/scripts/fig_orthogroup_polarity.py
# Output: results/figures/fig1_collapsed/Figure_orthogroup_polarity.pdf
import csv, collections, os, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
POL  = f"{BASE}/results/pangenome/orthogroup_polarity.tsv"
ENR  = f"{BASE}/results/pangenome/polarity_cog_enrichment.tsv"
OUTD = f"{BASE}/results/figures/fig1_collapsed"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

rows = list(csv.DictReader(open(POL), delimiter="\t"))
cnt = collections.Counter(r["category"] for r in rows)
B = cnt["B_lost_at_gut_ancestor"]; A = cnt["A_lost_in_Akkermansia"]
S = cnt["shared_gut"]; C = cnt["C_Akkermansia_enriched"]; D = cnt["D_novel_specific"]
print("orthogroups:", len(rows))

COL = {"B_lost_at_gut_ancestor": "#b03030", "shared_gut": "#9a9a9a",
       "A_lost_in_Akkermansia": "#3f7f9f", "D_novel_specific": "#c96a3f",
       "C_Akkermansia_enriched": "#7a5ea8"}
ORDER = ["shared_gut", "B_lost_at_gut_ancestor", "A_lost_in_Akkermansia",
         "D_novel_specific", "C_Akkermansia_enriched"]

fig = plt.figure(figsize=(14.5, 7.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.08], wspace=0.10)

# ---------------- panel A: occupancy
ax = fig.add_subplot(gs[0, 0])
rng = np.random.default_rng(20260822)
for cat in ORDER:
    sel = [r for r in rows if r["category"] == cat]
    if not sel:
        continue
    x = np.array([float(r["pct_novel"]) for r in sel]) + rng.normal(0, 1.05, len(sel))
    y = np.array([float(r["pct_akk"]) for r in sel]) + rng.normal(0, 0.85, len(sel))
    ax.scatter(x, y, s=13, c=COL[cat], alpha=0.55, linewidths=0, zorder=2)
for v in (10, 70):
    ax.axvline(v, color="#cccccc", lw=0.8, ls="--", zorder=1)
    ax.axhline(v, color="#cccccc", lw=0.8, ls="--", zorder=1)

for x, y, t, c in [(88, 60, "retained across\nboth gut genera\n%s" % "{:,}".format(S), "#6a6a6a"),
                   (30, 4,  "lost at the shared\ngut ancestor\n%s" % "{:,}".format(B), "#b03030"),
                   (88, 22, "present in the candidate genus,\nabsent from Akkermansia\n"
                            "%d branch losses, %d restricted" % (A, D), "#3f7f9f"),
                   (22, 92, "Akkermansia\nenriched\n%d" % C, "#7a5ea8")]:
    ax.text(x, y, t, fontsize=8.5, color=c, ha="center", va="center", linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.35", fc="#ffffff", ec="none", alpha=0.82), zorder=4)

ax.set_xlabel("candidate-genus representatives carrying the orthogroup (%, n = 17)", fontsize=9)
ax.set_ylabel("Akkermansia representatives carrying the orthogroup (%, n = 24)", fontsize=9)
ax.set_xlim(-9, 106); ax.set_ylim(-9, 106)
ax.set_xticks([0, 10, 25, 50, 70, 100]); ax.set_yticks([0, 10, 25, 50, 70, 100])
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
ax.tick_params(labelsize=8)
ax.set_title("A   Orthogroup occupancy across %s polarized orthogroups"
             % "{:,}".format(len(rows)), fontsize=11, loc="left", pad=10)

ins = ax.inset_axes([0.16, 0.30, 0.42, 0.34])
CG = "#333333"; LW = 1.6
for xs, ys in [([0.02,0.24],[0.50,0.50]), ([0.24,0.24],[0.14,0.78]),
               ([0.24,0.60],[0.14,0.14]), ([0.24,0.46],[0.78,0.78]),
               ([0.46,0.46],[0.62,0.92]), ([0.46,0.74],[0.92,0.92]),
               ([0.46,0.74],[0.62,0.62])]:
    ins.plot(xs, ys, color=CG, lw=LW)
ins.text(0.77, 0.92, "candidate\ngenus", fontsize=7.5, va="center", color="#c96a3f", linespacing=1.2)
ins.text(0.77, 0.62, "Akkermansia", fontsize=7.5, va="center", color="#3f7f9f", style="italic")
ins.text(0.62, 0.14, "free-living", fontsize=7.5, va="center", color="#6a6a6a")
ins.plot([0.35], [0.78], marker="o", ms=7, color="#ffffff", mec="#b03030", mew=1.7)
ins.text(0.35, 0.99, "%s lost" % "{:,}".format(B), fontsize=7.5, color="#b03030",
         ha="center", va="bottom")
ins.plot([0.60], [0.62], marker="o", ms=7, color="#ffffff", mec="#3f7f9f", mew=1.7)
ins.text(0.60, 0.47, "%d lost" % A, fontsize=7.5, color="#3f7f9f", ha="center", va="top")
ins.set_xlim(0, 1.18); ins.set_ylim(0, 1.14); ins.axis("off")
ins.patch.set_facecolor("#ffffff"); ins.patch.set_alpha(0.90)

# ---------------- panel B: COG comparison
ax2 = fig.add_subplot(gs[0, 1])
if not os.path.exists(ENR):
    print("REFUSED: missing", ENR); sys.exit(1)
enr = list(csv.DictReader(open(ENR), delimiter="\t"))
by = collections.defaultdict(dict)
for r in enr:
    by[r["category"]][r["COG"]] = r

Bd = by.get("B_lost_at_gut_ancestor", {})
Ad = by.get("A_lost_in_Akkermansia", {})
if not Bd:
    print("REFUSED: no category B rows in the enrichment table"); sys.exit(1)

# order by the background percentage so the panel reads as a profile
keys = sorted(Bd, key=lambda k: -float(Bd[k]["background_pct"]))
names = {k: Bd[k]["COG_name"] for k in keys}
yp = np.arange(len(keys))[::-1]
h = 0.27

bg = [float(Bd[k]["background_pct"]) for k in keys]
bb = [float(Bd[k]["set_pct"]) for k in keys]
aa = [float(Ad[k]["set_pct"]) if k in Ad else 0.0 for k in keys]

ax2.barh(yp + h, bg, height=h, color="#9a9a9a", label="retained across both gut genera")
ax2.barh(yp,     bb, height=h, color="#b03030", label="lost at the shared gut ancestor")
ax2.barh(yp - h, aa, height=h, color="#3f7f9f", label="lost on the Akkermansia branch")

def star(q):
    q = float(q)
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""

for i, k in enumerate(keys):
    y = yp[i]
    sB = star(Bd[k]["q_BH"])
    if sB:
        ax2.text(bb[i] + 0.6, y, sB, va="center", fontsize=8, color="#b03030")
    if k in Ad and int(Ad[k]["set_n"]) > 0:
        sA = star(Ad[k]["q_BH"])
        if sA:
            ax2.text(aa[i] + 0.6, y - h, sA, va="center", fontsize=8, color="#3f7f9f")

ax2.set_yticks(yp)
ax2.set_yticklabels(["%s  %s" % (k, names[k]) for k in keys], fontsize=7.5)
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
ax2.set_xlabel("percentage of annotated orthogroups in the set", fontsize=9)
ax2.set_xlim(0, max(max(bg), max(bb), max(aa)) * 1.16)
for s_ in ("top", "left"):
    ax2.spines[s_].set_visible(False)
ax2.tick_params(axis="x", labelsize=8)
ax2.tick_params(axis="y", length=0)
ax2.legend(fontsize=8, frameon=False, loc="upper left",
           bbox_to_anchor=(0.0, -0.075), ncol=1, handlelength=1.6)
ax2.set_title("B   Functional composition of the two loss events", fontsize=11, loc="left", pad=10)

import textwrap
CAPTION = (
 "Panel A: points jittered; occupancy is quantised by panel size. Dashed lines mark the "
 "10 and 70 percent thresholds. Categories A and B additionally require occurrence in "
 "free-living genera and are read as branch events; C and D do not. "
 "Panel B: Fisher exact test against the retained set, Benjamini-Hochberg corrected. "
 "* q < 0.05, ** q < 0.01, *** q < 0.001; stars are omitted where a category is absent "
 "from a set. The background is the orthogroups retained across both gut genera, not the "
 "genome, so these are differences between lost and retained genes rather than enrichment "
 "against a neutral expectation.")

# wrap to the figure width: roughly 2.05 characters per inch per point of font size
FS = 7.5
chars = int(fig.get_figwidth() * 72.0 / (FS * 0.52))
fig.text(0.0, -0.02, textwrap.fill(CAPTION, chars),
         fontsize=FS, color="#666666", linespacing=1.55, va="top", ha="left")

for ext in ("pdf", "png"):
    out = f"{OUTD}/Figure_orthogroup_polarity.{ext}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("wrote", out)
