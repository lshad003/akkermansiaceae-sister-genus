#!/usr/bin/env python3
# Orthogroup occupancy figure
# Output: results/figures/fig1_collapsed/Figure_orthogroup_polarity.pdf
import csv, collections, os, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
POL  = f"{BASE}/results/pangenome/orthogroup_polarity.tsv"
ENR  = f"{BASE}/results/pangenome/polarity_cog_enrichment.tsv"
OUTD = f"{BASE}/results/figures/fig1_collapsed"
os.makedirs(OUTD, exist_ok=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

rows = list(csv.DictReader(open(POL), delimiter="\t"))
cnt = collections.Counter(r["category"] for r in rows)
B = cnt["B_lost_at_gut_ancestor"]; A = cnt["A_lost_in_Akkermansia"]
S = cnt["shared_gut"]; C = cnt["C_Akkermansia_enriched"]; D = cnt["D_novel_specific"]
print("orthogroups:", len(rows), dict(cnt))

RED = "#D55E00"; GREY = "#999999"; BLUE = "#0072B2"
GREEN = "#009E73"; PURPLE = "#CC79A7"
COL = {"B_lost_at_gut_ancestor": RED, "shared_gut": GREY,
       "A_lost_in_Akkermansia": BLUE, "D_novel_specific": GREEN,
       "C_Akkermansia_enriched": PURPLE}
# small categories drawn last, larger and more opaque so single points are legible
ORDER = [("shared_gut", 13, 0.45), ("B_lost_at_gut_ancestor", 13, 0.45),
         ("C_Akkermansia_enriched", 20, 0.85), ("A_lost_in_Akkermansia", 20, 0.85),
         ("D_novel_specific", 20, 0.85)]

fig = plt.figure(figsize=(14.5, 7.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.08], wspace=0.10)

# ---------------- panel A: occupancy
ax = fig.add_subplot(gs[0, 0])
rng = np.random.default_rng(20260822)
for cat, sz, al in ORDER:
    sel = [r for r in rows if r["category"] == cat]
    if not sel:
        continue
    x = np.array([float(r["pct_novel"]) for r in sel]) + rng.normal(0, 1.05, len(sel))
    y = np.array([float(r["pct_akk"]) for r in sel]) + rng.normal(0, 0.85, len(sel))
    ax.scatter(x, y, s=sz, c=COL[cat], alpha=al, linewidths=0, zorder=2)
for v in (10, 70):
    ax.axvline(v, color="#cccccc", lw=0.8, ls="--", zorder=1)
    ax.axhline(v, color="#cccccc", lw=0.8, ls="--", zorder=1)

for x, y, t, c in [(88, 60, "retained across\nboth gut genera\n%s" % "{:,}".format(S), "#6a6a6a"),
                   (30, 4,  "lost at the shared\ngut ancestor\n%s" % "{:,}".format(B), RED),
                   (22, 92, "Akkermansia\nenriched\n%d" % C, PURPLE)]:
    ax.text(x, y, t, fontsize=8.5, color=c, ha="center", va="center", linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.35", fc="#ffffff", ec="none", alpha=0.82), zorder=4)

# the bottom-right corner holds two categories, so label each in its own colour
ax.text(88, 30, "present in the candidate genus,\nabsent from Akkermansia",
        fontsize=8.5, color="#4a4a4a", ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.30", fc="#ffffff", ec="none", alpha=0.82), zorder=4)
ax.text(88, 25, "%d lost on the Akkermansia branch" % A, fontsize=8.5, color=BLUE,
        ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.22", fc="#ffffff", ec="none", alpha=0.82), zorder=4)
ax.text(88, 20, "%d restricted to the candidate genus" % D, fontsize=8.5, color=GREEN,
        ha="center", va="center",
        bbox=dict(boxstyle="round,pad=0.22", fc="#ffffff", ec="none", alpha=0.82), zorder=4)

ax.set_xlabel("candidate-genus representatives carrying the orthogroup (%, n = 17)", fontsize=9)
ax.set_ylabel("Akkermansia representatives carrying the orthogroup (%, n = 24)", fontsize=9)
ax.set_xlim(-9, 106); ax.set_ylim(-9, 106)
ax.set_xticks([0, 10, 25, 50, 70, 100]); ax.set_yticks([0, 10, 25, 50, 70, 100])
for s_ in ("top", "right"):
    ax.spines[s_].set_visible(False)
ax.tick_params(labelsize=8)
ax.set_title("A   Orthogroup occupancy across %s polarized orthogroups"
             % "{:,}".format(len(rows)), fontsize=11, loc="left", pad=10)

handles = [Line2D([], [], marker="o", ls="", ms=5, color=GREY,   label="retained across both gut genera (%s)" % "{:,}".format(S)),
           Line2D([], [], marker="o", ls="", ms=5, color=RED,    label="lost at the shared gut ancestor (%s)" % "{:,}".format(B)),
           Line2D([], [], marker="o", ls="", ms=5, color=BLUE,   label="lost on the Akkermansia branch (%d)" % A),
           Line2D([], [], marker="o", ls="", ms=5, color=GREEN,  label="restricted to the candidate genus (%d)" % D),
           Line2D([], [], marker="o", ls="", ms=5, color=PURPLE, label="Akkermansia enriched (%d)" % C)]
ax.legend(handles=handles, fontsize=7.5, frameon=False, loc="upper left",
          bbox_to_anchor=(0.0, -0.075), ncol=2, handlelength=1.0, handletextpad=0.5)

ins = ax.inset_axes([0.16, 0.34, 0.42, 0.34])
CG = "#333333"; LW = 1.6
for xs, ys in [([0.02,0.24],[0.50,0.50]), ([0.24,0.24],[0.14,0.78]),
               ([0.24,0.60],[0.14,0.14]), ([0.24,0.46],[0.78,0.78]),
               ([0.46,0.46],[0.62,0.92]), ([0.46,0.74],[0.92,0.92]),
               ([0.46,0.74],[0.62,0.62])]:
    ins.plot(xs, ys, color=CG, lw=LW)
ins.text(0.77, 0.92, "candidate\ngenus", fontsize=7.5, va="center", color="#333333", linespacing=1.2)
ins.text(0.77, 0.62, "Akkermansia", fontsize=7.5, va="center", color="#333333", style="italic")
ins.text(0.62, 0.14, "free-living", fontsize=7.5, va="center", color="#6a6a6a")
ins.plot([0.35], [0.78], marker="o", ms=7, color="#ffffff", mec=RED, mew=1.7)
ins.text(0.35, 0.99, "%s lost" % "{:,}".format(B), fontsize=7.5, color=RED,
         ha="center", va="bottom")
ins.plot([0.60], [0.62], marker="o", ms=7, color="#ffffff", mec=BLUE, mew=1.7)
ins.text(0.60, 0.47, "%d lost" % A, fontsize=7.5, color=BLUE, ha="center", va="top")
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

keys = sorted(Bd, key=lambda k: -float(Bd[k]["background_pct"]))
names = {k: Bd[k]["COG_name"] for k in keys}
yp = np.arange(len(keys))[::-1]
h = 0.27

bg = [float(Bd[k]["background_pct"]) for k in keys]
bb = [float(Bd[k]["set_pct"]) for k in keys]
aa = [float(Ad[k]["set_pct"]) if k in Ad else 0.0 for k in keys]

ax2.barh(yp + h, bg, height=h, color=GREY, label="retained across both gut genera")
ax2.barh(yp,     bb, height=h, color=RED,  label="lost at the shared gut ancestor")
ax2.barh(yp - h, aa, height=h, color=BLUE, label="lost on the Akkermansia branch")

def star(q):
    q = float(q)
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""

for i, k in enumerate(keys):
    y = yp[i]
    sB = star(Bd[k]["q_BH"])
    if sB:
        ax2.text(bb[i] + 0.6, y, sB, va="center", fontsize=8, color=RED)
    if k in Ad and int(Ad[k]["set_n"]) > 0:
        sA = star(Ad[k]["q_BH"])
        if sA:
            ax2.text(aa[i] + 0.6, y - h, sA, va="center", fontsize=8, color=BLUE)

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

for ext in ("pdf", "png"):
    out = f"{OUTD}/Figure_orthogroup_polarity.{ext}"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    print("wrote", out)
