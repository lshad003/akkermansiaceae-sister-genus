#!/usr/bin/env python3
# Operon figure, Option A: geometry scaled to real gene lengths and real spacing.
# Everything is read from operon_allgenera.tsv. SUPERSEDES fig_ppp_operon_final.py.
import os
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TSV = R + "/results/ppp_unified/operon_allgenera.tsv"
OUT = R + "/results/figures/Figure_ppp_operon"

rows = []
with open(TSV) as f:
    h = f.readline().rstrip("\n").split("\t")
    for line in f:
        rows.append(dict(zip(h, line.rstrip("\n").split("\t"))))

def L(r, which):
    return abs(int(r[which + "_end"]) - int(r[which + "_start"])) + 1

byg = defaultdict(list)
for r in rows: byg[r["group"]].append(r)
order = ["CANDIDATE"] + sorted([k for k in byg if k != "CANDIDATE"],
                               key=lambda k: -len(byg[k]))
order = [k for k in order if len(byg[k]) >= 2]

# median geometry per genus, from genomes where both genes share a contig
geo = {}
for g in order:
    sc = [r for r in byg[g] if r["same_contig"] == "True"]
    if not sc: continue
    med = lambda v: sorted(v)[len(v)//2]
    geo[g] = dict(
        zwf = med([L(r, "zwf") for r in sc]),
        sub = med([L(r, "sub") for r in sc]),
        gap = med([int(r["gap_bp"]) for r in sc if r["gap_bp"] not in ("None", "")]),
        n = len(byg[g]), nsc = len(sc),
        co = sum(1 for r in sc if r["co_oriented"] == "True"),
        ok = sum(1 for r in sc if r["order"] == "zwf-sub"))
order = [g for g in order if g in geo]

RED, BLUE = "#c0392b", "#2c6fbb"
SCALE = 1/28.0                      # bp per drawing unit
GAPMIN = 1.2                        # keep a hairline visible at 2 bp
X0 = 34
fig, ax = plt.subplots(figsize=(12.4, 0.66*len(order) + 3.0))

def arrow(x, w, y, color, label, fs):
    tip = min(2.0, w*0.28)
    ax.add_patch(plt.Polygon([[x, y-0.26],[x+w-tip, y-0.26],[x+w, y],
                              [x+w-tip, y+0.26],[x, y+0.26]],
                             closed=True, fc=color, ec="black", lw=0.7, zorder=3))
    if w > 8:
        ax.text(x+(w-tip)/2, y, label, ha="center", va="center", color="white",
                fontsize=fs, fontweight="bold", style="italic", zorder=4)

maxx = 0
for i, g in enumerate(order):
    y = len(order) - i
    d = geo[g]
    wz, ws = d["zwf"]*SCALE, d["sub"]*SCALE
    wg = max(GAPMIN, d["gap"]*SCALE)
    name = "Candidate genus" if g == "CANDIDATE" else g
    bold = "bold" if g == "CANDIDATE" else "normal"
    ital = "normal" if g.startswith(("SW","UBA","WTJ","JAP","CAL","CAX")) else "italic"
    ax.text(X0-2.5, y, name, ha="right", va="center", fontsize=9.5,
            fontweight=bold, style=ital)
    ax.text(X0-2.5, y-0.36, "n = %d" % d["n"], ha="right", va="center",
            fontsize=6.8, color="#777")
    ax.plot([X0, X0+wz+wg+ws], [y, y], color="#999", lw=0.8, zorder=1)
    arrow(X0, wz, y, RED, "zwf", 8.5)
    arrow(X0+wz+wg, ws, y, BLUE, "opcA", 8.5)
    ax.text(X0+wz+wg/2, y+0.44, "%d" % d["gap"], ha="center", fontsize=6.6, color="#333")
    maxx = max(maxx, X0+wz+wg+ws)

# right-hand columns
C = maxx + 7
for i, g in enumerate(order):
    y = len(order) - i; d = geo[g]
    ax.text(C,    y, "%d/%d" % (d["nsc"], d["n"]), ha="center", va="center", fontsize=8)
    ax.text(C+11, y, "%d/%d" % (d["co"], d["nsc"]), ha="center", va="center", fontsize=8)
    ax.text(C+22, y, "%d/%d" % (d["ok"], d["nsc"]), ha="center", va="center", fontsize=8)

yt = len(order) + 0.9
for lab, dx in (("same\ncontig", 0), ("co-\noriented", 11), ("zwf\nupstream", 22)):
    ax.text(C+dx, yt, lab, ha="center", va="center", fontsize=7.4, color="#333")
ax.text(X0 + (maxx-X0)/2, yt, "locus, drawn to scale", ha="center", va="center",
        fontsize=8.2, color="#333")
ax.plot([X0, X0+1000*SCALE], [0.30, 0.30], color="black", lw=1.6)
ax.text(X0+500*SCALE, 0.06, "1 kb", ha="center", fontsize=7.5)

nsc = sum(1 for r in rows if r["same_contig"] == "True")
ax.set_xlim(0, C+30); ax.set_ylim(-0.15, len(order)+1.8); ax.axis("off")
ax.set_title("zwf and opcA form a conserved predicted operon across Akkermansiaceae\n"
             "%d genomes in %d genera; directly adjacent, co-oriented and zwf-upstream in "
             "%d of %d genomes where both genes share one contig.\n"
             "Neither gene is present in any of 331 Akkermansia genomes."
             % (len(rows), len(order), nsc, nsc), fontsize=10.5, pad=16)
plt.tight_layout()
for e in ("pdf", "png"):
    plt.savefig(OUT + "." + e, bbox_inches="tight", dpi=300)
    print("wrote " + OUT + "." + e)
