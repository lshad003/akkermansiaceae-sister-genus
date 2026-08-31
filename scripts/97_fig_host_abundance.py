#!/usr/bin/env python3
# Prevalence and abundance by host order figure
# Source: ch3-chitin-evolution/scripts/fig_host_abundance.py
# Output: results/figures/Figure_host_abundance.pdf
import csv, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
T = CH3 + "results/rurik_16s/per_sample_abundance.tsv"
OUT = CH3 + "results/figures/Figure_host_abundance"

ORD = ["Caudata", "Anura", "Squamata", "Testudines"]
LAB = ["Salamanders", "Frogs", "Lizards\nand snakes", "Turtles"]
COL = ["#2E6DA4", "#4A90C4", "#C0562B", "#8A6520"]

rows = list(csv.DictReader(open(T), delimiter="\t"))
pct = collections.defaultdict(list)
for r in rows:
    pct[r["host_order"]].append(float(r["percent"]))

fig, ax = plt.subplots(1, 2, figsize=(9.4, 4.3), gridspec_kw={"width_ratios": [1, 1.35]})

# left: prevalence
prev, ns = [], []
for o in ORD:
    v = pct[o]
    prev.append(100.0*sum(1 for x in v if x > 0)/len(v))
    ns.append((sum(1 for x in v if x > 0), len(v)))
b = ax[0].bar(range(4), prev, color=COL, width=0.62, edgecolor="none")
for i, (p, (a, t)) in enumerate(zip(prev, ns)):
    ax[0].text(i, p+1.6, "%.0f%%" % p, ha="center", fontsize=10, color="#15181B")
    ax[0].text(i, -5.5, "%d/%d" % (a, t), ha="center", fontsize=8, color="#8B9198")
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(LAB, fontsize=9)
ax[0].set_ylabel("samples where the genus is detected (%)", fontsize=9)
ax[0].set_ylim(0, 78); ax[0].set_xlim(-0.6, 3.6)
ax[0].spines[["top","right"]].set_visible(False)
ax[0].tick_params(labelsize=8)
ax[0].set_title("a   Prevalence", loc="left", fontsize=11, fontweight="600", pad=10)

# right: abundance among positives
for i, o in enumerate(ORD):
    v = [x for x in pct[o] if x > 0]
    if not v: continue
    x = np.random.default_rng(7+i).normal(i, 0.075, len(v))
    ax[1].scatter(x, v, s=13, color=COL[i], alpha=0.5, linewidths=0)
    med = float(np.median(v))
    ax[1].plot([i-0.26, i+0.26], [med, med], color="#15181B", lw=1.8, zorder=3)
    ax[1].text(i+0.31, med, "%.2f%%" % med, va="center", fontsize=8, color="#15181B")
ax[1].axhline(1.0, color="#8B9198", lw=1, ls="--")
ax[1].text(3.62, 1.0, "1%", va="center", fontsize=8, color="#8B9198")
for i, o in enumerate(ORD):
    n = sum(1 for x in pct[o] if x >= 1.0)
    ax[1].text(i, 78, str(n), ha="center", fontsize=9.5, color="#15181B")
ax[1].text(-0.62, 78, "above 1%:", fontsize=8.5, color="#8B9198", ha="left")
ax[1].set_yscale("log")
ax[1].set_ylim(0.0018, 130)
ax[1].set_xlim(-0.6, 3.6)
ax[1].set_xticks(range(4)); ax[1].set_xticklabels(LAB, fontsize=9)
ax[1].set_ylabel("share of the community where detected (%, log)", fontsize=9)
ax[1].spines[["top","right"]].set_visible(False)
ax[1].tick_params(labelsize=8)
ax[1].set_title("b   Abundance where present", loc="left", fontsize=11, fontweight="600", pad=10)

fig.tight_layout()
fig.savefig(OUT+".pdf")
fig.savefig(OUT+".png", dpi=300)
print("WROTE:", OUT+".pdf")
print("")
for i, o in enumerate(ORD):
    v = [x for x in pct[o] if x > 0]
    print("%-12s positive %3d of %3d   median %.3f%%   above 1%%: %d   max %.2f%%" % (
        o, len(v), len(pct[o]), float(np.median(v)) if v else 0,
        sum(1 for x in pct[o] if x >= 1.0), max(pct[o])))
