#!/usr/bin/env python3
# Figure: oxidative PPP across the family, from the single unified search.
# SUPERSEDES Figure_ppp_loss.pdf, which showed the dead 0/326 denominator.
import os
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W   = R + "/results/ppp_unified"
CEN = R + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
NOV = R + "/results/novel_akk_tree/novel_size_gc.tsv"
OUT = R + "/results/figures/Figure_ppp_loss.pdf"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

QN = {"UHM979.41089_R.bin.103_CDS_0654": "zwf",
      "UHM1210.23070_R.bin.101_CDS_0293": "gnd",
      "EHM058980_CDS_1432": "opcA"}

def load(p):
    d = defaultdict(set)
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            d[f[1].split("|")[0]].add(QN.get(f[0], f[0]))
    return d
ppp, ctl = load(W + "/ppp_trio.tsv"), load(W + "/control_gh20.tsv")

indb = set()
for line in open(W + "/manifest.tsv"):
    indb.add(line.split("\t")[0].strip())

with open(NOV) as f:
    f.readline(); cand = set(l.split("\t")[0].strip() for l in f if l.strip())

fam, gen = {}, {}
with open(CEN) as f:
    h = f.readline().rstrip("\n").split("\t")
    ifm, ign = h.index("family"), h.index("genus")
    s = lambda v: v[3:] if v[:3] in ("f__", "g__") else v
    for line in f:
        p = line.rstrip("\n").split("\t")
        if ifm < len(p): fam[p[0].strip()] = s(p[ifm].strip())
        if ign < len(p): gen[p[0].strip()] = s(p[ign].strip())

op = {}
with open(R + "/results/ppp_unified/operon_allgenera.tsv") as f:
    oh = f.readline().rstrip("\n").split("\t")
    for line in f:
        d = dict(zip(oh, line.rstrip("\n").split("\t")))
        if d["same_contig"] == "True" and d["order"] == "zwf-sub":
            op[d["genome"]] = True

cols = defaultdict(list)
for g in indb:
    if g in cand: cols["Candidate\ngenus"].append(g)
    elif gen.get(g) == "Akkermansia": cols["Akkermansia"].append(g)
    elif fam.get(g) == "Akkermansiaceae":
        gg = gen.get(g) or "unassigned"
        if gg: cols[gg].append(g)

keep = ["Candidate\ngenus", "Akkermansia"] + \
       sorted([k for k in cols if k not in ("Candidate\ngenus", "Akkermansia")],
              key=lambda k: -len(cols[k]))[:6]
rows = ["zwf", "gnd", "opcA", "zwf-opcA adjacent", "GH20 control"]

M = np.zeros((len(rows), len(keep))); lab = []
for j, c in enumerate(keep):
    gl = cols[c]; n = len(gl); col = []
    for i, r in enumerate(rows):
        if r == "GH20 control":
            k = sum(1 for g in gl if ctl.get(g))
        elif r == "zwf-opcA adjacent":
            k = sum(1 for g in gl if op.get(g))
        elif r == "opcA":
            k = sum(1 for g in gl if "opcA" in ppp.get(g, set()))
        else:
            k = sum(1 for g in gl if r in ppp.get(g, set()))
        M[i, j] = 100.0 * k / n if n else 0
        col.append("%d/%d" % (k, n))
    lab.append(col)

fig, ax = plt.subplots(figsize=(11, 3.6))
im = ax.imshow(M, cmap="Reds", vmin=0, vmax=100, aspect="auto")
ax.set_xticks(range(len(keep)))
ax.set_xticklabels(["%s\n(%d)" % (c.replace("\n", " "), len(cols[c])) for c in keep], fontsize=9)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels(rows, fontsize=10, style="italic")
for i in range(len(rows)):
    for j in range(len(keep)):
        ax.text(j, i, lab[j][i], ha="center", va="center", fontsize=8.5,
                color="white" if M[i, j] > 55 else "black",
                fontweight="bold" if keep[j] == "Akkermansia" else "normal")
ja = keep.index("Akkermansia")
ax.add_patch(plt.Rectangle((ja-0.5, -0.5), 1, len(rows), fill=False, ec="#1f4e79", lw=2.2))
ax.axhline(2.5, color="white", lw=3)
ax.axhline(3.5, color="white", lw=3)
ax.set_title("Oxidative pentose phosphate pathway across Akkermansiaceae\n"
             "single DIAMOND search; zwf-opcA adjacency and GH20 as positive control", fontsize=11)
plt.colorbar(im, ax=ax, label="% genomes", fraction=0.025, pad=0.015)
plt.tight_layout()
for ext in ("pdf", "png"):
    f = OUT[:-4] + "." + ext
    plt.savefig(f, bbox_inches="tight", dpi=300)
    print("wrote " + f)
