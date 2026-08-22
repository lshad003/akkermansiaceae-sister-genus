#!/usr/bin/env python3
# Enzyme family panels
# Source: ch3-chitin-evolution/scripts/fig_cazy_panels.py
# Output: results/figures/Fig_cazy_panels.png
import csv, os, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
M=f"{BASE}/results/novel_akk_tree/cazy_heatmap_matrix.tsv"
PNG=f"{BASE}/results/figures/Fig_cazy_panels.png"

rows=list(csv.reader(open(M),delimiter="\t")); hdr=rows[0]; data=rows[1:]
groups=hdr[1:]
V={r[0]:[float(x) for x in r[1:]] for r in data}
NOV,AKKa,AKKall = 0,1,2
free=list(range(3,len(groups)))
def fmax(f): return max(V[f][i] for i in free)

cats=collections.OrderedDict()
cats["A. Ancestral, lost in Akkermansia\n(novel high, Akk low, free-living present)"]=[
    f for f in V if V[f][NOV]>=50 and V[f][AKKall]<=15 and fmax(f)>=25]
cats["B. Lost in BOTH gut lineages\n(novel + Akk absent, free-living present)"]=[
    f for f in V if V[f][NOV]<=5 and V[f][AKKall]<=5 and fmax(f)>=50]
cats["C. Gained in Akkermansia\n(Akk high, novel low)"]=[
    f for f in V if V[f][AKKall]>=40 and V[f][NOV]<=15]
cats["D. Novel-lineage specific\n(novel high, Akk and free-living low)"]=[
    f for f in V if V[f][NOV]>=50 and V[f][AKKall]<=15 and fmax(f)<25]

for k in cats: cats[k].sort(key=lambda f:-(V[f][NOV]-V[f][AKKall]))
cats={k:v[:14] for k,v in cats.items() if v}
for k,v in cats.items(): print(f"{k.splitlines()[0]:55s} {len(v)}")

n=len(cats); hs=[len(v) for v in cats.values()]
fig,axes=plt.subplots(n,1,figsize=(1.05*len(groups)+3, sum(0.34*h+1.5 for h in hs)),
                      gridspec_kw={"height_ratios":hs})
if n==1: axes=[axes]
for ax,(title,fams) in zip(axes,cats.items()):
    A=np.array([V[f] for f in fams])
    im=ax.imshow(A,aspect="auto",cmap="RdYlBu_r",vmin=0,vmax=100)
    ax.set_yticks(range(len(fams))); ax.set_yticklabels(fams,fontsize=8)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g.replace(" (","\n(") for g in groups],fontsize=7)
    for i in range(len(fams)):
        for j in range(len(groups)):
            v=A[i,j]
            ax.text(j,i,"%.0f"%v,ha="center",va="center",fontsize=6,
                    color="white" if (v>70 or v<12) else "black")
    ax.axvline(2.5,color="k",lw=1.6)
    ax.set_title(title,fontsize=9,loc="left")
fig.suptitle("CAZyme repertoire: novel genus vs Akkermansia vs free-living relatives",fontsize=11)
fig.tight_layout(rect=[0,0,1,0.985]); fig.savefig(PNG,dpi=200)
print("wrote",PNG)
