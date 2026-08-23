#!/usr/bin/env python3
# Figure 2B
# Source: ch3-chitin-evolution/scripts/fig_ppp_operon_final.py
# Output: results/figures/Fig2B_ppp_operon.png
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
OUT=f"{BASE}/results/figures/Fig4B_ppp_operon.png"
os.makedirs(os.path.dirname(OUT),exist_ok=True)
RED="#c0392b"

fig,ax=plt.subplots(figsize=(11,2.7))
ax.set_xlim(0,100); ax.set_ylim(0,10); ax.axis("off")

def gene(x,w,label,fs=12):
    ax.add_patch(plt.Polygon([[x,3],[x+w-1.4,3],[x+w,4.5],[x+w-1.4,6],[x,6]],
                             closed=True,facecolor=RED,edgecolor="black",lw=1))
    ax.text(x+(w-1.4)/2,4.5,label,ha="center",va="center",color="white",
            fontsize=fs,fontweight="bold",style="italic")

gene(5,17,"zwf")
gene(24,30,"G6PD accessory\nsubunit",fs=10)
# gap marker
ax.annotate("",xy=(24,7),xytext=(22,7),arrowprops=dict(arrowstyle="<->",color="black",lw=1))
ax.text(23,7.9,"~20 bp",ha="center",fontsize=9)

# unlinked break
ax.plot([58,74],[4.5,4.5],ls=":",color="black",lw=1.6)
ax.text(66,5.6,"unlinked",ha="center",fontsize=10,style="italic")
gene(76,17,"gnd")

ax.text(1,8.7,"B",fontsize=15,fontweight="bold")
ax.text(50,9.5,"Genomic organization of the oxidative PPP",ha="center",fontsize=12)
ax.text(50,1.3,"Conserved co-oriented adjacency in 91/95 genomes containing both genes; "
               "median intergenic distance 20 bp (range 12 to 60 bp).",
        ha="center",fontsize=9,color="#333333")

fig.savefig(OUT,dpi=200,bbox_inches="tight")
print("wrote",OUT)
