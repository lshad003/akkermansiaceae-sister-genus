#!/usr/bin/env python3
# Figure 3
# Source: ch3-chitin-evolution/scripts/fig3_mucin_5group.py
# Output: results/figures/Fig3_mucin_5group.png
import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SRC=f"{BASE}/results/novel_akk_tree/novel_genus_function.tsv"
OUT=f"{BASE}/results/figures/Fig3_mucin_5group.png"
os.makedirs(os.path.dirname(OUT),exist_ok=True)

rows={r["family"]:r for r in csv.DictReader(open(SRC),delimiter="\t")}
mucin=[f for f in rows if rows[f].get("in_mucin_panel")=="1"]
mucin.sort(key=lambda f: float(rows[f]["pct_novel"]), reverse=True)

# EDIT these keys if the header check shows different names:
cols=[("pct_novel","Ca. Novel genus\namphibian (105)"),
      ("pct_akk_amph","Akkermansia\namphibian (94)"),
      ("pct_akk_podarcis","Akkermansia\nPodarcis (137)"),
      ("pct_akk_mammal","Akkermansia\nmammal (71)"),
      ("pct_akk_gtdb","Akkermansia\nGTDB (42)")]

M=np.array([[float(rows[f][c]) for c,_ in cols] for f in mucin])

fig,ax=plt.subplots(figsize=(9,8))
im=ax.imshow(M,cmap="RdYlBu_r",vmin=0,vmax=100,aspect="auto")
ax.set_xticks(range(len(cols))); ax.set_xticklabels([lbl for _,lbl in cols],fontsize=9)
ax.set_yticks(range(len(mucin))); ax.set_yticklabels(mucin,fontsize=10)
for i in range(len(mucin)):
    for j in range(len(cols)):
        v=M[i,j]
        ax.text(j,i,f"{v:.0f}",ha="center",va="center",
                color="white" if (v>75 or v<25) else "black",fontsize=9)
cb=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.03); cb.set_label("% genomes")
ax.set_title("Mucin-degradation machinery is conserved across all gut-associated Akkermansiaceae",
             fontsize=12,pad=12)
fig.tight_layout()
fig.savefig(OUT,dpi=200,bbox_inches="tight")
print("wrote",OUT)
print("families:",mucin)
print("median per group:")
for j,(c,lbl) in enumerate(cols):
    print(f"  {lbl.split(chr(10))[0]:<18} median {np.median(M[:,j]):.1f}")
