#!/usr/bin/env python3
# Figure 2A
# Source: ch3-chitin-evolution/scripts/fig_ppp_final.py
# Output: results/figures/Fig2_PPP_loss_verified.png
import os, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"
OUT=f"{BASE}/results/figures/Fig2_PPP_loss_verified.png"

trio={"UHM979.41089_R.bin.103_CDS_0654":"zwf",
      "UHM1210.23070_R.bin.101_CDS_0293":"gnd",
      "EHM058980_CDS_1432":"g6pd_sub"}
# exact denominators from outgroups_all.faa (verified this session)
FREE_DENOM={"Haloferula":15,"Luteolibacter":115,"Rubritalea":11,"Roseibacillus":7}  # Rosei = 2 + Rosei_B 5

def parse(path,pid=30.0):
    h=collections.defaultdict(set)
    for l in open(path,errors="ignore"):
        x=l.rstrip("\n").split("\t")
        if len(x)<3: continue
        try: p=float(x[2])
        except: continue
        if p<pid: continue
        s=x[1]; genus=s.split("::")[0] if "::" in s else "?"
        acc=s.split("::")[-1].split("|")[0]
        if genus=="Roseibacillus_B": genus="Roseibacillus"
        h[x[0]].add((genus,acc))
    return h

og=parse(f"{P}/ppp_outgroups.tsv")
cand=parse(f"{P}/n22_all199.tsv")
novel=set()
for l in open(f"{P}/genome_groups.tsv"):
    x=l.rstrip("\n").split("\t")
    if len(x)>1 and x[1]=="NOVEL": novel.add(x[0])
def cand_n(cds):
    return len({a for (g,a) in cand.get(cds,set()) if a in novel})

freegenera=["Haloferula","Luteolibacter","Rubritalea","Roseibacillus"]
cols=["Ca. Novel\ngenus (105)","Akkermansia\n(0/326)"]+[f"{g}\n({FREE_DENOM[g]})" for g in freegenera]
M=[]; ann=[]
for cds,gene in trio.items():
    row=[]; a=[]
    cn=cand_n(cds); row.append(100*cn/105); a.append(f"{cn}/105")
    row.append(0.0); a.append("0/326")
    for g in freegenera:
        n=len({acc for (gg,acc) in og.get(cds,set()) if gg==g})
        d=FREE_DENOM[g]; row.append(100*n/d); a.append(f"{n}/{d}")
    M.append(row); ann.append(a)
M=np.array(M)

fig,ax=plt.subplots(figsize=(12,4.3))
im=ax.imshow(M,cmap="Reds",vmin=0,vmax=100,aspect="auto")
ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols,fontsize=9)
ax.set_yticks(range(3)); ax.set_yticklabels(["zwf","gnd","G6PD subunit"],fontsize=12,style="italic")
for i in range(3):
    for j in range(len(cols)):
        v=M[i,j]
        ax.text(j,i,ann[i][j],ha="center",va="center",
                color="white" if v>55 else "black",
                fontweight="bold" if j==1 else "normal",fontsize=10)
ax.add_patch(plt.Rectangle((0.5,-0.5),1,3,fill=False,edgecolor="#2166ac",lw=2.5))
cb=fig.colorbar(im,ax=ax,fraction=0.02,pad=0.02); cb.set_label("% genomes")
ax.set_title("Oxidative pentose phosphate pathway: present in the candidate lineage and every\nfree-living genus, absent from all 326 Akkermansia",fontsize=12,pad=10)
fig.tight_layout()
fig.savefig(OUT,dpi=200,bbox_inches="tight")
print("wrote",OUT)
for cds,gene in trio.items():
    freestr=" ".join(f"{g}:{len({a for (gg,a) in og.get(cds,set()) if gg==g})}/{FREE_DENOM[g]}" for g in freegenera)
    print(f"  {gene}: cand {cand_n(cds)}/105  akk 0/326  {freestr}")
