#!/usr/bin/env python3
# Enzyme family heatmap
# Source: ch3-chitin-evolution/scripts/fig_cazy_heatmap.py
# Output: results/figures/Fig_cazy_heatmap.png
import csv, os, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUTD=f"{BASE}/results/figures"; os.makedirs(OUTD,exist_ok=True)
PNG=f"{OUTD}/Fig_cazy_heatmap.png"; TSV=f"{BASE}/results/novel_akk_tree/cazy_heatmap_matrix.tsv"
DIRS=[f"{BASE}/results/dbcan_verru", f"{BASE}/results/dbcan_bact_refs",
      f"{BASE}/results/dbcan_ehi_amphibian", f"{BASE}/results/dbcan_ehi_nonamph",
      f"{BASE}/results/dbcan_endo_allphyla", f"{BASE}/results/dbcan_flavo_refs",
      f"{BASE}/results/dbcan_scaffold", RUM]
E=1e-10; COV=0.30
FREE=["Luteolibacter","Haloferula","Rubritalea","Roseibacillus_B","SW10","Oceaniferula","UBA956"]

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t") if r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()
G=collections.OrderedDict()
G["Novel genus\n(n=105)"]=[r for r in rows if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian" and gen(r) in ("","unknown","NO_GENUS")]
G["Akkermansia\namphibian (n=94)"]=[r for r in rows if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian" and gen(r)=="Akkermansia"]
G["Akkermansia\nall (n=346)"]=[r for r in rows if r["family"]=="Akkermansiaceae" and gen(r)=="Akkermansia"]
for g in FREE:
    v=[r for r in rows if gen(r)==g]
    if len(v)>=10: G[f"{g}\n(n={len(v)})"]=v

idx={}
for d in DIRS:
    if not os.path.isdir(d): continue
    for fn in os.listdir(d):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            idx.setdefault(fn[:-4], os.path.join(d,fn))
def res(a):
    if a in idx: return idx[a]
    if a[:3] in ("GB_","RS_") and a[3:] in idx: return idx[a[3:]]
    return None
def fam_of(h): return h[:-4].split("_")[0] if h.endswith(".hmm") else h.split("_")[0]

FG={}
for k,v in G.items():
    for r in v:
        a=r["accession"]
        if a in FG: continue
        p=res(a)
        if not p: continue
        s=set()
        for line in open(p):
            f=line.rstrip("\n").split("\t")
            if len(f)<10: continue
            try: ev=float(f[4]); cv=float(f[9])
            except: continue
            if cv<COV or ev>E: continue
            s.add(fam_of(f[0]))
        FG[a]=s

groups=list(G)
def prev(k,f):
    v=[r["accession"] for r in G[k] if r["accession"] in FG]
    return 100.0*sum(1 for a in v if f in FG[a])/len(v) if v else 0.0

allf=sorted(set().union(*FG.values()))
M={f:[prev(k,f) for k in groups] for f in allf}
# keep informative families: >=40% somewhere AND >=30 point spread
sel=[f for f in allf if max(M[f])>=40 and (max(M[f])-min(M[f]))>=30]
# order by novel-minus-Akkermansia_all
sel.sort(key=lambda f: -(M[f][0]-M[f][2]))
print("families plotted:",len(sel))

with open(TSV,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["family"]+[g.replace("\n"," ") for g in groups])
    for f in sel: w.writerow([f]+["%.1f"%x for x in M[f]])

A=np.array([M[f] for f in sel])
fig,ax=plt.subplots(figsize=(1.15*len(groups)+3, 0.30*len(sel)+2.2))
im=ax.imshow(A, aspect="auto", cmap="RdYlBu_r", vmin=0, vmax=100)
ax.set_xticks(range(len(groups))); ax.set_xticklabels(groups, fontsize=8)
ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel, fontsize=7)
for i in range(len(sel)):
    for j in range(len(groups)):
        v=A[i,j]
        ax.text(j,i,"%.0f"%v,ha="center",va="center",fontsize=5.5,
                color="white" if (v>70 or v<12) else "black")
ax.axvline(2.5,color="k",lw=1.6)
ax.set_title("CAZyme family prevalence: novel genus, Akkermansia, and free-living relatives\n"
             "(dbCAN, E<=1e-10, coverage>=0.30; % of genomes carrying >=1 domain)",fontsize=10)
cb=fig.colorbar(im,ax=ax,shrink=0.5); cb.set_label("% genomes",fontsize=8)
fig.tight_layout(); fig.savefig(PNG,dpi=200)
print("wrote",PNG); print("wrote",TSV)
