#!/usr/bin/env python3
# Family prevalence by free-living genus
# Source: ch3-chitin-evolution/scripts/cazy_by_freeliving_genus.py
# Output: results/novel_akk_tree/cazy_by_genus.tsv
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUT=f"{BASE}/results/novel_akk_tree/cazy_by_genus.tsv"
DIRS=[f"{BASE}/results/dbcan_verru", f"{BASE}/results/dbcan_bact_refs",
      f"{BASE}/results/dbcan_ehi_amphibian", f"{BASE}/results/dbcan_ehi_nonamph",
      f"{BASE}/results/dbcan_endo_allphyla", f"{BASE}/results/dbcan_flavo_refs",
      f"{BASE}/results/dbcan_scaffold", RUM]
E=1e-10; COV=0.30
FOCUS=["GH92","GH139","GH120","PL33","GH154","GH38","CBM91","GT92","GH73","GH97",
       "GH105","GH74","GH75","GH18","GH163","GH141","GH151"]
FREE=["Luteolibacter","Haloferula","Rubritalea","Roseibacillus","Roseibacillus_B",
      "SW10","Oceaniferula","UBA956"]

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t") if r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()
G=collections.OrderedDict()
G["NOVEL"]=[r for r in rows if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian" and gen(r) in ("","unknown","NO_GENUS")]
G["AKK_all"]=[r for r in rows if r["family"]=="Akkermansiaceae" and gen(r)=="Akkermansia"]
for g in FREE:
    v=[r for r in rows if gen(r)==g]
    if len(v)>=4: G[g]=v

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

groups=[k for k in G]
print("%-8s" % "family", "".join("%12s"%k[:11] for k in groups))
out=[]
for f in FOCUS:
    line=[]
    for k in groups:
        v=[r["accession"] for r in G[k] if r["accession"] in FG]
        p=100.0*sum(1 for a in v if f in FG[a])/len(v) if v else float("nan")
        line.append(p)
    print("%-8s"%f, "".join("%11.1f%%"%p for p in line))
    out.append([f]+["%.1f"%p for p in line])

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["family"]+groups); w.writerows(out)
print("\nn per group:", {k:len([r for r in G[k] if r['accession'] in FG]) for k in groups})
print("wrote",OUT)
