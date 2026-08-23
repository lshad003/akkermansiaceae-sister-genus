#!/usr/bin/env python3
# Orthogroups polarized against free-living outgroups
# Source: ch3-chitin-evolution/scripts/orthogroup_polarity.py
# Output: results/pangenome/orthogroup_polarity.tsv
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
R=f"{BASE}/results/pangenome/orthofinder_out/run1/Results_novelakk"
GC=f"{R}/Orthogroups/Orthogroups.GeneCount.tsv"
OUT=f"{BASE}/results/pangenome/orthogroup_polarity.tsv"

r=csv.reader(open(GC),delimiter="\t"); hdr=next(r)
sp=hdr[1:-1] if hdr[-1].lower().startswith("total") else hdr[1:]
def grp(s):
    if s.startswith("NOVEL"): return "NOVEL"
    if s.startswith("AKK"): return "AKK"
    return "FREE"
def genus(s): return s.split("_rep")[0]
cols={g:[i for i,s in enumerate(sp) if grp(s)==g] for g in ("NOVEL","AKK","FREE")}
fgen=sorted(set(genus(s) for s in sp if grp(s)=="FREE"))
gidx={g:[i for i,s in enumerate(sp) if grp(s)=="FREE" and genus(s)==g] for g in fgen}
print(f"species: NOVEL {len(cols['NOVEL'])} | AKK {len(cols['AKK'])} | FREE {len(cols['FREE'])} ({len(fgen)} genera)")

rows=[]
for row in r:
    og=row[0]
    vals=[int(x) for x in row[1:1+len(sp)]]
    pres=[1 if v>0 else 0 for v in vals]
    pn=100.0*sum(pres[i] for i in cols["NOVEL"])/len(cols["NOVEL"])
    pa=100.0*sum(pres[i] for i in cols["AKK"])/len(cols["AKK"])
    ngen=sum(1 for g in fgen if any(pres[i] for i in gidx[g]))
    rows.append((og,pn,pa,ngen))
print(f"orthogroups: {len(rows)}")

A=[x for x in rows if x[1]>=70 and x[2]<=10 and x[3]>=3]
B=[x for x in rows if x[1]<=5  and x[2]<=5  and x[3]>=5]
C=[x for x in rows if x[2]>=70 and x[1]<=10]
D=[x for x in rows if x[1]>=70 and x[2]<=10 and x[3]==0]
SH=[x for x in rows if x[1]>=70 and x[2]>=70]

print("\n=== GENOME-WIDE ORTHOGROUP POLARITY ===")
print(f"A  ancestral, lost on Akkermansia branch  (novel>=70, Akk<=10, >=3 free genera): {len(A)}")
print(f"B  lost at the shared gut ancestor        (both gut <=5, >=5 free genera):       {len(B)}")
print(f"C  Akkermansia-enriched, novel low        (Akk>=70, novel<=10):                  {len(C)}")
print(f"D  novel-lineage specific                 (novel>=70, Akk<=10, no free-living):  {len(D)}")
print(f"   shared by both gut lineages            (both >=70):                           {len(SH)}")

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["orthogroup","pct_novel","pct_akk","n_free_genera","category"])
    for lab,S in (("A_lost_in_Akkermansia",A),("B_lost_at_gut_ancestor",B),
                  ("C_Akkermansia_enriched",C),("D_novel_specific",D),("shared_gut",SH)):
        for og,pn,pa,ng in S: w.writerow([og,f"{pn:.1f}",f"{pa:.1f}",ng,lab])
print(f"\nwrote {OUT}")
