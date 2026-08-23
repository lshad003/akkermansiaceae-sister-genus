#!/usr/bin/env python3
# Operon co-orientation checked
# Source: ch3-chitin-evolution/scripts/check_ppp_strand.py
# Output: stdout
import os, re, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"; GFF=f"{P}/gff199"
zwf_q="UHM979.41089_R.bin.103_CDS_0654"; g6pd_q="EHM058980_CDS_1432"

def subj(qid,pid=30.0):
    best={}
    for l in open(f"{P}/ppp_all199.tsv",errors="ignore"):
        x=l.rstrip("\n").split("\t")
        if len(x)<3 or x[0]!=qid: continue
        try: p=float(x[2])
        except: continue
        if p<pid: continue
        s=x[1]; g=s.split("|")[0].split("_CDS_")[0]
        if g not in best: best[g]=s
    return best
zwf=subj(zwf_q); g6pd=subj(g6pd_q)

def contig_gene(s):
    tail=s.split("|")[1]; m=re.match(r"(.+)_(\d+)$",tail)
    return (m.group(1),int(m.group(2))) if m else (None,None)

def gff_strand(genome):
    fp=f"{GFF}/{genome}.gff"; d={}
    if not os.path.isfile(fp): return d
    for l in open(fp,errors="ignore"):
        if "\tCDS\t" not in l: continue
        f=l.split("\t"); contig=f[0]; strand=f[6]
        m=re.search(r"ID=\d+_(\d+)",f[8])
        if m: d[(contig,int(m.group(1)))]=strand
    return d

same_strand=0; opp=0; total=0
for genome in set(zwf)&set(g6pd):
    zc,zi=contig_gene(zwf[genome]); gc,gi=contig_gene(g6pd[genome])
    if zc!=gc: continue
    st=gff_strand(genome)
    zs=st.get((zc,zi)); gs=st.get((gc,gi))
    if zs and gs:
        total+=1
        if zs==gs: same_strand+=1
        else: opp+=1

print(f"same-contig pairs with strand data: {total}")
print(f"  same strand (co-oriented): {same_strand}")
print(f"  opposite strand: {opp}")
if total:
    print(f"  -> {'CO-ORIENTED, arrows OK' if same_strand/total>0.9 else 'MIXED, use rectangles'}")
