#!/usr/bin/env python3
# Operon adjacency and intergenic gap measured
# Source: ch3-chitin-evolution/scripts/measure_ppp_operon2.py
# Output: stdout
import os, re
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
        s=x[1]; genome=s.split("|")[0].split("_CDS_")[0]
        if genome not in best: best[genome]=s
    return best
zwf=subj(zwf_q); g6pd=subj(g6pd_q)

# subject: GENOME|CONTIG_geneidx  where CONTIG may itself contain underscores
# e.g. EHM058709|EHA04264_bin.23^_185_5  -> contig 'EHA04264_bin.23^_185', gene 5
def contig_gene(s):
    tail=s.split("|")[1]
    m=re.match(r"(.+)_(\d+)$", tail)
    return (m.group(1), int(m.group(2))) if m else (None,None)

# GFF: seqid col0 = contig, ID=<contigidx>_<geneidx>. Build (contig, geneidx)->(start,end)
def gff_map(genome):
    fp=f"{GFF}/{genome}.gff"; d={}
    if not os.path.isfile(fp): return d
    for l in open(fp,errors="ignore"):
        if "\tCDS\t" not in l: continue
        f=l.split("\t")
        contig=f[0]; start=int(f[3]); end=int(f[4])
        m=re.search(r"ID=\d+_(\d+)",f[8])
        if m: d[(contig,int(m.group(1)))]=(start,end)
    return d

gaps=[]; same=0; both=0
for genome in set(zwf)&set(g6pd):
    both+=1
    zc,zi=contig_gene(zwf[genome]); gc,gi=contig_gene(g6pd[genome])
    gm=gff_map(genome)
    z=gm.get((zc,zi)); g=gm.get((gc,gi))
    if z and g and zc==gc:
        same+=1
        gap=max(g[0]-z[1], z[0]-g[1])   # bp between them
        gaps.append(abs(gap))
    elif zc==gc and zc is not None:
        same+=1  # same contig even if coord lookup failed

print(f"genomes with both genes: {both}")
print(f"same contig (adjacent operon): {same}/{both}")
if gaps:
    gaps.sort(); n=len(gaps); med=gaps[n//2] if n%2 else (gaps[n//2-1]+gaps[n//2])//2
    print(f"intergenic gap: median {med} bp  range {min(gaps)}-{max(gaps)}  n={n}")
    print(f"  under 100 bp: {sum(1 for x in gaps if x<100)}/{n}")
else:
    print("coord lookup still failing; sample:")
    g0=list(set(zwf)&set(g6pd))[0]
    print(f"  {g0}: zwf={zwf[g0]} -> {contig_gene(zwf[g0])}")
    print(f"          g6pd={g6pd[g0]} -> {contig_gene(g6pd[g0])}")
    gm=gff_map(g0)
    print(f"  gff keys sample: {list(gm.keys())[:3]}")
