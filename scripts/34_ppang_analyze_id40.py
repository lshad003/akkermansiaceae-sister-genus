#!/usr/bin/env python3
# Partitions summarized at the corrected identity threshold
# Source: ch3-chitin-evolution/scripts/ppang_analyze_id40.py
# Output: results/pangenome/novel_diagnostic_families_id40.txt
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"; O=f"{P}/ppang_out_id40"
grp={}; animal={}
for r in csv.DictReader(open(f"{P}/genome_groups.tsv"),delimiter="\t"):
    grp[r["genome"]]=r["group"]; animal[r["genome"]]=r["host_animal"]

print("=== partitions dir ===")
for f in sorted(os.listdir(f"{O}/partitions")):
    n=sum(1 for _ in open(f"{O}/partitions/{f}"))
    print(f"  {f:<28}{n}")

rt=f"{O}/gene_presence_absence.Rtab"
with open(rt) as fh:
    hdr=fh.readline().rstrip("\n").split("\t")
    genomes=hdr[1:]
    fam={}
    for line in fh:
        f=line.rstrip("\n").split("\t")
        fam[f[0]]={genomes[i] for i,v in enumerate(f[1:]) if v=="1"}
print(f"\nRtab: {len(fam)} families x {len(genomes)} genomes")
bad=[g for g in genomes if g not in grp]
if bad: print("  NOT IN GROUPS:",bad[:5])

N={g for g in grp if grp[g]=="NOVEL"}; A={g for g in grp if grp[g]=="AKK_AMPH"}
print(f"NOVEL {len(N)} / {len(set(animal[g] for g in N))} animals")
print(f"AKK_AMPH {len(A)} / {len(set(animal[g] for g in A))} animals")

part={}
for f in os.listdir(f"{O}/partitions"):
    lab=f.replace(".txt","")
    for line in open(f"{O}/partitions/{f}"):
        part.setdefault(line.strip(),set()).add(lab)

def summary(gs,lab):
    n=len(gs); occ=collections.Counter()
    tot=pers=0
    for k,v in fam.items():
        c=len(v&gs)
        if c: tot+=1; occ[c]+=1
    print(f"\n=== {lab} (n={n}) ===")
    print(f"  families present in >=1: {tot}")
    for lo in (0.95,0.90,0.50,0.15):
        c=sum(1 for k,v in fam.items() if len(v&gs)>=lo*n)
        print(f"  present in >={int(lo*100)}%: {c}")
summary(N,"NOVEL"); summary(A,"AKK_AMPH")

print("\n=== PPanGGOLiN partition x group ===")
for lab in sorted({x for v in part.values() for x in v}):
    fams=[k for k,v in part.items() if lab in v]
    nn=sum(1 for k in fams if k in fam and (fam[k]&N))
    aa=sum(1 for k in fams if k in fam and (fam[k]&A))
    both=sum(1 for k in fams if k in fam and (fam[k]&N) and (fam[k]&A))
    print(f"  {lab:<22} total {len(fams):>6}  inNOVEL {nn:>6}  inAKK {aa:>6}  shared {both:>6}")

print("\n=== EXCLUSIVE ===")
non=[k for k,v in fam.items() if (v&N) and not (v&A)]
aon=[k for k,v in fam.items() if (v&A) and not (v&N)]
print(f"  NOVEL-exclusive: {len(non)}   AKK-exclusive: {len(aon)}   shared: {sum(1 for k,v in fam.items() if (v&N) and (v&A))}")
nd=[k for k in non if len(fam[k]&N)>=0.90*len(N)]
ad=[k for k in aon if len(fam[k]&A)>=0.90*len(A)]
print(f"  NOVEL-exclusive AND in >=90% of NOVEL: {len(nd)}   <- genus-diagnostic")
print(f"  AKK-exclusive AND in >=90% of AKK:     {len(ad)}")
with open(f"{P}/novel_diagnostic_families_id40.txt","w") as fh:
    for k in nd: fh.write(k+"\n")
print(f"\nwrote {P}/novel_diagnostic_families_id40.txt")
