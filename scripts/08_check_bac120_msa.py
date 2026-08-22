#!/usr/bin/env python3
# Reference-containing bac120 alignment confirmed as the placement input
# Source: ch3-chitin-evolution/scripts/check_bac120_msa.py
# Output: stdout
import gzip, os, csv, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
A=f"{BASE}/results/gtdbtk_amphibia/align/align"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
cen={r["accession"]:r for r in csv.DictReader(open(CEN),delimiter="\t")}
def g(r): return (r["genus"] or "").strip()
novel={a for a,r in cen.items() if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian"
       and r["annotated"]=="1" and g(r) in ("","unknown","NO_GENUS")}
named={a for a,r in cen.items() if r["family"]=="Akkermansiaceae" and g(r)=="Akkermansia"}
print("target novel:",len(novel)," named Akkermansia in census:",len(named))
for fn in ("gtdbtk.bac120.msa.fasta.gz","gtdbtk.bac120.user_msa.fasta.gz"):
    p=f"{A}/{fn}"
    if not os.path.exists(p): print(f"\n{fn}: MISSING"); continue
    ids=[]; L=None
    with gzip.open(p,"rt") as fh:
        cur=None; n=0
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0]); n+=1
                if n==2 and L is None: L=len(cur or "")
                cur=""
            else: cur=(cur or "")+line.strip()
    print(f"\n{fn}: {len(ids)} seqs, aln len ~{L}")
    print("  first 3 ids:",ids[:3])
    idset=set(ids)
    def strip(x):
        return x[3:] if x[:3] in ("GB_","RS_") else x
    stripped={strip(i) for i in ids}
    hitn=sum(1 for a in novel if a in idset or strip(a) in stripped)
    hitk=sum(1 for a in named if a in idset or strip(a) in stripped)
    print(f"  contains novel: {hitn}/{len(novel)}")
    print(f"  contains named Akkermansia: {hitk}/{len(named)}")
    pref=collections.Counter(i[:3] for i in ids)
    print("  id prefixes:",dict(pref.most_common(5)))
