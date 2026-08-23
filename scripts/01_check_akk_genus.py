#!/usr/bin/env python3
# Genus composition of the Akkermansiaceae reference set
# Source: ch3-chitin-evolution/scripts/check_akk_genus.py
# Output: stdout
import csv, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"]
print("annotated Akkermansiaceae:",len(rows))
for hc in ("amphibian","reptile","mammal","unknown"):
    sub=[r for r in rows if r["host_class"]==hc]
    c=collections.Counter(r["genus"] or "NO_GENUS" for r in sub)
    print(f"\n{hc}  n={len(sub)}")
    for g,n in c.most_common(10): print(f"    {g:<30} {n}")
print("\n=== genus x host_class overlap ===")
gh=collections.defaultdict(collections.Counter)
for r in rows: gh[r["genus"] or "NO_GENUS"][r["host_class"]]+=1
for g,c in sorted(gh.items(), key=lambda x:-sum(x[1].values())):
    print(f"  {g:<30} {dict(c)}")
