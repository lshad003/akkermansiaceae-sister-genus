#!/usr/bin/env python3
# Collection denominators verified against database directories
# Source: ch3-chitin-evolution/scripts/check_akk_denominator.py
# Output: stdout
import csv, collections
CEN="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
rows=[r for r in csv.DictReader(open(CEN),delimiter="\t") if r["annotated"]=="1"]
akk=[r for r in rows if r["family"]=="Akkermansiaceae" and (r["genus"] or "").strip()=="Akkermansia"]
print("Akkermansia annotated total:",len(akk))
hc=collections.Counter(r["host_class"] or "(blank)" for r in akk)
print("by host_class:")
for k,v in hc.most_common(): print(f"   {v:4d}  {k}")
print("\nsum of the four groups you cite (amph94 + Podarcis137 + mammal71 + GTDB42):",94+137+71+42)
print("difference:",len(akk)-(94+137+71+42))
