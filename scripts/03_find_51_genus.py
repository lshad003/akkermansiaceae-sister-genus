#!/usr/bin/env python3
# Genus field checked for the EHI subset
# Source: ch3-chitin-evolution/scripts/find_51_genus.py
# Output: stdout
import csv, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"
      and r["host_class"]=="amphibian"]
print("amphibian Akkermansiaceae:",len(rows))
c=collections.Counter((r["from_dataset"], (r["genus"] or "EMPTY"), r["family_source"]) for r in rows)
for k,n in c.most_common():
    print(f"  {k[0]:<16}{k[1]:<26}{k[2]:<12}{n}")
miss=[r for r in rows if (r["genus"] or "").strip()=="unknown"]
print(f"\nliteral 'unknown' genus: {len(miss)}")
print("  datasets:", collections.Counter(r["from_dataset"] for r in miss))
print("  examples:", [r["accession"] for r in miss[:5]])
