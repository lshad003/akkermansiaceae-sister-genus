#!/usr/bin/env python3
# GTDB taxonomy strings checked for the EHI subset
# Source: ch3-chitin-evolution/scripts/check_51_gtdb_string.py
# Output: stdout
import csv, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CAT=f"{BASE}/data/unified_annotation/unified_genome_annotation.tsv"
rows=[r for r in csv.DictReader(open(CAT),delimiter="\t")
      if r.get("from_dataset")=="EHI_2025" and r.get("family")=="Akkermansiaceae"]
unk=[r for r in rows if (r.get("genus") or "")=="unknown"]
akk=[r for r in rows if (r.get("genus") or "")=="Akkermansia"]
print(f"EHI Akkermansiaceae: unknown-genus {len(unk)}, Akkermansia {len(akk)}")
print("\n=== full GTDB string, 3 unknown-genus ===")
for r in unk[:3]: print("  ",r["gtdb_taxonomy_full"])
print("=== full GTDB string, 2 Akkermansia ===")
for r in akk[:2]: print("  ",r["gtdb_taxonomy_full"])
print("\n=== g__ slot for all 51 ===")
def g(s):
    for t in (s or "").split(";"):
        t=t.strip()
        if t.startswith("g__"): return t[3:] or "<EMPTY g__>"
    return "<NO g__ FIELD>"
print(collections.Counter(g(r["gtdb_taxonomy_full"]) for r in unk))
print("\n=== the 51: host species, captivity, completeness ===")
print("  hosts:",collections.Counter(r["host_animal_type"] for r in unk))
print("  captivity:",collections.Counter(r["captivity_status"] for r in unk))
comp=sorted(float(r["completeness"]) for r in unk if r["completeness"])
print(f"  completeness median {comp[len(comp)//2]:.1f}  min {comp[0]:.1f}  max {comp[-1]:.1f}")
print("  animals (sample_id):",len(set(r["sample_id"] for r in unk)))
print("\n=== compare: host species of the 204 Akkermansia EHI ===")
print("  hosts:",collections.Counter(r["host_animal_type"] for r in akk).most_common(8))
