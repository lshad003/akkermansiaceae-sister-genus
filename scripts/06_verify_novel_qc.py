#!/usr/bin/env python3
# Completeness and contamination verified for the candidate set
# Source: ch3-chitin-evolution/scripts/verify_novel_qc.py
# Output: stdout
import csv, statistics

ANN="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/unified_annotation/unified_genome_annotation.tsv"
TIPS="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_tree/tree_tip_labels.txt"

novel=[]
for l in open(TIPS):
    l=l.strip()
    if l and l.split("|")[0]=="NOVEL":
        novel.append(l.split("|")[-1])
print("NOVEL tips found:", len(novel))

rows={}
r=csv.reader(open(ANN),delimiter="\t")
hdr=next(r); ci={c:i for i,c in enumerate(hdr)}
for row in r:
    rows[row[ci["accession"]]]=row

def match(acc):
    if acc in rows: return acc
    if acc[:3] in ("GB_","RS_") and acc[3:] in rows: return acc[3:]
    for pre in ("GB_","RS_"):
        if pre+acc in rows: return pre+acc
    return None

matched=[]; unmatched=[]
for a in novel:
    m=match(a)
    (matched if m else unmatched).append(m or a)
print("matched:", len(matched), " unmatched:", len(unmatched))
if unmatched: print("  unmatched examples:", unmatched[:6])

comps=[]; animals=set(); by_ds={}
for m in matched:
    row=rows[m]; ds=row[ci["from_dataset"]]
    by_ds[ds]=by_ds.get(ds,0)+1
    try: comps.append(float(row[ci["completeness"]]))
    except: pass
    acc=row[ci["accession"]]
    unit=acc.split(".")[0] if ds=="herptile_MAG" else (row[ci["sample_id"]] or acc)
    animals.add((ds,unit))

print("by from_dataset:", by_ds)
print("host animals (corrected unit):", len(animals))
if comps:
    print("completeness: n=%d median=%.1f min=%.1f max=%.1f" % (
        len(comps), statistics.median(comps), min(comps), max(comps)))
