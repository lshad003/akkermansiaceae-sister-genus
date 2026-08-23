#!/usr/bin/env python3
# Wild and captive status reconciled across metadata sources
# Source: ch3-chitin-evolution/scripts/reconcile_captivity.py
# Output: stdout
import csv
ANN="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/unified_annotation/unified_genome_annotation.tsv"
SRC="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/akkermansia_animal_source.tsv"
TIPS="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_tree/tree_tip_labels.txt"

novel=set(l.strip().split("|")[-1] for l in open(TIPS) if l.strip() and l.split("|")[0]=="NOVEL")
ann={}
r=csv.reader(open(ANN),delimiter="\t"); h=next(r); ci={c:i for i,c in enumerate(h)}
for row in r: ann[row[ci["accession"]]]=row
def m(a):
    if a in ann: return a
    if a[:3] in ("GB_","RS_") and a[3:] in ann: return a[3:]
    return None
src={}
for row in csv.DictReader(open(SRC),delimiter="\t"):
    src[row["mag"]]=row["source"]

print("novel genomes where catalogue and deposit-source DISAGREE:")
dis=0
for a in sorted(novel):
    k=m(a)
    cat = ann[k][ci["captivity_status"]] if k else "?"
    dep = src.get(a,"(not in src file)")
    cat_captive = cat.lower()=="captive"
    dep_captive = dep in ("WF22","WF23","WF24","zoo","captive")
    if cat_captive != dep_captive:
        dis+=1
        print(f"  {a}  catalogue={cat}  deposit={dep}")
print(f"\ndisagreements: {dis} of {len(novel)} novel genomes")
