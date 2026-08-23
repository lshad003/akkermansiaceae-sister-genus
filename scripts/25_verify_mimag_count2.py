#!/usr/bin/env python3
# High-quality genome count verified against the MIMAG standard
# Source: ch3-chitin-evolution/scripts/verify_mimag_count2.py
# Output: stdout
import os, re, csv
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
MD=f"{BASE}/results/mimag"
ANN=f"{BASE}/data/unified_annotation/unified_genome_annotation.tsv"
TIPS=f"{BASE}/results/novel_akk_tree/tree_tip_labels.txt"

print("--- sample rrna.gff ---")
for i,l in enumerate(open(f"{MD}/EHM034674.rrna.gff")):
    if i<6: print("  ",l.rstrip()[:140])
print("--- sample trna.txt ---")
for i,l in enumerate(open(f"{MD}/EHM034674.trna.txt")):
    if i<8: print("  ",l.rstrip()[:140])

novel=set(l.strip().split("|")[-1] for l in open(TIPS) if l.strip() and l.split("|")[0]=="NOVEL")
comp={};cont={}
r=csv.reader(open(ANN),delimiter="\t"); h=next(r); ci={c:i for i,c in enumerate(h)}
for row in r:
    a=row[ci["accession"]]
    if a in novel:
        try: comp[a]=float(row[ci["completeness"]]); cont[a]=float(row[ci["contamination"]])
        except: pass

def rrna(acc):
    p=f"{MD}/{acc}.rrna.gff"
    s=set()
    if not os.path.isfile(p): return s
    for l in open(p, errors="ignore"):
        if l.startswith("#"): continue
        m=re.search(r'Name=(\d+S)', l)
        if m: s.add(m.group(1))
    return s

def trna(acc):
    p=f"{MD}/{acc}.trna.txt"
    aas=set()
    if not os.path.isfile(p): return aas
    for l in open(p, errors="ignore"):
        f=l.split()
        if len(f)<5: continue
        if f[0].startswith(("Sequence","Name","----")): continue
        for tok in f[4:6]:
            if tok.isalpha() and len(tok)==3 and tok[0].isupper():
                aas.add(tok); break
    return aas

hq=[]; parsed=0
for acc in sorted(comp):
    rr=rrna(acc); tr=trna(acc)
    if not rr and not tr: continue
    parsed+=1
    if comp[acc]>90 and cont[acc]<5 and {"5S","16S","23S"}<=rr and len(tr)>=18:
        hq.append((acc,comp[acc],cont[acc],sorted(rr),len(tr)))

print(f"\nnovel genomes with output parsed: {parsed} of {len(comp)}")
print(f"MIMAG HIGH-QUALITY: {len(hq)}")
for a,c,ct,rr,nt in hq:
    print(f"   {a}  compl {c}  contam {ct}  rRNA {rr}  tRNA_aa {nt}")
print("\nEHM058340 in HQ set:", any(a=="EHM058340" for a,_,_,_,_ in hq))
e=[x for x in hq if x[0]=="EHM058340"]
if not e:
    print("  EHM058340 detail: rRNA",sorted(rrna("EHM058340")),"tRNA_aa",len(trna("EHM058340")),
          "compl",comp.get("EHM058340"),"contam",cont.get("EHM058340"))
