#!/usr/bin/env python3
# Representative proteins extracted for annotation
# Source: ch3-chitin-evolution/scripts/extract_polarity_reps.py
# Output: results/pangenome/polarity_reps
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
R=f"{BASE}/results/pangenome/orthofinder_out/run1/Results_novelakk"
POL=f"{BASE}/results/pangenome/orthogroup_polarity.tsv"
OGS=f"{R}/Orthogroups/Orthogroups.tsv"
IN=f"{BASE}/results/pangenome/orthofinder_in"
OUTD=f"{BASE}/results/pangenome/polarity_reps"
os.makedirs(OUTD,exist_ok=True)

want=collections.defaultdict(set)
for row in csv.DictReader(open(POL),delimiter="\t"):
    c=row["category"]
    if c in ("A_lost_in_Akkermansia","C_Akkermansia_enriched","D_novel_specific"):
        want[c].add(row["orthogroup"])
print({k:len(v) for k,v in want.items()})

need=set()
og2gene={}
r=csv.reader(open(OGS),delimiter="\t"); hdr=next(r)
for row in r:
    og=row[0]
    allog=set().union(*[set(x.strip() for x in c.split(",") if x.strip()) for c in row[1:]]) if len(row)>1 else set()
    for c,s in want.items():
        if og in s and allog:
            g=sorted(allog)[0]
            og2gene[og]=(c,g); need.add(g)
print("representatives selected:",len(og2gene))

seqs={}
cur=None
for fn in os.listdir(IN):
    if not fn.endswith(".faa"): continue
    for line in open(os.path.join(IN,fn)):
        if line.startswith(">"):
            cur=line[1:].split()[0]
            keep = cur in need
            if keep: seqs[cur]=[line]
        elif cur in seqs:
            seqs[cur].append(line)

files={}
for og,(c,g) in og2gene.items():
    if g not in seqs: continue
    files.setdefault(c,[]).append((og,seqs[g]))
for c,items in files.items():
    p=f"{OUTD}/{c}.faa"
    with open(p,"w") as fh:
        for og,lines in items:
            fh.write(f">{og}\n")
            fh.write("".join(lines[1:]))
    print(f"wrote {p}  ({len(items)} sequences)")

print("\n--- annotation tools available ---")
