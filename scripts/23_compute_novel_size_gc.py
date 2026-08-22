#!/usr/bin/env python3
# Genome size and GC computed from assemblies
# Source: ch3-chitin-evolution/scripts/compute_novel_size_gc.py
# Output: results/novel_akk_tree/novel_size_gc.tsv
import os,csv
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
LIST=f"{BASE}/results/akk_composition/novel_akk_genome_list.txt"
FADIR=f"{BASE}/data/amphibia_gtdbtk_input"
OUT=f"{BASE}/results/novel_akk_tree/novel_size_gc.tsv"

def stats(fa):
    size=0; gc=0
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"): continue
            s=line.strip().upper()
            size+=len(s); gc+=s.count("G")+s.count("C")
    return size, (100.0*gc/size if size else 0.0)

paths=[l.strip() for l in open(LIST) if l.strip()]
rows=[]
miss=[]
for p in paths:
    fa=p if os.path.isabs(p) else os.path.join(FADIR,p)
    if not os.path.isfile(fa):
        base=os.path.basename(p)
        fa=os.path.join(FADIR,base)
    if not os.path.isfile(fa):
        miss.append(p); continue
    gid=os.path.basename(fa)[:-3] if fa.endswith(".fa") else os.path.basename(fa)
    sz,gc=stats(fa)
    rows.append((gid,sz,round(gc,2)))

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["genome","genome_size_bp","gc_percent"])
    w.writerows(rows)

sizes=sorted(r[1] for r in rows); gcs=sorted(r[2] for r in rows)
def med(x): 
    n=len(x); return x[n//2] if n%2 else (x[n//2-1]+x[n//2])/2
print(f"computed {len(rows)}/105 genomes  (missing {len(miss)})")
if miss: print("  MISSING:",miss[:5])
print(f"  size median {med(sizes)/1e6:.2f} Mb  range {min(sizes)/1e6:.2f}-{max(sizes)/1e6:.2f}")
print(f"  GC median {med(gcs):.1f}%  range {min(gcs):.1f}-{max(gcs):.1f}")
print(f"wrote {OUT}")
