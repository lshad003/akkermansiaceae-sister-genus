#!/usr/bin/env python3
# Type genome checked against the set median
# Source: ch3-chitin-evolution/scripts/check_type_genome_median.py
# Output: stdout
import csv, statistics

F="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_tree/novel_size_gc.tsv"
TYPE="EHM058340"

rows=list(csv.DictReader(open(F),delimiter="\t"))
print("columns:", rows[0].keys() if rows else "EMPTY")
print("n rows:", len(rows))

def pick(keys, *cands):
    for c in cands:
        for k in keys:
            if c in k.lower(): return k
    return None
keys=list(rows[0].keys())
kacc=pick(keys,"accession","genome","id")
ksize=pick(keys,"size","length","bp","mb")
kgc=pick(keys,"gc")
print("using columns -> acc:",kacc," size:",ksize," gc:",kgc)

sizes=[]; gcs=[]; recs=[]
for r in rows:
    try:
        s=float(r[ksize]); g=float(r[kgc])
    except: continue
    sizes.append(s); gcs.append(g); recs.append((r[kacc],s,g))

ms=statistics.median(sizes); mg=statistics.median(gcs)
print("\nmedian size: %.4f   median GC: %.3f   (n=%d)" % (ms,mg,len(recs)))
print("size range: %.4f - %.4f | GC range: %.2f - %.2f" % (min(sizes),max(sizes),min(gcs),max(gcs)))

# rank every genome by combined normalized distance to the medians
srange=(max(sizes)-min(sizes)) or 1.0
grange=(max(gcs)-min(gcs)) or 1.0
scored=[]
for a,s,g in recs:
    d=abs(s-ms)/srange + abs(g-mg)/grange
    scored.append((d,a,s,g))
scored.sort()

print("\nTOP 10 CLOSEST TO MEDIAN (combined normalized distance):")
for i,(d,a,s,g) in enumerate(scored[:10],1):
    tag=" <-- stated type" if TYPE in a else ""
    print("  %2d. %-22s size=%.3f gc=%.2f  dist=%.4f%s" % (i,a,s,g,d,tag))

hit=[(i,d,a,s,g) for i,(d,a,s,g) in enumerate(scored,1) if TYPE in a]
print()
if hit:
    i,d,a,s,g = hit[0]
    print("TYPE %s: rank %d of %d | size=%.4f (median %.4f, diff %+.4f) | gc=%.2f (median %.2f, diff %+.2f)"
          % (a,i,len(scored),s,ms,s-ms,g,mg,g-mg))
    print("VERDICT:", "closest-to-median claim OK" if i<=3 else
          "NOT closest to median -- rewrite the type-selection sentence")
else:
    print("TYPE", TYPE, "NOT FOUND in this file. Check the accession or the file.")
