#!/usr/bin/env python3
# Representative proteins extracted for any polarity category
# Source: ch3-chitin-evolution/scripts/extract_polarity_reps_any.py
# Output: results/pangenome/polarity_reps/shared_gut.faa
import csv, os, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
R    = f"{BASE}/results/pangenome/orthofinder_out/run1/Results_novelakk"
POL  = f"{BASE}/results/pangenome/orthogroup_polarity.tsv"
OGS  = f"{R}/Orthogroups/Orthogroups.tsv"
IN   = f"{BASE}/results/pangenome/orthofinder_in"
OUTD = f"{BASE}/results/pangenome/polarity_reps"

if len(sys.argv) != 2:
    print("usage: extract_polarity_reps_any.py <category>")
    print("categories in the polarity file:")
    import collections
    for k, v in collections.Counter(
            r["category"] for r in csv.DictReader(open(POL), delimiter="\t")).most_common():
        print("   %-26s %d" % (k, v))
    sys.exit(1)
CAT = sys.argv[1]

# identical selection logic to scripts/extract_polarity_reps.py
want = set(row["orthogroup"] for row in csv.DictReader(open(POL), delimiter="\t")
           if row["category"] == CAT)
if not want:
    print("REFUSED: no orthogroups in category %r" % CAT); sys.exit(1)
print("orthogroups in %s: %d" % (CAT, len(want)))

og2gene, need = {}, set()
r = csv.reader(open(OGS), delimiter="\t")
next(r)
for row in r:
    og = row[0]
    if og not in want:
        continue
    allog = set()
    for c in row[1:]:
        allog |= set(x.strip() for x in c.split(",") if x.strip())
    if allog:
        g = sorted(allog)[0]
        og2gene[og] = g
        need.add(g)
print("representatives selected:", len(og2gene))

seqs, cur = {}, None
for fn in sorted(os.listdir(IN)):
    if not fn.endswith(".faa"):
        continue
    for line in open(os.path.join(IN, fn)):
        if line.startswith(">"):
            cur = line[1:].split()[0]
            if cur in need:
                seqs[cur] = [line]
        elif cur in seqs:
            seqs[cur].append(line)
print("sequences recovered:", len(seqs))

p = f"{OUTD}/{CAT}.faa"
if os.path.exists(p):
    print("REFUSED: %s already exists, refusing to overwrite" % p); sys.exit(1)
n = 0
with open(p, "w") as fh:
    for og, g in sorted(og2gene.items()):
        if g not in seqs:
            continue
        fh.write(">%s\n" % og)
        fh.write("".join(seqs[g][1:]))
        n += 1
print("wrote %s  (%d sequences of %d orthogroups)" % (p, n, len(want)))
