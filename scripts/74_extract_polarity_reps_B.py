#!/usr/bin/env python3
# Representative proteins extracted for the ancestral losses
# Source: ch3-chitin-evolution/scripts/extract_polarity_reps_B.py
# Output: results/pangenome/polarity_reps/B_lost_at_gut_ancestor.faa
import csv, os, collections

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
R    = f"{BASE}/results/pangenome/orthofinder_out/run1/Results_novelakk"
POL  = f"{BASE}/results/pangenome/orthogroup_polarity.tsv"
OGS  = f"{R}/Orthogroups/Orthogroups.tsv"
IN   = f"{BASE}/results/pangenome/orthofinder_in"
OUTD = f"{BASE}/results/pangenome/polarity_reps"
CAT  = "B_lost_at_gut_ancestor"

# same selection logic as scripts/extract_polarity_reps.py, category B only.
want = set(row["orthogroup"] for row in csv.DictReader(open(POL), delimiter="\t")
           if row["category"] == CAT)
print("orthogroups in %s: %d" % (CAT, len(want)))

og2gene = {}
need = set()
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

seqs = {}
cur = None
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
n = 0
with open(p, "w") as fh:
    for og, g in sorted(og2gene.items()):
        if g not in seqs:
            continue
        fh.write(">%s\n" % og)
        fh.write("".join(seqs[g][1:]))
        n += 1
print("wrote %s  (%d sequences)" % (p, n))
if n < 0.9 * len(want):
    print("NOTE: %d of %d orthogroups had no recoverable sequence." % (len(want) - n, len(want)))
