#!/usr/bin/env python3
# Non-amphibian genomes excluded by amino acid identity to the type genome
# Source: ch3-chitin-evolution/scripts/tortoise_exclusion_aai.py
# Output: results/tortoise/tortoise_exclusion_aai.tsv
import csv, os, statistics, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TDIR = f"{BASE}/results/tortoise"
OUT  = f"{TDIR}/tortoise_exclusion_aai.tsv"

TYPE = "EHM058340"
QUERIES = ["UHM893.23051_R.bin.46", "UHM893.41075_R.bin.111"]

def load(path):
    if not os.path.exists(path):
        print("MISSING:", path); sys.exit(1)
    best = {}
    nfield = None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if nfield is None:
            nfield = len(f)
            if nfield < 3:
                print("REFUSED: %s has %d columns, expected at least 3 "
                      "(qseqid sseqid pident ...)" % (path, nfield))
                sys.exit(1)
        try:
            pid = float(f[2])
        except ValueError:
            print("REFUSED: column 3 of %s is not numeric: %r" % (path, f[2]))
            sys.exit(1)
        q, s = f[0], f[1]
        if q not in best or pid > best[q][1]:
            best[q] = (s, pid)
    return best

def aai(fwd_path, rev_path):
    fwd = load(fwd_path)
    rev = load(rev_path)
    ids = [pid for q, (s, pid) in fwd.items() if rev.get(s, (None,))[0] == q]
    if not ids:
        return float("nan"), 0, len(fwd), len(rev)
    return statistics.mean(ids), len(ids), len(fwd), len(rev)

rows = []
print("%-30s %-10s %-8s %-10s %-10s" % ("query genome", "AAI_pct", "n_RBH", "n_fwd_hit", "n_rev_hit"))
for q in QUERIES:
    a, n, nf, nr = aai(f"{TDIR}/{q}_fwd.tsv", f"{TDIR}/{q}_rev.tsv")
    print("%-30s %-10.2f %-8d %-10d %-10d" % (q, a, n, nf, nr))
    rows.append([q, "tortoise", TYPE, "%.2f" % a, n])

with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["query_genome", "host_animal_type", "reference_type_genome", "AAI_pct", "n_RBH"])
    for r in rows:
        w.writerow(r)

print()
print("wrote", OUT)
print()
print("INTERPRETATION REQUIRES A WITHIN-GENUS COMPARISON.")
print("These AAI values are only meaningful against the AAI among the 105 amphibian")
print("genomes and the same type genome. That within-genus distribution is NOT computed")
print("here and must come from its own run before any exclusion is claimed.")
print()
print("NOTE: both query genomes come from host animal UHM893, so this test rests on")
print("two bins from ONE tortoise, not two independent reptile hosts.")
