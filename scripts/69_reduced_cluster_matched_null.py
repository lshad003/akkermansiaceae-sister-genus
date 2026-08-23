#!/usr/bin/env python3
# Loss overlap tested against a prevalence-matched null
# Source: ch3-chitin-evolution/scripts/reduced_cluster_matched_null.py
# Output: results/novel_akk_tree/reduced_cluster_matched_null.txt
import csv, random, statistics, sys
from collections import defaultdict

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
MAT  = f"{BASE}/results/novel_akk_tree/cazy_per_genome_matrix.tsv"
DER  = f"{BASE}/results/akkfam_derep/clusters.tsv"
OUT  = f"{BASE}/results/novel_akk_tree/reduced_cluster_matched_null.txt"
NPERM = 499
random.seed(20260822)

red = None
for r in csv.DictReader(open(DER), delimiter="\t"):
    if r["genus"] == "NOVEL" and int(r["n_genomes"]) == 18:
        red = set(r["members"].split(","))
if red is None:
    print("REFUSED: no NOVEL n=18 cluster"); sys.exit(1)

rows = list(csv.DictReader(open(MAT), delimiter="\t"))
fams = [k for k in rows[0].keys() if k not in ("genome", "group")]
nov  = {r["genome"]: r for r in rows if r["group"] == "novel"}
akkd = {r["genome"]: r for r in rows if r["group"] == "Akkermansia"}

R = [g for g in nov if g in red]
O = [g for g in nov if g not in red]

def pct(sel, f, src):
    return 100.0 * sum(1 for g in sel if int(src[g][f]) > 0) / len(sel)

pR = {f: pct(R, f, nov) for f in fams}
pO = {f: pct(O, f, nov) for f in fams}
pA = {f: pct(list(akkd), f, akkd) for f in fams}

lost_red = set(f for f in fams if pO[f] >= 50.0 and pR[f] <= 10.0)
lost_akk = set(f for f in fams if pO[f] >= 50.0 and pA[f] <= 10.0)
obs = len(lost_red & lost_akk)

# ---- the eligible pool: only families that COULD have been called lost,
# i.e. those present in >=50% of the 87 non-reduced novel genomes.
pool = [f for f in fams if pO[f] >= 50.0]
print("families total: %d" % len(fams))
print("ELIGIBLE pool (present in >=50%% of the 87 non-reduced): %d" % len(pool))
print("lost in reduced cluster: %d   lost on Akkermansia branch: %d   overlap: %d"
      % (len(lost_red), len(lost_akk), obs))
print()
print("prevalence in the 87 non-reduced novel genomes:")
for lab, s in (("reduced-lost", lost_red), ("Akkermansia-lost", lost_akk),
               ("eligible pool", set(pool))):
    v = sorted(pO[f] for f in s)
    print("  %-18s n=%-4d median %.1f%%  range %.1f to %.1f" % (
        lab, len(v), statistics.median(v), v[0], v[-1]))

# ---- NULL 1: uniform over the eligible pool
n1 = []
for _ in range(NPERM):
    d = set(random.sample(pool, len(lost_red)))
    n1.append(len(d & lost_akk))
p1 = (sum(1 for x in n1 if x >= obs) + 1.0) / (NPERM + 1.0)

# ---- NULL 2: prevalence-matched. Bin the eligible pool by prevalence in the
# 87 non-reduced genomes, then draw the same number from each bin as observed.
EDGES = [50, 60, 70, 80, 90, 95, 100.01]
def binof(f):
    for i in range(len(EDGES) - 1):
        if EDGES[i] <= pO[f] < EDGES[i + 1]:
            return i
    return len(EDGES) - 2

bins = defaultdict(list)
for f in pool:
    bins[binof(f)].append(f)
need = defaultdict(int)
for f in lost_red:
    need[binof(f)] += 1

print()
print("prevalence bins (in the 87 non-reduced), eligible pool vs the reduced-lost set:")
short = False
for i in sorted(bins):
    print("  %5.0f-%-5.0f%%  pool %-4d  reduced-lost %-3d" % (
        EDGES[i], EDGES[i + 1] - 0.01, len(bins[i]), need.get(i, 0)))
    if need.get(i, 0) > len(bins[i]):
        short = True
if short:
    print("REFUSED: a bin has fewer eligible families than the observed set requires.")
    sys.exit(1)

n2 = []
for _ in range(NPERM):
    d = set()
    for i, k in need.items():
        d |= set(random.sample(bins[i], k))
    n2.append(len(d & lost_akk))
p2 = (sum(1 for x in n2 if x >= obs) + 1.0) / (NPERM + 1.0)

print()
print("%-34s %-8s %-10s %-8s %s" % ("null", "obs", "null mean", "null max", "p"))
print("%-34s %-8d %-10.2f %-8d %.4f" % ("uniform over eligible pool", obs, statistics.mean(n1), max(n1), p1))
print("%-34s %-8d %-10.2f %-8d %.4f" % ("prevalence-matched", obs, statistics.mean(n2), max(n2), p2))
print()
print("READ IT: if the prevalence-matched p is not significant, the overlap is explained")
print("by both lineages losing their least-common families, which is expected under any")
print("reduction and is NOT a convergence result.")

with open(OUT, "w") as fh:
    fh.write("families_total\t%d\n" % len(fams))
    fh.write("eligible_pool_pO_ge50\t%d\n" % len(pool))
    fh.write("lost_reduced\t%d\n" % len(lost_red))
    fh.write("lost_akkermansia\t%d\n" % len(lost_akk))
    fh.write("overlap_observed\t%d\n" % obs)
    fh.write("null_uniform_mean\t%.3f\nnull_uniform_max\t%d\np_uniform\t%.4f\n"
             % (statistics.mean(n1), max(n1), p1))
    fh.write("null_matched_mean\t%.3f\nnull_matched_max\t%d\np_matched\t%.4f\n"
             % (statistics.mean(n2), max(n2), p2))
    fh.write("n_perm\t%d\nseed\t20260822\n" % NPERM)
print("wrote", OUT)
