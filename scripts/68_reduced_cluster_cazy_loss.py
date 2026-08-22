#!/usr/bin/env python3
# Enzyme families lost in the reduced species cluster
# Source: ch3-chitin-evolution/scripts/reduced_cluster_cazy_loss.py
# Output: results/novel_akk_tree/reduced_cluster_cazy_loss.tsv
import csv, random, statistics, sys
from collections import defaultdict

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
MAT  = f"{BASE}/results/novel_akk_tree/cazy_per_genome_matrix.tsv"
DER  = f"{BASE}/results/akkfam_derep/clusters.tsv"
SGC  = f"{BASE}/results/novel_akk_tree/novel_size_gc.tsv"
OUT  = f"{BASE}/results/novel_akk_tree/reduced_cluster_cazy_loss.tsv"
NPERM = 499
random.seed(20260822)

red = None
for r in csv.DictReader(open(DER), delimiter="\t"):
    if r["genus"] == "NOVEL" and int(r["n_genomes"]) == 18:
        red = set(r["members"].split(","))
if red is None:
    print("REFUSED: no NOVEL n=18 cluster in", DER); sys.exit(1)

sz = {r["genome"]: float(r["genome_size_bp"]) / 1e6
      for r in csv.DictReader(open(SGC), delimiter="\t")}

rows = list(csv.DictReader(open(MAT), delimiter="\t"))
fams = [k for k in rows[0].keys() if k not in ("genome", "group")]

nov = {r["genome"]: r for r in rows if r["group"] == "novel"}
akk = [r for r in rows if r["group"] == "Akkermansia"]

R = [g for g in nov if g in red]
O = [g for g in nov if g not in red]
print("reduced n=%d   other novel n=%d   Akkermansia n=%d" % (len(R), len(O), len(akk)))
print("median size: reduced %.2f Mb, other %.2f Mb" % (
    statistics.median([sz[g] for g in R if g in sz]),
    statistics.median([sz[g] for g in O if g in sz])))

def nfam(r):
    return sum(1 for f in fams if int(r[f]) > 0)
print("median families/genome: reduced %.1f, other novel %.1f, Akkermansia %.1f" % (
    statistics.median([nfam(nov[g]) for g in R]),
    statistics.median([nfam(nov[g]) for g in O]),
    statistics.median([nfam(r) for r in akk])))

def pct(sel, f, src):
    n = sum(1 for g in sel if int(src[g][f]) > 0)
    return 100.0 * n / len(sel)

pR = {f: pct(R, f, nov) for f in fams}
pO = {f: pct(O, f, nov) for f in fams}
akkd = {r["genome"]: r for r in akk}
pA = {f: pct(list(akkd), f, akkd) for f in fams}

# lost in the reduced cluster relative to the rest of the genus
lost_red = set(f for f in fams if pO[f] >= 50.0 and pR[f] <= 10.0)
# lost on the Akkermansia branch relative to the novel genus
lost_akk = set(f for f in fams if pO[f] >= 50.0 and pA[f] <= 10.0)

obs = len(lost_red & lost_akk)
pool = list(fams)
null = []
for _ in range(NPERM):
    d = set(random.sample(pool, len(lost_red)))
    null.append(len(d & lost_akk))
p = (sum(1 for x in null if x >= obs) + 1.0) / (NPERM + 1.0)

print()
print("CRITERION for both sets: present in >=50%% of the 87 non-reduced novel genomes,")
print("<=10%% in the focal set.")
print("families total: %d" % len(fams))
print("lost in reduced cluster: %d" % len(lost_red))
print("lost on Akkermansia branch: %d" % len(lost_akk))
print("overlap: %d" % obs)
print("null over %d perms: mean %.2f max %d" % (NPERM, statistics.mean(null), max(null)))
print("p = %.4f" % p)
print()
print("shared:", ", ".join(sorted(lost_red & lost_akk)) or "none")
print("reduced-only:", ", ".join(sorted(lost_red - lost_akk)) or "none")

with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["family", "pct_reduced_n18", "pct_other_novel_n87", "pct_akkermansia_n346",
                "lost_in_reduced", "lost_in_akkermansia"])
    for f in sorted(fams):
        w.writerow([f, "%.1f" % pR[f], "%.1f" % pO[f], "%.1f" % pA[f],
                    "yes" if f in lost_red else "no",
                    "yes" if f in lost_akk else "no"])
print()
print("wrote", OUT)
