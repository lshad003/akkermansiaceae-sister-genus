#!/usr/bin/env python3
# Enzyme repertoire of the reduced cluster tested against genome size
# Source: ch3-chitin-evolution/scripts/reduced_cluster_proportionality.py
# Output: results/novel_akk_tree/reduced_cluster_proportionality.txt
import csv, statistics, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CNT  = f"{BASE}/results/novel_akk_tree/cazy_per_genome_counts.tsv"
SGC  = f"{BASE}/results/novel_akk_tree/novel_size_gc.tsv"
DER  = f"{BASE}/results/akkfam_derep/clusters.tsv"
OUT  = f"{BASE}/results/novel_akk_tree/reduced_cluster_proportionality.txt"

red = None
for r in csv.DictReader(open(DER), delimiter="\t"):
    if r["genus"] == "NOVEL" and int(r["n_genomes"]) == 18:
        red = set(r["members"].split(","))
if red is None:
    print("REFUSED: no NOVEL n=18 cluster"); sys.exit(1)

sz = {r["genome"]: float(r["genome_size_bp"]) / 1e6
      for r in csv.DictReader(open(SGC), delimiter="\t")}

fam, prot = {}, {}
for r in csv.DictReader(open(CNT), delimiter="\t"):
    if r["group"] != "novel":
        continue
    fam[r["genome"]] = int(r["n_families"])
    prot[r["genome"]] = int(r["n_proteins"])

both = [g for g in fam if g in sz]
R = [g for g in both if g in red]
O = [g for g in both if g not in red]
print("reduced n=%d   other novel n=%d" % (len(R), len(O)))
if len(R) < 5 or len(O) < 5:
    print("REFUSED: too few genomes with both size and family counts"); sys.exit(1)

def med(v): return statistics.median(v)

mR_sz, mO_sz = med([sz[g] for g in R]), med([sz[g] for g in O])
mR_fa, mO_fa = med([fam[g] for g in R]), med([fam[g] for g in O])
mR_pr, mO_pr = med([prot[g] for g in R]), med([prot[g] for g in O])

print()
print("%-26s %-12s %-12s %s" % ("", "reduced", "other novel", "reduced as % of other"))
for lab, a, b in (("genome size (Mb)", mR_sz, mO_sz),
                  ("CAZyme families", mR_fa, mO_fa),
                  ("CAZyme proteins", mR_pr, mO_pr)):
    print("%-26s %-12.2f %-12.2f %.1f%%" % (lab, a, b, 100.0 * a / b))

print()
print("DENSITY, the proportionality test:")
dR_fa = [fam[g] / sz[g] for g in R]
dO_fa = [fam[g] / sz[g] for g in O]
dR_pr = [prot[g] / sz[g] for g in R]
dO_pr = [prot[g] / sz[g] for g in O]
print("%-26s %-14s %-14s" % ("", "reduced", "other novel"))
print("%-26s %-14.2f %-14.2f" % ("families per Mb", med(dR_fa), med(dO_fa)))
print("%-26s %-14.2f %-14.2f" % ("CAZyme proteins per Mb", med(dR_pr), med(dO_pr)))

# Mann-Whitney U on families per Mb, exact-ish via rank sum
def mwu_p(a, b):
    allv = [(v, 0) for v in a] + [(v, 1) for v in b]
    allv.sort()
    ranks = {}
    i = 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks.setdefault(k, r)
        i = j + 1
    R1 = sum(ranks[k] for k in range(len(allv)) if allv[k][1] == 0)
    n1, n2 = len(a), len(b)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    U = min(U1, n1 * n2 - U1)
    mu = n1 * n2 / 2.0
    sd = (n1 * n2 * (n1 + n2 + 1) / 12.0) ** 0.5
    if sd == 0:
        return float("nan"), U
    z = (U - mu) / sd
    # two-sided normal approximation
    import math
    p = math.erfc(abs(z) / (2 ** 0.5))
    return p, U

p_fa, _ = mwu_p(dR_fa, dO_fa)
p_pr, _ = mwu_p(dR_pr, dO_pr)
print()
print("Mann-Whitney U, two-sided, normal approximation:")
print("  families per Mb        p = %.4f" % p_fa)
print("  CAZyme proteins per Mb p = %.4f" % p_pr)

print()
print("READ IT: if reduced CAZyme families as a percentage of the other novel genomes is close")
print("to the genome-size percentage, the CAZyme loss simply tracks genome size and is not a")
print("targeted reduction. If families per Mb is significantly LOWER in the reduced cluster,")
print("CAZymes were lost DISPROPORTIONATELY, which is the stronger claim.")

with open(OUT, "w") as fh:
    fh.write("n_reduced\t%d\nn_other_novel\t%d\n" % (len(R), len(O)))
    fh.write("median_size_Mb_reduced\t%.3f\nmedian_size_Mb_other\t%.3f\n" % (mR_sz, mO_sz))
    fh.write("median_families_reduced\t%.1f\nmedian_families_other\t%.1f\n" % (mR_fa, mO_fa))
    fh.write("median_proteins_reduced\t%.1f\nmedian_proteins_other\t%.1f\n" % (mR_pr, mO_pr))
    fh.write("size_pct\t%.1f\nfamilies_pct\t%.1f\nproteins_pct\t%.1f\n"
             % (100.0*mR_sz/mO_sz, 100.0*mR_fa/mO_fa, 100.0*mR_pr/mO_pr))
    fh.write("median_families_per_Mb_reduced\t%.3f\nmedian_families_per_Mb_other\t%.3f\n"
             % (med(dR_fa), med(dO_fa)))
    fh.write("median_proteins_per_Mb_reduced\t%.3f\nmedian_proteins_per_Mb_other\t%.3f\n"
             % (med(dR_pr), med(dO_pr)))
    fh.write("p_families_per_Mb\t%.4f\np_proteins_per_Mb\t%.4f\n" % (p_fa, p_pr))
print()
print("wrote", OUT)
