#!/usr/bin/env python3
# Species cluster representatives selected by completeness
# Source: ch3-chitin-evolution/scripts/pick_cluster_reps.py
# Output: results/cluster_aai/cluster_representatives.tsv
import csv

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
CLU = CH3 + "results/akkfam_derep/clusters.tsv"
UA = CH3 + "data/unified_annotation/unified_genome_annotation.tsv"
OUT = CH3 + "results/cluster_aai/cluster_representatives.tsv"

g2c = {}
size = {}
for r in csv.DictReader(open(CLU), delimiter="\t"):
    if r["genus"] not in ("NOVEL", "novel"):
        continue
    size[r["cluster_id"]] = int(r["n_genomes"])
    for m in r["members"].replace(",", " ").split():
        g2c[m.strip()] = r["cluster_id"]

comp = {}
for r in csv.DictReader(open(UA), delimiter="\t"):
    if r["accession"] in g2c:
        try:
            comp[r["accession"]] = float(r["completeness"])
        except (ValueError, TypeError):
            pass

best = {}
for g, c in g2c.items():
    v = comp.get(g, -1.0)
    if c not in best or v > best[c][1]:
        best[c] = (g, v)

print("NOVEL clusters: %d   genomes: %d" % (len(size), sum(size.values())))
print("")
print("%-8s %5s %-32s %s" % ("cluster", "n", "representative", "completeness"))
lines = ["cluster_id\tn_genomes\trepresentative\tcompleteness"]
for c in sorted(best):
    g, v = best[c]
    print("%-8s %5d %-32s %.2f" % (c, size[c], g, v))
    lines.append("%s\t%d\t%s\t%.2f" % (c, size[c], g, v))

open(OUT, "w").write("\n".join(lines) + "\n")
print("")
print("WROTE:", OUT)
print("Representative = highest CheckM completeness within each cluster.")
print("Lowest-completeness representative carries a caveat in the README.")
