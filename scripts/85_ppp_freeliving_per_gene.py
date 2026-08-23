#!/usr/bin/env python3
# Pathway prevalence per gene across the free-living genera
# Source: ch3-chitin-evolution/scripts/ppp_freeliving_per_gene.py
# Output: results/pangenome/ppp_freeliving_per_gene.tsv
import os
import collections

PG = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/"
HITS = PG + "ppp_outgroups.tsv"
OUTDIR = PG + "outgroups"
OUT = PG + "ppp_freeliving_per_gene.tsv"

# Identified by sequence motif in ppp_trio.faa, not assumed:
#   CDS_0654  518 aa  VIFGATGDL   = G6PD Rossmann motif        -> zwf
#   CDS_0293  476 aa  GLVGLAVMGRNL = 6PGD NADP-binding motif   -> gnd
#   CDS_1432  367 aa  neither motif, by elimination            -> g6pd_sub
QGENE = {
    "UHM979.41089_R.bin.103_CDS_0654": "zwf",
    "UHM1210.23070_R.bin.101_CDS_0293": "gnd",
    "EHM058980_CDS_1432": "g6pd_sub",
}

acc2genus = {}
rows = []
for line in open(HITS):
    f = line.rstrip("\n").split("\t")
    if len(f) < 5 or "::" not in f[1]:
        continue
    genus, rest = f[1].split("::", 1)
    acc = rest.split("|")[0]
    acc2genus[acc] = genus
    rows.append((f[0], genus, acc))

faa = [x[:-4] for x in os.listdir(OUTDIR) if x.endswith(".faa")]
denom = collections.Counter(acc2genus[a] for a in faa if a in acc2genus)

hit = collections.defaultdict(set)
for q, genus, acc in rows:
    hit[(QGENE.get(q, "UNMAPPED"), genus)].add(acc)

gen = sorted(denom)
lines = ["gene\t" + "\t".join(gen) + "\ttotal"]
print("%-10s %s %s" % ("gene", "  ".join("%-16s" % g for g in gen), "total"))
for gene in ("zwf", "gnd", "g6pd_sub"):
    cells, tn, td = [], 0, 0
    for g in gen:
        n = len(hit.get((gene, g), ()))
        cells.append("%d/%d" % (n, denom[g]))
        tn += n
        td += denom[g]
    print("%-10s %s %d/%d" % (gene, "  ".join("%-16s" % c for c in cells), tn, td))
    lines.append(gene + "\t" + "\t".join(cells) + "\t%d/%d" % (tn, td))

open(OUT, "w").write("\n".join(lines) + "\n")
print("")
print("WROTE:", OUT)
