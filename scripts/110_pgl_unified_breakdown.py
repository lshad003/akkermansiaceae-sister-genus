#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
import csv, os, sys
from collections import Counter, defaultdict
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
U = B + "/results/pgl_unified"
cen = list(csv.DictReader(open(B + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"), delimiter="\t"))
grp, gen = {}, {}
for r in cen:
    if r["family"] != "Akkermansiaceae": continue
    a = r["accession"]; gen[a] = r["genus"]
    if r["genus"] == "Akkermansia": grp[a] = "Akkermansia"
    elif not r["genus"].strip() or r["genus"] in ("unknown", "NO_GENUS"):
        if r["host_class"] == "amphibian": grp[a] = "candidate"
    else: grp[a] = "free-living"
man = [l.split("\t")[0].strip() for l in open(B + "/results/ppp_unified/manifest.tsv") if l.strip()]
tot = Counter(grp[g] for g in man if g in grp)
print("manifest genomes: %d   by group: %s" % (len(man), dict(tot)))
for name, f in (("pgl", "pgl_hits.tsv"), ("zwf CONTROL", "zwf_hits.tsv")):
    p = U + "/" + f
    if not os.path.exists(p):
        print("\n%s: MISSING %s" % (name, p)); continue
    hit = defaultdict(set)
    for line in open(p):
        g = line.split("\t")[1].split("|", 1)[0].split("::")[-1]
        if g in grp: hit[grp[g]].add(g)
    print("\n%s:" % name)
    for k in ("candidate", "Akkermansia", "free-living"):
        n, t = len(hit.get(k, ())), tot.get(k, 0)
        print("  %-12s %3d / %3d  %5.1f%%" % (k, n, t, 100.0 * n / t if t else 0))
    fl = Counter(gen[g] for g in hit.get("free-living", ()))
    print("  free-living genera with a hit: %d of 17" % len(fl))
    print("   ", dict(fl))
