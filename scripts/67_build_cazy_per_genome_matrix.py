#!/usr/bin/env python3
# Enzyme family presence assembled per genome from the annotation output
# Source: ch3-chitin-evolution/scripts/build_cazy_per_genome_matrix.py
# Output: results/novel_akk_tree/cazy_per_genome_matrix.tsv
import csv, os, re, sys
from collections import defaultdict

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUT  = f"{BASE}/results/novel_akk_tree/cazy_per_genome_matrix.tsv"
CNT  = f"{BASE}/results/novel_akk_tree/cazy_per_genome_counts.tsv"

# identical to gh75_verru_census_v3.py lines 46-58
DIRS = ["dbcan_verru", "dbcan_bact_refs", "dbcan_ehi_amphibian", "dbcan_ehi_nonamph",
        "dbcan_endo_allphyla", "dbcan_flavo_refs", "dbcan_scaffold",
        "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"]
EVAL_MAX = 1e-10
COV_MIN  = 0.30
FAM_RE = re.compile(r"^([A-Za-z]+\d+)(_\d+)?\.hmm$")

def bare(a):
    return a.split("_", 1)[1] if a.startswith(("GB_", "RS_")) else a
def variants(a):
    b = bare(a)
    return [a, b, "GB_" + b, "RS_" + b]

idx = {}
for d in DIRS:
    p = d if d.startswith("/") else f"{BASE}/results/{d}"
    if not os.path.isdir(p):
        continue
    for fn in os.listdir(p):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            idx.setdefault(fn[:-4], f"{p}/{fn}")
print("indexed dbCAN tsv files:", len(idx))

rows = list(csv.DictReader(open(CEN), delimiter="\t"))
novel = [r["accession"] for r in rows
         if r["family"] == "Akkermansiaceae"
         and (not r["genus"].strip() or r["genus"] in ("unknown", "NO_GENUS"))
         and r["host_class"] == "amphibian"]
akk = [r["accession"] for r in rows
       if r["family"] == "Akkermansiaceae" and r["genus"] == "Akkermansia"]
print("novel: %d   Akkermansia: %d" % (len(novel), len(akk)))

def families(acc):
    tsv = next((idx[v] for v in variants(acc) if v in idx), None)
    if tsv is None:
        return None
    fams = defaultdict(set)
    with open(tsv) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            m = FAM_RE.match(f[0])
            if not m:
                continue
            try:
                ev = float(f[4]); cv = float(f[9])
            except ValueError:
                continue
            if ev > EVAL_MAX or cv < COV_MIN:
                continue
            fams[m.group(1)].add(f[2])
    return fams

data, missing = {}, []
for group, accs in (("novel", novel), ("Akkermansia", akk)):
    for a in accs:
        fams = families(a)
        if fams is None:
            missing.append(a); continue
        data[a] = (group, fams)
print("genomes with a parsed matrix: %d   missing: %d" % (len(data), len(missing)))
if missing:
    print("missing (first 5):", missing[:5])

allfam = sorted({f for _, fams in data.values() for f in fams})
print("families observed: %d" % len(allfam))

with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome", "group"] + allfam)
    for a in sorted(data):
        g, fams = data[a]
        w.writerow([a, g] + [len(fams.get(f, ())) for f in allfam])

with open(CNT, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome", "group", "n_families", "n_proteins"])
    for a in sorted(data):
        g, fams = data[a]
        w.writerow([a, g, len(fams), sum(len(v) for v in fams.values())])

import statistics
print()
print("%-14s %-6s %-22s %s" % ("group", "n", "median families/genome", "IQR"))
for g in ("novel", "Akkermansia"):
    v = sorted(len(f) for a, (gg, f) in data.items() if gg == g)
    if not v:
        continue
    print("%-14s %-6d %-22.1f %.1f to %.1f" % (
        g, len(v), statistics.median(v), v[len(v)//4], v[3*len(v)//4]))
print()
print("wrote", OUT)
print("wrote", CNT)
print()
print("FILTER: dbCAN hmmscan-parser chain, evalue <= 1e-10, coverage >= 0.30,")
print("identical to gh75_verru_census_v3.py. Counts are PROTEINS per family per genome.")
