#!/usr/bin/env python3
# Denominator = the searched database. Taxonomy joined on version-stripped accession.
import os, sys
from collections import defaultdict, Counter

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W   = R + "/results/ppp_baserate"
CAT = R + "/data/unified_annotation/unified_genome_annotation.tsv"
CEN = R + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUT = W + "/ppp_baserate_v3.txt"

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

QN = {"UHM979.41089_R.bin.103_CDS_0654": "zwf",
      "UHM1210.23070_R.bin.101_CDS_0293": "gnd",
      "EHM058980_CDS_1432": "g6pd_sub"}

def load(p):
    d = defaultdict(set)
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            d[f[1].split("|")[0]].add(QN.get(f[0], f[0]))
    return d
ppp = load(W + "/ppp_trio.tsv")

indb = set()
for line in open(W + "/all.faa"):
    if line.startswith(">"):
        indb.add(line[1:].split("|", 1)[0])
say("database genomes: %d   with a hit: %d" % (len(indb), len(ppp)))

def key(a):
    a = a.strip()
    for p in ("GB_", "RS_"):
        if a.startswith(p): a = a[3:]
    return a.split(".")[0]

# taxonomy from both sources, version-stripped
fam = {}
for path, accfield, famfield in ((CEN, 0, None), (CAT, None, None)):
    if not os.path.exists(path): continue
    with open(path) as f:
        h = f.readline().rstrip("\n").split("\t")
        ia = 0 if accfield == 0 else (h.index("accession") if "accession" in h else 0)
        i_f = next((h.index(c) for c in ("family", "Family") if c in h), None)
        if i_f is None: continue
        for line in f:
            p = line.rstrip("\n").split("\t")
            if ia >= len(p) or i_f >= len(p): continue
            v = p[i_f].strip()
            if v.startswith("f__"): v = v[3:]
            k = key(p[ia])
            if v and k not in fam: fam[k] = v
say("taxonomy keys loaded: %d" % len(fam))

tot = Counter(); has = defaultdict(Counter); unk = 0
for g in indb:
    fm = fam.get(key(g))
    if fm is None:
        unk += 1; fm = "(no taxonomy)"
    tot[fm] += 1
    s = ppp.get(g, set())
    if "zwf" in s: has[fm]["zwf"] += 1
    if "gnd" in s: has[fm]["gnd"] += 1
    if "zwf" in s and "gnd" in s: has[fm]["both"] += 1
say("database genomes without taxonomy: %d" % unk)

say(""); say("=" * 70)
say("OXIDATIVE PPP ACROSS VERRUCOMICROBIOTA (diamond blastp, evalue 1e-5)")
say("=" * 70)
say("  %-30s %6s %8s %8s %8s" % ("family", "n", "zwf%", "gnd%", "both%"))
for fm, n in tot.most_common():
    if n < 25: continue
    say("  %-30s %6d %7.1f%% %7.1f%% %7.1f%%"
        % (fm[:30], n, 100.0*has[fm]["zwf"]/n, 100.0*has[fm]["gnd"]/n, 100.0*has[fm]["both"]/n))

N = sum(tot.values())
say("")
for k in ("zwf", "gnd", "both"):
    v = sum(has[f][k] for f in tot)
    say("  ALL  %-5s %5d / %d = %.1f%%" % (k, v, N, 100.0*v/N))
say("")
say("Compare: 0 of 187 Akkermansia carry any of the three.")
say(""); say("written to " + OUT)
fh.close()
