#!/usr/bin/env python3
# One table: candidate genus, Akkermansia, free-living Akkermansiaceae, all from one search.
import os, sys
from collections import defaultdict, Counter

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W   = R + "/results/ppp_unified"
CEN = R + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
NOV = R + "/results/novel_akk_tree/novel_size_gc.tsv"
OUT = W + "/ppp_unified.txt"

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

QN = {"UHM979.41089_R.bin.103_CDS_0654": "zwf",
      "UHM1210.23070_R.bin.101_CDS_0293": "gnd",
      "EHM058980_CDS_1432": "g6pd_sub"}

def load(p):
    d = defaultdict(set)
    if not os.path.exists(p): return d
    for line in open(p):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2:
            d[f[1].split("|")[0]].add(QN.get(f[0], f[0]))
    return d
ppp, ctl = load(W + "/ppp_trio.tsv"), load(W + "/control_gh20.tsv")

indb = set()
for line in open(W + "/all.faa"):
    if line.startswith(">"):
        indb.add(line[1:].split("|", 1)[0])
say("genomes searched: %d" % len(indb))

with open(NOV) as f:
    f.readline(); cand = set(l.split("\t")[0].strip() for l in f if l.strip())

fam, gen = {}, {}
with open(CEN) as f:
    h = f.readline().rstrip("\n").split("\t")
    ifm, ign = h.index("family"), h.index("genus")
    for line in f:
        p = line.rstrip("\n").split("\t")
        if ifm >= len(p): continue
        s = lambda v: v[3:] if v[:3] in ("f__", "g__") else v
        fam[p[0].strip()] = s(p[ifm].strip())
        gen[p[0].strip()] = s(p[ign].strip()) if ign < len(p) else ""

def group(g):
    if g in cand: return "1. CANDIDATE GENUS"
    if gen.get(g) == "Akkermansia": return "2. Akkermansia"
    if fam.get(g) == "Akkermansiaceae": return "3. other Akkermansiaceae"
    return None

tot = Counter(); has = defaultdict(Counter); bygen = defaultdict(Counter)
for g in indb:
    grp = group(g)
    if grp is None: continue
    tot[grp] += 1
    s = ppp.get(g, set())
    for k in ("zwf", "gnd", "g6pd_sub"):
        if k in s: has[grp][k] += 1
    if {"zwf", "gnd"} <= s: has[grp]["both"] += 1
    if ctl.get(g): has[grp]["gh20"] += 1
    if grp == "3. other Akkermansiaceae":
        bygen[gen.get(g) or "(unassigned)"][bool(s & {"zwf", "gnd"})] += 1

say(""); say("=" * 78)
say("OXIDATIVE PPP, SINGLE SEARCH, diamond blastp evalue 1e-5")
say("=" * 78)
say("  %-26s %6s %8s %8s %8s %8s %8s" % ("group", "n", "zwf", "gnd", "sub", "zwf+gnd", "GH20"))
for grp in sorted(tot):
    n = tot[grp]
    say("  %-26s %6d %7d %8d %8d %8d %8d"
        % (grp, n, has[grp]["zwf"], has[grp]["gnd"], has[grp]["g6pd_sub"],
           has[grp]["both"], has[grp]["gh20"]))
say("")
say("  %-26s %6s %8s %8s %8s %8s %8s" % ("", "", "zwf%", "gnd%", "sub%", "both%", "GH20%"))
for grp in sorted(tot):
    n = tot[grp]
    say("  %-26s %6d %7.1f%% %7.1f%% %7.1f%% %7.1f%% %7.1f%%"
        % (grp, n, 100.0*has[grp]["zwf"]/n, 100.0*has[grp]["gnd"]/n,
           100.0*has[grp]["g6pd_sub"]/n, 100.0*has[grp]["both"]/n, 100.0*has[grp]["gh20"]/n))

say(""); say("other Akkermansiaceae, by genus:")
say("  %-22s %6s %6s" % ("genus", "hit", "nohit"))
for g in sorted(bygen):
    say("  %-22s %6d %6d" % (g[:22], bygen[g][True], bygen[g][False]))

say("")
say("GH20 is the positive control: it must be high in groups 1 and 2.")
say("If GH20 is high and zwf/gnd is zero in Akkermansia, absence is real.")
say(""); say("written to " + OUT)
fh.close()
