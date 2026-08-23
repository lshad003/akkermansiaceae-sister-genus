#!/usr/bin/env python3
# Operon adjacency, intergenic gap and co-orientation written to file
# Source: ch3-chitin-evolution/scripts/ppp_operon_tee.py
# Output: results/pangenome/ppp_operon_verified.txt
import os
import re

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P = BASE + "/results/pangenome"
GFF = P + "/gff199"
OUT = P + "/ppp_operon_verified.txt"

zwf_q = "UHM979.41089_R.bin.103_CDS_0654"
g6pd_q = "EHM058980_CDS_1432"

L = []
def say(s):
    print(s)
    L.append(s)

say("PPP OPERON, zwf and the G6PD accessory subunit")
say("query zwf      = %s  (518 aa, VIFGATGDL Rossmann motif)" % zwf_q)
say("query g6pd_sub = %s  (367 aa, accessory subunit)" % g6pd_q)
say("identity floor 30.0, best hit per genome, from ppp_all199.tsv")
say("")

def subj(qid, pid=30.0):
    best = {}
    for l in open(P + "/ppp_all199.tsv", errors="ignore"):
        x = l.rstrip("\n").split("\t")
        if len(x) < 3 or x[0] != qid:
            continue
        try:
            p = float(x[2])
        except ValueError:
            continue
        if p < pid:
            continue
        s = x[1]
        g = s.split("|")[0].split("_CDS_")[0]
        if g not in best:
            best[g] = s
    return best

zwf = subj(zwf_q)
g6pd = subj(g6pd_q)

def contig_gene(s):
    tail = s.split("|")[1]
    m = re.match(r"(.+)_(\d+)$", tail)
    return (m.group(1), int(m.group(2))) if m else (None, None)

def gff_map(genome):
    fp = "%s/%s.gff" % (GFF, genome)
    d = {}
    if not os.path.isfile(fp):
        return d
    for l in open(fp, errors="ignore"):
        if "\tCDS\t" not in l:
            continue
        f = l.split("\t")
        m = re.search(r"ID=\d+_(\d+)", f[8])
        if m:
            d[(f[0], int(m.group(1)))] = (int(f[3]), int(f[4]), f[6])
    return d

gaps = []
same = 0
both = 0
co = 0
opp = 0
strand_n = 0

for genome in sorted(set(zwf) & set(g6pd)):
    both += 1
    zc, zi = contig_gene(zwf[genome])
    gc, gi = contig_gene(g6pd[genome])
    if zc != gc or zc is None:
        continue
    same += 1
    gm = gff_map(genome)
    z = gm.get((zc, zi))
    g = gm.get((gc, gi))
    if z and g:
        gaps.append(abs(max(g[0] - z[1], z[0] - g[1])))
        strand_n += 1
        if z[2] == g[2]:
            co += 1
        else:
            opp += 1

say("genomes carrying both genes: %d" % both)
say("same contig: %d of %d" % (same, both))
if gaps:
    gaps.sort()
    n = len(gaps)
    med = gaps[n // 2] if n % 2 else (gaps[n // 2 - 1] + gaps[n // 2]) // 2
    say("intergenic gap: median %d bp, range %d to %d, n = %d" % (med, min(gaps), max(gaps), n))
    say("  under 100 bp: %d of %d" % (sum(1 for x in gaps if x < 100), n))
say("")
say("co-orientation, pairs with strand data: %d" % strand_n)
say("  same strand: %d" % co)
say("  opposite strand: %d" % opp)
say("")
say("CAVEAT: different contigs in a MAG is assembly fragmentation, not evidence against linkage.")

open(OUT, "w").write("\n".join(L) + "\n")
print("")
print("WROTE:", OUT)
