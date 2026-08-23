#!/usr/bin/env python3
# Loss categories tested against the retained background
# Source: ch3-chitin-evolution/scripts/polarity_cog_enrichment.py
# Output: results/pangenome/polarity_cog_enrichment.tsv
import collections, math, os, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
EGG  = f"{BASE}/results/pangenome/eggnog_polarity"
OUT  = f"{BASE}/results/pangenome/polarity_cog_enrichment.tsv"

COG = {
 "J":"translation, ribosome","A":"RNA processing","K":"transcription","L":"replication and repair",
 "B":"chromatin","D":"cell cycle and division","V":"defence","T":"signal transduction",
 "M":"cell wall, membrane, envelope","N":"cell motility","Z":"cytoskeleton",
 "W":"extracellular structures","U":"trafficking and secretion",
 "O":"protein turnover and chaperones","X":"mobilome","C":"energy production",
 "G":"carbohydrate transport and metabolism","E":"amino acid transport and metabolism",
 "F":"nucleotide transport and metabolism","H":"coenzyme transport and metabolism",
 "I":"lipid transport and metabolism","P":"inorganic ion transport",
 "Q":"secondary metabolites","R":"general function only","S":"function unknown"}

def cogs(cat):
    p = f"{EGG}/{cat}.emapper.annotations"
    if not os.path.exists(p):
        return None
    hdr, c, n = None, collections.Counter(), 0
    for line in open(p):
        if line.startswith("#query"):
            hdr = line.lstrip("#").rstrip("\n").split("\t"); continue
        if line.startswith("#") or hdr is None:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < len(hdr):
            continue
        n += 1
        v = dict(zip(hdr, f)).get("COG_category", "-").strip()
        if v in ("", "-"):
            continue
        for ch in set(v):
            c[ch] += 1
    return n, c

def lfact(n):
    return math.lgamma(n + 1)

def fisher_two_sided(a, b, c, d):
    # a,b top row; c,d bottom row
    n = a + b + c + d
    def p(x):
        y = a + b - x
        z = a + c - x
        w = d - (a - x)
        if min(x, y, z, w) < 0:
            return 0.0
        return math.exp(lfact(a+b) + lfact(c+d) + lfact(a+c) + lfact(b+d)
                        - lfact(n) - lfact(x) - lfact(y) - lfact(z) - lfact(w))
    obs = p(a)
    tot = 0.0
    lo = max(0, a - d)
    hi = min(a + b, a + c)
    for x in range(lo, hi + 1):
        px = p(x)
        if px <= obs * (1 + 1e-9):
            tot += px
    return min(1.0, tot)

BG = cogs("shared_gut")
if BG is None:
    print("REFUSED: shared_gut annotations not found. Has job 27711022 finished?")
    sys.exit(1)
nBG, cBG = BG
print("background: shared_gut, %d annotated orthogroups" % nBG)

rows = []
for cat in ("B_lost_at_gut_ancestor", "A_lost_in_Akkermansia",
            "C_Akkermansia_enriched", "D_novel_specific"):
    r = cogs(cat)
    if r is None:
        print("skipping, not annotated:", cat); continue
    nS, cS = r
    print()
    print("=== %s vs retained background (n=%d annotated) ===" % (cat, nS))
    res = []
    for k in sorted(set(cS) | set(cBG)):
        a = cS.get(k, 0); b = nS - a
        c = cBG.get(k, 0); d = nBG - c
        if a + c < 5:
            continue
        pv = fisher_two_sided(a, b, c, d)
        fs = (a / nS) / (c / nBG) if c else float("inf")
        res.append((k, a, 100.0*a/nS, c, 100.0*c/nBG, fs, pv))
    res.sort(key=lambda x: x[6])
    m = len(res)
    print("%-4s %-36s %6s %7s %6s %7s %7s %9s %9s"
          % ("COG","category","set n","set %","bg n","bg %","ratio","p","q"))
    for i, (k, a, ap, c, cp, fs, pv) in enumerate(res, 1):
        q = min(1.0, pv * m / i)
        print("%-4s %-36s %6d %6.1f%% %6d %6.1f%% %7.2f %9.2g %9.2g"
              % (k, COG.get(k, "?"), a, ap, c, cp, fs, pv, q))
        rows.append([cat, k, COG.get(k, "?"), a, "%.1f" % ap, c, "%.1f" % cp,
                     "%.2f" % fs, "%.3g" % pv, "%.3g" % q])

with open(OUT, "w") as fh:
    fh.write("category\tCOG\tCOG_name\tset_n\tset_pct\tbackground_n\tbackground_pct\t"
             "ratio\tp_fisher\tq_BH\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")
print()
print("wrote", OUT)
print()
print("READ IT: the background is the 1,026 orthogroups RETAINED across both gut genera,")
print("not the genome as a whole. Retained genes are by definition those constrained enough")
print("to survive in both lineages, so they skew toward core functions. State results as")
print("'lost genes differ from retained genes in X', never as enrichment against the genome.")
