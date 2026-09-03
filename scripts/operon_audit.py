#!/usr/bin/env python3
# Three audits: literal adjacency, the 91/95 vs 93/95 discrepancy,
# and oriented distance-to-contig-edge for the split pairs.
import os, glob
from collections import defaultdict

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TSV = R + "/results/ppp_unified/operon_allgenera.tsv"
OLD = R + "/results/pangenome/ppp_operon_verified.txt"
OUT = R + "/results/ppp_unified/operon_audit.txt"

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

faa = {}
for p in glob.glob(R + "/results/**/*.faa", recursive=True):
    faa.setdefault(os.path.basename(p)[:-4], p)

rows = []
with open(TSV) as f:
    h = f.readline().rstrip("\n").split("\t")
    for line in f: rows.append(dict(zip(h, line.rstrip("\n").split("\t"))))

def index(g):
    d = {}
    for line in open(faa[g]):
        if line.startswith(">"):
            p = line[1:].split(" # ")
            if len(p) >= 4:
                n = p[0].strip()
                d[(n.rsplit("_",1)[0], int(p[1]))] = dict(
                    idx=int(n.rsplit("_",1)[1]), start=int(p[1]), end=int(p[2]),
                    strand="+" if p[3].strip()=="1" else "-", contig=n.rsplit("_",1)[0])
    return d

# ---------- AUDIT 1: literal adjacency ----------
say("=" * 70); say("AUDIT 1. ZERO INTERVENING CDS?"); say("=" * 70)
byg = defaultdict(lambda: defaultdict(int))
detail = defaultdict(list)
for r in rows:
    if r["same_contig"] != "True" or r["genome"] not in faa: continue
    idx = index(r["genome"])
    z = idx.get((r["ctg_zwf"], int(r["zwf_start"])))
    s = idx.get((r["ctg_sub"], int(r["sub_start"])))
    if not z or not s: byg[r["group"]]["unresolved"] += 1; continue
    iv = abs(z["idx"] - s["idx"]) - 1
    byg[r["group"]]["adjacent" if iv == 0 else "gapped"] += 1
    if iv: detail[r["group"]].append((r["genome"], iv, r["gap_bp"]))
say("  %-18s %9s %8s %11s" % ("group", "adjacent", "gapped", "unresolved"))
ta = tg = 0
for g in sorted(byg, key=lambda k: -sum(byg[k].values())):
    a, gp, u = byg[g]["adjacent"], byg[g]["gapped"], byg[g]["unresolved"]
    ta += a; tg += gp
    say("  %-18s %9d %8d %11d" % (g[:18], a, gp, u))
say("  %-18s %9d %8d" % ("TOTAL", ta, tg))
for g, v in detail.items():
    if v:
        say(""); say("  gapped in %s:" % g)
        for gn, iv, bp in v[:6]:
            say("     %-30s %d intervening CDS, %s bp" % (gn[:30], iv, bp))

# ---------- AUDIT 2: 91/95 vs 93/95 ----------
say(""); say("=" * 70); say("AUDIT 2. ARCHIVED 91/95 VS NEW 93/95"); say("=" * 70)
if os.path.exists(OLD):
    for line in open(OLD): say("  OLD | " + line.rstrip())
c = [r for r in rows if r["group"] == "CANDIDATE"]
say(""); say("  NEW | candidate genomes with both genes: %d" % len(c))
say("  NEW | same contig: %d" % sum(1 for r in c if r["same_contig"] == "True"))
say("  NEW | split: %s" % ", ".join(r["genome"] for r in c if r["same_contig"] != "True"))
say("")
say("  The old run used ppp_all199.tsv with a 30 percent identity floor.")
say("  The new run used the unified search with no identity filter, so it can")
say("  recover a pair the old one missed. That is the likely cause. Confirm by")
say("  checking whether the two extra genomes have a sub hit below 30 percent.")

# ---------- AUDIT 3: contig edges for split pairs ----------
say(""); say("=" * 70); say("AUDIT 3. ORIENTED DISTANCE TO CONTIG EDGE, SPLIT PAIRS"); say("=" * 70)
def clen(g):
    d = {}
    for line in open(faa[g]):
        if line.startswith(">"):
            p = line[1:].split(" # ")
            if len(p) >= 4:
                c = p[0].strip().rsplit("_",1)[0]
                d[c] = max(d.get(c, 0), int(p[2]))
    return d
say("  %-16s %-28s %10s %10s" % ("group", "genome", "zwf->edge", "sub->edge"))
n = 0
for r in rows:
    if r["same_contig"] == "True" or r["genome"] not in faa: continue
    L = clen(r["genome"])
    lz, ls = L.get(r["ctg_zwf"], 0), L.get(r["ctg_sub"], 0)
    dz = min(int(r["zwf_start"]), lz - int(r["zwf_end"]))
    ds = min(int(r["sub_start"]), ls - int(r["sub_end"]))
    say("  %-16s %-28s %10d %10d" % (r["group"][:16], r["genome"][:28], dz, ds))
    n += 1
    if n >= 40: say("  ... more"); break
say("")
say("  Small values mean the genes abut contig ends, i.e. fragmentation.")
say("  Large values mean the pair is genuinely separated.")
say(""); say("written to " + OUT)
fh.close()
