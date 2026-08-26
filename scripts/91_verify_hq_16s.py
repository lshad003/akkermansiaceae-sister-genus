#!/usr/bin/env python3
# Ribosomal RNA of the high-quality genomes checked for provenance
# Source: ch3-chitin-evolution/scripts/verify_hq_16s.py
# Output: results/rrna16s/hq_16s_verdict.txt
import os

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
MIM = CH3 + "results/mimag"
FA = CH3 + "data/amphibia_gtdbtk_input"
OUT = CH3 + "results/rrna16s/hq_16s_verdict.txt"

HQ = ["EHM058340", "EHM059280", "EHM059433", "EHM060011", "UHM973.23044_R.bin.53"]
CORE = ["EHM060011", "EHM060058", "UHM1327.23073_R.bin.116", "UHM973.23044_R.bin.53"]

L = []
def say(s):
    print(s)
    L.append(s)

def readfa(p):
    d = {}
    n = None
    b = []
    for line in open(p):
        if line.startswith(">"):
            if n:
                d[n] = "".join(b)
            n = line[1:].split()[0]
            b = []
        else:
            b.append(line.strip())
    if n:
        d[n] = "".join(b)
    return d

def rc(s):
    return s[::-1].translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))

def get16s(acc):
    g = MIM + "/" + acc + ".rrna.gff"
    if not os.path.exists(g):
        return []
    p = None
    for e in (".fa", ".fna", ".fasta"):
        if os.path.exists(FA + "/" + acc + e):
            p = FA + "/" + acc + e
            break
    if p is None:
        return []
    seqs = readfa(p)
    out = []
    for line in open(g):
        if line.startswith("#") or "\t" not in line:
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or "16S" not in f[8]:
            continue
        if f[0] not in seqs:
            continue
        s = seqs[f[0]][int(f[3]) - 1:int(f[4])]
        if f[6] == "-":
            s = rc(s)
        out.append((f[0], len(seqs[f[0]]), s.upper()))
    return out

K = 12
def km(s):
    return set(s[i:i + K] for i in range(len(s) - K + 1))

core = []
for a in CORE:
    for ctg, cl, s in get16s(a):
        if len(s) >= 1200 and cl >= 20000:
            core.append(km(s))
say("confirmed genus 16S used as the reference set: %d sequences" % len(core))
say("")
say("CONTAINMENT = fraction of the query's 12-mers found in a reference 16S.")
say("Length-robust, so partial sequences are scored fairly.")
say("")

for acc in HQ:
    say("=== %s" % acc)
    hits = get16s(acc)
    if not hits:
        say("    no 16S found")
        say("")
        continue
    for i, (ctg, cl, s) in enumerate(hits, 1):
        k = km(s)
        if not k or not core:
            continue
        best = max(len(k & c) / float(len(k)) for c in core)
        if best >= 0.35:
            v = "GENUS"
        elif best <= 0.10:
            v = "FOREIGN, misbinned"
        else:
            v = "AMBIGUOUS"
        part = "" if len(s) >= 1200 else "   PARTIAL"
        say("    copy %d  %6d bp on a %7d bp contig   containment %.3f   %s%s"
            % (i, len(s), cl, best, v, part))
    say("")

say("MIMAG asks whether the 16S gene is present. It does not ask whether the")
say("gene belongs to the genome, and it is not usually read as satisfied by a")
say("partial sequence. Both matter for the count of 5.")

open(OUT, "w").write("\n".join(L) + "\n")
print("")
print("WROTE:", OUT)
