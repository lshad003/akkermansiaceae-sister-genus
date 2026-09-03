#!/usr/bin/env python3
# Is "66 wild animals" right? Count by animal unit, not by genome.
import os, sys
from collections import Counter, defaultdict

ROOT = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
NOV  = ROOT + "/results/novel_akk_tree/novel_size_gc.tsv"
EHI  = ROOT + "/data/ehi_2025_annotated/ehi_mags_annotated_v2.tsv"
CAT  = ROOT + "/data/unified_annotation/unified_genome_annotation.tsv"
OUT  = ROOT + "/results/novel_akk_tree/animal_count_verified.txt"

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

with open(NOV) as f:
    f.readline()
    n105 = [l.split("\t")[0].strip() for l in f if l.strip()]
say("genomes in the candidate set: %d" % len(n105))

ehi = [g for g in n105 if g.startswith("EHM")]
herp = [g for g in n105 if not g.startswith("EHM")]
say("  EHI: %d   herptile: %d" % (len(ehi), len(herp)))

# EHI animal = ncbi_biosample
bios, spec = {}, {}
with open(EHI) as f:
    h = f.readline().rstrip("\n").split("\t")
    ia = h.index("accession")
    ib = h.index("ncbi_biosample") if "ncbi_biosample" in h else None
    isrc = h.index("ncbi_isolation_source") if "ncbi_isolation_source" in h else None
    for line in f:
        p = line.rstrip("\n").split("\t")
        if ia < len(p) and p[ia].strip() in set(ehi):
            bios[p[ia].strip()] = p[ib].strip() if ib is not None and ib < len(p) else ""
            spec[p[ia].strip()] = p[isrc].strip() if isrc is not None and isrc < len(p) else ""

say("")
say("EHI: %d genomes -> %d distinct biosamples" % (len(ehi), len(set(bios.values()))))
dup = [k for k, v in Counter(bios.values()).items() if v > 1]
say("biosamples yielding more than one genome: %d" % len(dup))
for b in sorted(dup):
    gs = [g for g in ehi if bios.get(g) == b]
    say("   %-16s %s" % (b, ", ".join(gs)))
say("")
for k, v in Counter(spec.values()).most_common():
    say("   %-36s %3d genomes" % (k[:36], v))

# herptile animal = token before first period
ha = defaultdict(list)
for g in herp:
    ha[g.split(".")[0]].append(g)
say("")
say("herptile: %d genomes -> %d distinct animals" % (len(herp), len(ha)))
multi = {k: v for k, v in ha.items() if len(v) > 1}
say("animals yielding more than one genome: %d" % len(multi))
for k in sorted(multi):
    say("   %-14s %d: %s" % (k, len(multi[k]), ", ".join(multi[k])))

total = len(set(bios.values())) + len(ha)
say("")
say("=" * 60)
say("ANIMAL TOTAL: %d EHI + %d herptile = %d" % (len(set(bios.values())), len(ha), total))
say("manuscript currently says 66")
say("=" * 60)
if total != 66:
    say("MISMATCH. The manuscript number must change to %d, or the animal" % total)
    say("unit definition must be stated explicitly in Methods.")
else:
    say("MATCHES.")
say(""); say("written to " + OUT)
fh.close()
