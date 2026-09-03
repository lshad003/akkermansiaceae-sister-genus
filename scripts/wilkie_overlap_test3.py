#!/usr/bin/env python3
# Fixed field split: match genus and source against known vocabularies, longest first.
import os, re, sys
from collections import Counter

ROOT = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TSV  = ROOT + "/data/wilkie444.tsv"
TREE = ROOT + "/results/novel_akk_tree/akk_placement.nwk"
OUT  = ROOT + "/results/novel_akk_tree/wilkie_overlap.txt"
TODO = ROOT + "/results/novel_akk_tree/wilkie_to_download.txt"

GENERA = sorted(["Akkermansia","Luteolibacter","SW10","Oceaniferula","Roseibacillus_B",
                 "Haloferula","Rubritalea","UBA956","Roseibacillus","WTJZ01","UBA4581",
                 "Persicirhabdus","JAPKDA01","CALFDI01","CALFBI01"], key=len, reverse=True)
SOURCES = sorted(["gut","marine","freshwater","other","estuarine","animal_(undefined)",
                  "animal_","sediment","engineered","rhizosphere","soil","sponge",
                  "wetland","brackish","human","hypersaline","groundwater","coral"],
                 key=len, reverse=True)

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

for p in (TSV, TREE):
    if not os.path.exists(p):
        say("MISSING: " + p); fh.close(); sys.exit(1)

rows, bad = [], []
with open(TSV) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("accession"):
            continue
        m = re.match(r'^(GC[AF]_\d{9}\.\d+)', line)
        acc = m.group(1) if m else ""
        i = line.find("g__")
        if i < 0:
            bad.append(line[:60]); continue
        rest = line[i+3:]
        gen = next((g for g in GENERA if rest.startswith(g)), None)
        if gen is None:
            bad.append(line[:60]); continue
        rest = rest[len(gen):]
        src = next((s for s in SOURCES if rest.startswith(s)), None)
        if src is None:
            bad.append(line[:60]); continue
        rows.append((acc, gen, src, rest[len(src):]))

say("parsed %d   unparseable %d" % (len(rows), len(bad)))
for b in bad[:5]:
    say("   BAD: " + b)
say("")
say("SANITY, expect gut 198 / marine 117 / freshwater 74:")
for k, v in Counter(r[2] for r in rows).most_common():
    say("   %-22s %4d" % (k, v))
say("")
say("SANITY, expect Akkermansia 206 / Luteolibacter 82:")
for k, v in Counter(r[1] for r in rows).most_common(6):
    say("   %-22s %4d" % (k, v))

tree = open(TREE).read()
base = set(a.split(".")[0] for a in re.findall(r'GC[AF]_\d{9}\.\d+', tree))
have  = [r for r in rows if r[0] and r[0].split(".")[0] in base]
miss  = [r for r in rows if r[0] and r[0].split(".")[0] not in base]
noacc = [r for r in rows if not r[0]]

say(""); say("IN OUR TREE %d   NOT IN TREE %d   NO ACCESSION %d"
             % (len(have), len(miss), len(noacc)))

host = [r for r in miss if r[2] in ("gut", "human", "animal_", "animal_(undefined)")]
say("")
say("=" * 72)
say("HOST-ASSOCIATED AND NOT IN OUR TREE: %d" % len(host))
say("Every one of these must be AAI-tested against the type genome.")
say("=" * 72)
for k, v in Counter("%s / %s / %s" % (r[1], r[2], r[3] or "none") for r in host).most_common():
    say("   %-52s %4d" % (k[:52], v))

with open(TODO, "w") as f:
    for a, g, s, o in host:
        f.write("%s\t%s\t%s\t%s\n" % (a, g, s, o))
say(""); say("%d accessions written to %s" % (len(host), TODO))
say(""); say("written to " + OUT)
fh.close()
