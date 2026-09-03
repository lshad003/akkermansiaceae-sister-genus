#!/usr/bin/env python3
# Prodigal + DIAMOND reciprocal best hits: 135 Wilkie genomes vs the type genome.
import os, sys, glob, subprocess, shutil
from statistics import mean

ROOT = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W    = ROOT + "/results/wilkie_aai"
GEN  = W + "/genomes"
FAA  = W + "/faa"
OUT  = W + "/aai_vs_type.tsv"
LOG  = W + "/run.log"
TYPE = "EHM058340"

os.makedirs(FAA, exist_ok=True)
fh = open(LOG, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

# ---- tools
def find(name):
    p = shutil.which(name)
    if p: return p
    for d in glob.glob("/bigdata/stajichlab/lshad003/condaenvs/*/bin/" + name):
        return d
    for d in glob.glob("/bigdata/stajichlab/shared/condaenvs/*/bin/" + name):
        return d
    return None

PROD, DIAM = find("prodigal"), find("diamond")
say("prodigal: %s" % PROD)
say("diamond : %s" % DIAM)
if not PROD or not DIAM:
    say(""); say("MISSING A TOOL. Run this and paste the output:")
    say("  ls /bigdata/stajichlab/lshad003/condaenvs/")
    say("  module avail 2>&1 | grep -i -E 'prodigal|diamond'")
    fh.close(); sys.exit(1)

# ---- type genome proteome
cand = []
for pat in ("/results/pangenome/gff199/*%s*.faa", "/results/pangenome/*%s*.faa",
            "/data/**/%s*.faa"):
    cand += glob.glob(ROOT + pat % TYPE, recursive=True)
if not cand:
    say(""); say("type proteome for %s not found. paste:" % TYPE)
    say("  ls %s/results/pangenome/gff199 | head" % ROOT)
    fh.close(); sys.exit(1)
TFAA = cand[0]
say("type proteome: %s" % TFAA)

# ---- prodigal
fnas = sorted(glob.glob(GEN + "/*.fna"))
say(""); say("genomes: %d" % len(fnas))
for i, f in enumerate(fnas, 1):
    b = os.path.basename(f)[:-4]
    o = "%s/%s.faa" % (FAA, b)
    if os.path.exists(o) and os.path.getsize(o) > 0:
        continue
    subprocess.run([PROD, "-i", f, "-a", o, "-p", "meta", "-q"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if i % 25 == 0:
        say("  prodigal %d/%d" % (i, len(fnas)))
say("proteomes: %d" % len(glob.glob(FAA + "/*.faa")))

# ---- diamond
TDB = W + "/_type"
subprocess.run([DIAM, "makedb", "--in", TFAA, "-d", TDB, "--quiet"])

def best(q, d):
    r = subprocess.run([DIAM, "blastp", "-q", q, "-d", d, "--quiet", "--max-target-seqs", "1",
                        "-e", "1e-5", "--id", "30", "--query-cover", "50", "-p", "8",
                        "-f", "6", "qseqid", "sseqid", "pident"],
                       capture_output=True, text=True)
    h = {}
    for line in r.stdout.splitlines():
        p = line.split("\t")
        if len(p) == 3 and p[0] not in h:
            h[p[0]] = (p[1], float(p[2]))
    return h

res = []
faas = sorted(glob.glob(FAA + "/*.faa"))
say(""); say("running AAI on %d proteomes" % len(faas))
for i, q in enumerate(faas, 1):
    b = os.path.basename(q)[:-4]
    qdb = W + "/_q"
    subprocess.run([DIAM, "makedb", "--in", q, "-d", qdb, "--quiet"])
    fwd, rev = best(q, TDB + ".dmnd"), best(TFAA, qdb + ".dmnd")
    rbh = [v[1] for k, v in fwd.items() if rev.get(v[0], ("", 0))[0] == k]
    res.append((b, len(rbh), mean(rbh) if rbh else 0.0))
    if i % 20 == 0:
        say("  %d/%d" % (i, len(faas)))

res.sort(key=lambda x: -x[2])
with open(OUT, "w") as f:
    f.write("genome\tn_RBH\tAAI_pct\n")
    for b, n, a in res:
        f.write("%s\t%d\t%.2f\n" % (b, n, a))

say(""); say("=" * 62); say("TOP 15 BY AAI TO THE TYPE GENOME"); say("=" * 62)
for b, n, a in res[:15]:
    say("  %-22s %5d RBH   %6.2f%%" % (b, n, a))
vals = [a for _, _, a in res if a > 0]
say("")
say("n = %d   min %.2f   max %.2f" % (len(vals), min(vals), max(vals)))
say("")
say("READ IT: within-genus is 73.7 and above. Akkermansia is 56.6.")
say("  Anything above ~70 is our genus in a non-amphibian host and changes the paper.")
say("  All clustered near 56-57 means no Wilkie genome belongs to our genus.")
for t in (TDB + ".dmnd", W + "/_q.dmnd"):
    if os.path.exists(t): os.remove(t)
say(""); say("written to " + OUT)
fh.close()
