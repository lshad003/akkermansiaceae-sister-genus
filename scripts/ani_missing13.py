#!/usr/bin/env python3
# Are the 13 catalogue-missing bins redundant with the 105, or extra genomes?
import os, sys, glob, subprocess, shutil
ROOT = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W    = ROOT + "/results/novel_akk_tree"
QL   = W + "/missing_bin_paths_for_ani.txt"
NOV  = W + "/novel_size_gc.tsv"
FDIR = ROOT + "/data/amphibia_gtdbtk_input"
OUT  = W + "/missing13_ani.txt"

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

SK = shutil.which("skani")
if not SK:
    for d in glob.glob("/bigdata/stajichlab/lshad003/condaenvs/*/bin/skani") + \
             glob.glob("/bigdata/stajichlab/shared/condaenvs/*/bin/skani"):
        SK = d; break
say("skani: %s" % SK)
if not SK:
    say("skani not found. paste: ls /bigdata/stajichlab/lshad003/condaenvs/"); fh.close(); sys.exit(1)
if not os.path.exists(QL):
    say("MISSING: " + QL); fh.close(); sys.exit(1)

q = [l.strip() for l in open(QL) if l.strip() and os.path.exists(l.strip())]
say("query FASTAs found: %d" % len(q))
with open(NOV) as f:
    f.readline(); n105 = [l.split("\t")[0].strip() for l in f if l.strip()]
EXT = (".fa", ".fna", ".fasta")
fmap = {}
for fn in os.listdir(FDIR):
    b = fn
    for e in EXT:
        if b.endswith(e): b = b[:-len(e)]; break
    fmap.setdefault(b, os.path.join(FDIR, fn))
refs = [fmap[g] for g in n105 if g in fmap]
say("reference genomes: %d" % len(refs))

qf, rf, dist = W + "/_q.txt", W + "/_r.txt", W + "/missing13_skani.tsv"
open(qf, "w").write("\n".join(q) + "\n")
open(rf, "w").write("\n".join(refs) + "\n")
r = subprocess.run([SK, "dist", "--ql", qf, "--rl", rf, "-o", dist, "-t", "8", "--min-af", "15"],
                   capture_output=True, text=True)
if r.returncode != 0:
    say("skani failed:"); say(r.stderr[:1200]); fh.close(); sys.exit(1)

best = {}
with open(dist) as f:
    h = f.readline().rstrip("\n").split("\t")
    ia, ir, iq = h.index("ANI"), h.index("Ref_file"), h.index("Query_file")
    for line in f:
        p = line.rstrip("\n").split("\t")
        try: a = float(p[ia])
        except (ValueError, IndexError): continue
        qn = os.path.basename(p[iq])
        if qn not in best or a > best[qn][0]:
            best[qn] = (a, os.path.basename(p[ir]))

say(""); say("  %-40s %-8s %s" % ("missing bin", "best ANI", "closest of the 105"))
dupe = extra = novel = nohit = 0
for p in q:
    b = os.path.basename(p)
    v = best.get(b)
    if not v:
        say("  %-40s %-8s %s" % (b, "no hit", "(below floor)")); nohit += 1; continue
    a, ref = v
    say("  %-40s %-8.2f %s" % (b, a, ref))
    if a >= 99: dupe += 1
    elif a >= 95: extra += 1
    else: novel += 1

say("")
say("=" * 66)
say("  >=99 duplicate, nothing lost          : %d" % dupe)
say("  95-99 same species, an extra genome   : %d" % extra)
say("  <95  distinct species of the genus    : %d" % novel)
say("  no hit (expect the 2 tortoise bins)   : %d" % nohit)
say("=" * 66)
say("If all are >=95, the counts stand and one sentence covers it.")
say("If any are <95, the clade holds unsampled species diversity.")
for t in (qf, rf):
    if os.path.exists(t): os.remove(t)
say(""); say("written to " + OUT)
fh.close()
