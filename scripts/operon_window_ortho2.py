#!/usr/bin/env python3
# v2: reciprocal best hits instead of single linkage, so zwf and sub cannot merge
# through a shared Rossmann fold. One RBH per genus pair per protein.
import os, subprocess
from collections import defaultdict

R = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = R + "/results/operon_window"
D = "/bigdata/stajichlab/lshad003/condaenvs/MYNEWENV/bin/diamond"
OUT = W + "/window_groups.tsv"
PID, COV = 30.0, 50.0

subprocess.run([D, "makedb", "--in", W + "/window.faa", "-d", W + "/_w", "--quiet"], check=True)
r = subprocess.run([D, "blastp", "-q", W + "/window.faa", "-d", W + "/_w", "--quiet",
                    "--very-sensitive", "-e", "1e-3", "--max-target-seqs", "500",
                    "-p", "4", "-f", "6", "qseqid", "sseqid", "pident", "qcovhsp", "bitscore"],
                   capture_output=True, text=True)
print("alignments:", len(r.stdout.splitlines()))

# best hit from each protein into each OTHER genus, by bitscore
best = defaultdict(dict)
score = {}
for line in r.stdout.splitlines():
    q, s, p, c, b = line.split("\t")
    if q == s: continue
    if float(p) < PID or float(c) < COV: continue
    gq, gs = q.split("|")[0], s.split("|")[0]
    if gq == gs: continue
    b = float(b)
    if gs not in best[q] or b > best[q][gs][1]:
        best[q][gs] = (s, b)
    score[(q, s)] = max(score.get((q, s), 0), float(p))

# reciprocal pairs only
par = {}
def find(x):
    par.setdefault(x, x)
    while par[x] != x:
        par[x] = par[par[x]]; x = par[x]
    return x
def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb: par[ra] = rb

rbh = set()
for q, d in best.items():
    for gs, (s, _) in d.items():
        gq = q.split("|")[0]
        back = best.get(s, {}).get(gq)
        if back and back[0] == q:
            rbh.add(tuple(sorted([q, s])))
            union(q, s)
print("reciprocal best-hit pairs:", len(rbh))

rows = [l.rstrip("\n").split("\t") for l in open(W + "/window_coords.tsv")]
hdr = rows[0]
for r2 in rows[1:]:
    d = dict(zip(hdr, r2)); find("%s|%s" % (d["group"], d["protein"]))

grp = defaultdict(list)
for x in par: grp[find(x)].append(x)
sizes = sorted(grp.values(), key=len, reverse=True)

# sanity: does any group hold both zwf and sub?
role = {}
for r2 in rows[1:]:
    d = dict(zip(hdr, r2)); role["%s|%s" % (d["group"], d["protein"])] = d["role"]
bad = [g for g in sizes if len({role.get(m) for m in g} & {"zwf", "sub"}) == 2]
print("groups containing BOTH zwf and sub:", len(bad), "(must be 0)")

lab = {}
for i, g in enumerate(sizes):
    for m in g: lab[m] = "OG%02d" % (i+1)
with open(OUT, "w") as f:
    f.write("\t".join(hdr + ["og", "og_size", "og_genera"]) + "\n")
    for r2 in rows[1:]:
        d = dict(zip(hdr, r2)); tag = "%s|%s" % (d["group"], d["protein"])
        mem = grp.get(find(tag), [tag])
        f.write("\t".join(r2 + [lab.get(tag, "singleton"), str(len(mem)),
                                str(len({m.split('|')[0] for m in mem}))]) + "\n")
with open(W + "/window_pid.tsv", "w") as f:
    for (a, b) in rbh:
        f.write("%s\t%s\t%.1f\n" % (a, b, max(score.get((a,b),0), score.get((b,a),0))))

multi = [g for g in sizes if len({m.split("|")[0] for m in g}) > 1]
print("ortholog groups:", len(sizes), " spanning 2+ genera:", len(multi))
print()
print("  %-8s %5s %7s  %s" % ("group", "n", "genera", "roles"))
for i, g in enumerate(sizes[:16]):
    print("  %-8s %5d %7d  %s" % ("OG%02d" % (i+1), len(g),
          len({m.split('|')[0] for m in g}),
          ",".join(sorted({role.get(m, "?") for m in g}))))
if os.path.exists(W + "/_w.dmnd"): os.remove(W + "/_w.dmnd")
print(); print("wrote " + OUT)
