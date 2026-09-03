#!/usr/bin/env python3
# v2: rel is computed in TRANSCRIPTION direction, not contig coordinate order.
# On a minus-strand operon the downstream gene has a LOWER index, so rel must flip.
import os, glob
from collections import defaultdict

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TSV = R + "/results/ppp_unified/operon_allgenera.tsv"
W   = R + "/results/operon_window"
FLANK = 4
os.makedirs(W, exist_ok=True)

faa = {}
for p in glob.glob(R + "/results/**/*.faa", recursive=True):
    faa.setdefault(os.path.basename(p)[:-4], p)

rows = []
with open(TSV) as f:
    h = f.readline().rstrip("\n").split("\t")
    for line in f:
        rows.append(dict(zip(h, line.rstrip("\n").split("\t"))))

byg = defaultdict(list)
for r in rows:
    if r["same_contig"] == "True" and r["order"] == "zwf-sub":
        byg[r["group"]].append(r)

def genome_index(g):
    recs, name = {}, None
    for line in open(faa[g]):
        if line.startswith(">"):
            p = line[1:].split(" # ")
            name = p[0].strip()
            if len(p) >= 4:
                recs[name] = dict(contig=name.rsplit("_", 1)[0],
                                  idx=int(name.rsplit("_", 1)[1]),
                                  start=int(p[1]), end=int(p[2]),
                                  strand="+" if p[3].strip() == "1" else "-", seq=[])
            else: name = None
        elif name and name in recs:
            recs[name]["seq"].append(line.strip())
    return recs

pick = {}
for grp, v in byg.items():
    best = None
    for r in v:
        if r["genome"] not in faa: continue
        recs = genome_index(r["genome"])
        z = next((k for k, d in recs.items()
                  if d["contig"] == r["ctg_zwf"] and d["start"] == int(r["zwf_start"])), None)
        if z is None: continue
        onctg = sorted([d for d in recs.values() if d["contig"] == recs[z]["contig"]],
                       key=lambda d: d["idx"])
        i = next(j for j, d in enumerate(onctg) if d["start"] == recs[z]["start"])
        room = min(i, len(onctg) - 1 - i)
        if best is None or room > best[0]:
            best = (room, r["genome"], recs, recs[z]["contig"], recs[z]["start"], recs[z]["strand"])
    if best: pick[grp] = best

order = ["CANDIDATE"] + sorted([k for k in pick if k != "CANDIDATE"],
                               key=lambda k: -len(byg[k]))
order = [k for k in order if k in pick]

fa = open(W + "/window.faa", "w")
tb = open(W + "/window_coords.tsv", "w")
tb.write("group\tgenome\tprotein\tcontig\tidx\tstart\tend\tstrand\tzwf_strand\trel\tlen_aa\trole\n")
print("  %-18s %-26s %-8s %s" % ("group", "representative", "zwf str", "window"))
for grp in order:
    room, g, recs, ctg, zstart, zstrand = pick[grp]
    onctg = sorted([(k, d) for k, d in recs.items() if d["contig"] == ctg],
                   key=lambda x: x[1]["idx"])
    i = next(j for j, (k, d) in enumerate(onctg) if d["start"] == zstart)
    lo, hi = max(0, i - FLANK - 1), min(len(onctg), i + FLANK + 2)
    win = onctg[lo:hi]
    print("  %-18s %-26s %-8s %d" % (grp[:18], g[:26], zstrand, len(win)))
    for k, d in win:
        raw = d["idx"] - onctg[i][1]["idx"]
        rel = raw if zstrand == "+" else -raw          # THE FIX
        role = "zwf" if rel == 0 else ("sub" if rel == 1 else "flank")
        fa.write(">%s|%s\n%s\n" % (grp, k, "".join(d["seq"])))
        tb.write("\t".join(map(str, [grp, g, k, d["contig"], d["idx"], d["start"], d["end"],
                                     d["strand"], zstrand, rel,
                                     len("".join(d["seq"])), role])) + "\n")
fa.close(); tb.close()
print(); print("wrote " + W + "/window.faa and window_coords.tsv")
