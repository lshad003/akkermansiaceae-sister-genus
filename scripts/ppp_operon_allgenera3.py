#!/usr/bin/env python3
# v3: coordinates from prodigal FASTA headers, so all genomes with both genes are used.
import os, glob
from collections import defaultdict

R   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
HIT = R + "/results/ppp_unified/ppp_trio.tsv"
CEN = R + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
NOV = R + "/results/novel_akk_tree/novel_size_gc.tsv"
OUT = R + "/results/ppp_unified/operon_allgenera.tsv"
TXT = R + "/results/ppp_unified/operon_allgenera.txt"

fh = open(TXT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

QN = {"UHM979.41089_R.bin.103_CDS_0654": "zwf", "EHM058980_CDS_1432": "sub"}

best = defaultdict(dict)
for line in open(HIT):
    f = line.rstrip("\n").split("\t")
    if len(f) < 5: continue
    q = QN.get(f[0])
    if q is None: continue
    g, prot = f[1].split("|", 1)
    pid = float(f[2])
    if q not in best[g] or pid > best[g][q][1]:
        best[g][q] = (prot, pid)

faa = {}
for p in glob.glob(R + "/results/**/*.faa", recursive=True):
    faa.setdefault(os.path.basename(p)[:-4], p)
say("proteomes indexed: %d" % len(faa))

fam, gen = {}, {}
with open(CEN) as f:
    h = f.readline().rstrip("\n").split("\t")
    ifm, ign = h.index("family"), h.index("genus")
    s = lambda v: v[3:] if v[:3] in ("f__", "g__") else v
    for line in f:
        p = line.rstrip("\n").split("\t")
        if ifm < len(p): fam[p[0].strip()] = s(p[ifm].strip())
        if ign < len(p): gen[p[0].strip()] = s(p[ign].strip())
with open(NOV) as f:
    f.readline(); cand = set(l.split("\t")[0].strip() for l in f if l.strip())

def coords(path, want):
    """prodigal header: >contig_N # start # end # strand # ID=..."""
    out = {}
    for line in open(path):
        if not line.startswith(">"): continue
        p = line[1:].split(" # ")
        if len(p) < 4: continue
        name = p[0].strip()
        if name not in want: continue
        contig = name.rsplit("_", 1)[0]
        try:
            out[name] = (contig, int(p[1]), int(p[2]), "+" if p[3].strip() == "1" else "-")
        except ValueError:
            pass
    return out

rows = []; skip = 0
for g, hits in best.items():
    if "zwf" not in hits or "sub" not in hits: continue
    if g not in faa: skip += 1; continue
    grp = ("CANDIDATE" if g in cand else
           gen.get(g) if fam.get(g) == "Akkermansiaceae" else None)
    if grp is None: continue
    c = coords(faa[g], {hits["zwf"][0], hits["sub"][0]})
    z, s2 = c.get(hits["zwf"][0]), c.get(hits["sub"][0])
    if not z or not s2: skip += 1; continue
    same = z[0] == s2[0]
    gap = None; order = ""
    if same:
        first, second = (z, s2) if z[1] < s2[1] else (s2, z)
        gap = second[1] - first[2] - 1
        up_zwf = (first is z)
        if z[3] == "-": up_zwf = not up_zwf
        order = "zwf-sub" if up_zwf else "sub-zwf"
    rows.append((grp, g, same, gap, z[3], s2[3], z[3] == s2[3], order,
                 z[0], z[1], z[2], s2[0], s2[1], s2[2]))

say("pairs resolved: %d   skipped: %d" % (len(rows), skip))
with open(OUT, "w") as f:
    f.write("group\tgenome\tsame_contig\tgap_bp\tstrand_zwf\tstrand_sub\tco_oriented\torder\t"
            "ctg_zwf\tzwf_start\tzwf_end\tctg_sub\tsub_start\tsub_end\n")
    for r in rows: f.write("\t".join(str(x) for x in r) + "\n")

byg = defaultdict(list)
for r in rows: byg[r[0]].append(r)
say(""); say("  %-18s %5s %9s %10s %11s %s" % ("group","n","same-ctg","co-orient","median gap","order"))
for grp in sorted(byg, key=lambda k: (k != "CANDIDATE", -len(byg[k]))):
    v = byg[grp]; sc = [x for x in v if x[2]]; co = [x for x in sc if x[6]]
    gaps = sorted(x[3] for x in sc if x[3] is not None)
    od = defaultdict(int)
    for x in sc: od[x[7]] += 1
    say("  %-18s %5d %9s %10s %11s %s"
        % (grp[:18], len(v), "%d/%d" % (len(sc), len(v)),
           "%d/%d" % (len(co), len(sc)) if sc else "-",
           "%d bp" % gaps[len(gaps)//2] if gaps else "-",
           ", ".join("%s %d" % (k, w) for k, w in sorted(od.items(), key=lambda x: -x[1]))))
say(""); say("written to " + OUT)
fh.close()
