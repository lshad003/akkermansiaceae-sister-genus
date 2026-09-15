#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Does partner completeness or proteome size drive POCP to the type genome (CAND_C288)?
import sys, statistics
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
S1 = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus/tables/TableS1_genome_quality.tsv"
META = B + "/data/gtdb_r226/bac120_metadata_r226.tsv"
OUT = B + "/results/delim58/pocp_completeness_check.tsv"
TYPE = "CAND_C288"

def matrix(p):
    with open(p) as f:
        hdr = f.readline().rstrip("\n").split("\t")[1:]
        for line in f:
            x = line.rstrip("\n").split("\t")
            if x[0] == TYPE: return dict(zip(hdr, map(float, x[1:])))
    sys.exit("REFUSED: %s row missing in %s" % (TYPE, p))

pocp = matrix(B + "/results/delim58/pocp_matrix.tsv")
aai  = matrix(B + "/results/delim58/aai_matrix.tsv")

members = {}
with open(B + "/results/delim58/members.tsv") as f:
    f.readline()
    for line in f:
        lab, gen, path = line.rstrip("\n").split("\t")[:3]
        members[lab] = (gen, path)

comp, cont = {}, {}
with open(S1) as f:
    f.readline()
    for line in f:
        x = line.rstrip("\n").split("\t")
        comp[x[0]] = float(x[2]); cont[x[0]] = float(x[3])

need = {g for l, (g, p) in members.items() if not l.startswith("CAND_")}
alt = {g.replace("RS_", "").replace("GB_", ""): g for g in need}
with open(META) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ci = [i for i, h in enumerate(hdr) if "completeness" in h.lower()]
    ki = [i for i, h in enumerate(hdr) if "contamination" in h.lower()]
    cc = [i for i in ci if "checkm2" in hdr[i]] or ci
    kk = [i for i in ki if "checkm2" in hdr[i]] or ki
    print("metadata completeness column:", hdr[cc[0]])
    for line in f:
        acc = line.split("\t", 1)[0]
        key = acc if acc in need else alt.get(acc.replace("RS_", "").replace("GB_", ""))
        if key:
            x = line.rstrip("\n").split("\t")
            comp[key] = float(x[cc[0]]); cont[key] = float(x[kk[0]])

def nprot(p):
    n = 0
    with open(p) as f:
        for line in f:
            if line.startswith(">"): n += 1
    return n

rows = []
for lab, (gen, path) in members.items():
    if lab == TYPE: continue
    grp = "candidate" if lab.startswith("CAND_") else ("akkermansia" if lab.startswith("AKK_") else "free")
    rows.append((grp, lab, gen, comp.get(gen), cont.get(gen), nprot(path), pocp[lab], aai[lab]))

def rank(v):
    s = sorted(range(len(v)), key=lambda i: v[i]); r = [0]*len(v)
    for k, i in enumerate(s): r[i] = k
    return r
def spearman(x, y):
    rx, ry = rank(x), rank(y); n = len(x)
    mx, my = sum(rx)/n, sum(ry)/n
    num = sum((a-mx)*(b-my) for a, b in zip(rx, ry))
    den = (sum((a-mx)**2 for a in rx) * sum((b-my)**2 for b in ry)) ** 0.5
    return num/den if den else float("nan")

with open(OUT, "w") as o:
    o.write("group\tlabel\tgenome\tcompleteness\tcontamination\tn_proteins\tPOCP_to_type\tAAI_to_type\n")
    for r in sorted(rows):
        o.write("\t".join("NA" if v is None else str(v) for v in r) + "\n")

print("\n%-12s %-14s %-24s %6s %6s %7s %7s %7s" % ("group", "label", "genome", "compl", "cont", "nprot", "POCP", "AAI"))
for r in sorted(rows, key=lambda r: (r[0], -r[6])):
    print("%-12s %-14s %-24s %6s %6s %7d %7.2f %7.2f" % (r[0], r[1], r[2][:24],
          "NA" if r[3] is None else "%.1f" % r[3], "NA" if r[4] is None else "%.1f" % r[4], r[5], r[6], r[7]))

print("\nSpearman rho of POCP-to-type within group:")
for grp in ("candidate", "akkermansia", "free"):
    g = [r for r in rows if r[0] == grp and r[3] is not None]
    if len(g) < 4: print("  %-12s n=%d, too few with completeness" % (grp, len(g))); continue
    print("  %-12s n=%2d  vs completeness %+.2f   vs n_proteins %+.2f   completeness range %.1f to %.1f" % (
        grp, len(g), spearman([r[3] for r in g], [r[6] for r in g]),
        spearman([r[5] for r in g], [r[6] for r in g]), min(r[3] for r in g), max(r[3] for r in g)))
low = [r for r in rows if r[3] is not None and r[3] < 70]
print("\nreps below 70 percent completeness: %d" % len(low))
for r in low: print("  %s %s %.1f  POCP %.2f" % (r[0], r[2], r[3], r[6]))
missing = [r[2] for r in rows if r[3] is None]
print("reps with no completeness found: %d %s" % (len(missing), missing))
print("\nwrote", OUT)
