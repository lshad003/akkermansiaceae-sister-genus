#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Stage 2: type-genome (CAND_C288) rows from delim35 and delim58 matrices, plus TableS1 tRNA column
import statistics
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
R = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus"
T = "UHM1073.23039_R.bin.125"
ROW = "CAND_C288"

def load(path):
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")[1:]
        rows = {}
        for line in f:
            p = line.rstrip("\n").split("\t")
            rows[p[0]] = dict(zip(hdr, [float(x) for x in p[1:]]))
    return hdr, rows

def rng(vals, tag):
    vals = sorted(vals)
    print("  %-38s n=%2d  min %.2f  median %.2f  max %.2f" % (tag, len(vals), vals[0], statistics.median(vals), vals[-1]))

for name in ("delim35", "delim58"):
    for metric in ("aai", "pocp"):
        path = "%s/results/%s/%s_matrix.tsv" % (B, name, metric)
        hdr, rows = load(path)
        if ROW not in rows:
            print("REFUSE: %s not in %s" % (ROW, path)); continue
        r = rows[ROW]
        cand = [r[c] for c in hdr if c.startswith("CAND_") and c != ROW]
        akk  = [r[c] for c in hdr if c.startswith("AKK_") or c == "Akkermansia"]
        free = [r[c] for c in hdr if not c.startswith("CAND_") and not c.startswith("AKK_") and c != "Akkermansia"]
        print("%s %s  (%s)" % (name, metric.upper(), path))
        rng(cand, "type vs other 16 candidate clusters")
        rng(akk,  "type vs Akkermansia reps")
        rng(free, "type vs free-living genus reps")
        if name == "delim58":
            print("  free-living columns: " + " ".join(c for c in hdr if not c.startswith("CAND_") and not c.startswith("AKK_")))

print()
p = "%s/results/delim58/TableS4_genus_boundary_metrics.tsv" % B
print("TableS4 (%s):" % p)
with open(p) as f:
    for line in f: print("  " + line.rstrip("\n")[:160])

print()
p = "%s/tables/TableS1_genome_quality.tsv" % R
with open(p) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    idx = [i for i, h in enumerate(hdr) if "trna" in h.lower()]
    print("TableS1 tRNA columns: " + ", ".join(hdr[i] for i in idx))
    for line in f:
        if line.startswith(T):
            p2 = line.rstrip("\n").split("\t")
            print("  %s -> " % T + ", ".join("%s=%s" % (hdr[i], p2[i]) for i in idx))
