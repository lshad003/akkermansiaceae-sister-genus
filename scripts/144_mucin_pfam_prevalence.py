import csv
from collections import defaultdict

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/mucin_pfam"

mem = {r["genome"]: r["group"] for r in csv.DictReader(open(W + "/members.tsv"), delimiter="\t")}
N = defaultdict(int)
for g in mem.values(): N[g] += 1
print("denominators:", dict(N))

pres = defaultdict(lambda: defaultdict(set))
kept = 0
for line in open(W + "/panel.domtbl", errors="replace"):
    if line.startswith("#"): continue
    f = line.split()
    if len(f) < 23: continue
    fam, tlen = f[0], int(f[2])
    qname = f[3]
    ev_dom = float(f[12])
    hf, ht = int(f[15]), int(f[16])
    cov = (ht - hf + 1) / tlen if tlen else 0
    if ev_dom > 1e-10 or cov < 0.30: continue
    kept += 1
    genome = qname.split("|")[0]
    g = mem.get(genome)
    if g: pres[fam][g].add(genome)
print("domains passing E<=1e-10 and coverage>=0.30:", kept)

print()
print("%-16s %10s %12s %10s" % ("family", "candidate", "Akkermansia", "free"))
for fam in ("Sulfatase", "FGE-sulfatase", "Peptidase_M60", "Peptidase_M43"):
    row = []
    for g in ("candidate", "akkermansia", "free"):
        n = len(pres[fam][g])
        row.append("%d/%d (%.1f%%)" % (n, N[g], 100.0 * n / N[g]))
    print("%-16s %14s %18s %16s" % (fam, row[0], row[1], row[2]))

with open(W + "/mucin_pfam_prevalence.tsv", "w") as o:
    o.write("family\tgroup\tn_present\tn_total\tpct\n")
    for fam in ("Sulfatase", "FGE-sulfatase", "Peptidase_M60", "Peptidase_M43"):
        for g in ("candidate", "akkermansia", "free"):
            n = len(pres[fam][g])
            o.write("%s\t%s\t%d\t%d\t%.1f\n" % (fam, g, n, N[g], 100.0 * n / N[g]))
print()
print("WROTE", W + "/mucin_pfam_prevalence.tsv")
