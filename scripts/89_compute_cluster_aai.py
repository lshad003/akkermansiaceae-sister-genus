#!/usr/bin/env python3
# Cluster-to-cluster identity matrix computed from reciprocal best hits
# Source: ch3-chitin-evolution/scripts/compute_cluster_aai.py
# Output: results/cluster_aai/cluster_aai_report.txt
import os
import itertools
import collections

OUT = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/cluster_aai"
H = OUT + "/hits"
RES = OUT + "/cluster_aai_matrix.tsv"

CL = ["C286","C287","C288","C289","C290","C291","C292","C293","C294",
      "C295","C296","C297","C298","C299","C300","C301","C302"]

LOW = {"C286","C294","C296","C298","C300"}

files = os.listdir(H)
print("hit files:", len(files), " expected 272")
missing = []
for a, b in itertools.permutations(CL, 2):
    if not os.path.exists("%s/%s_vs_%s.tsv" % (H, a, b)):
        missing.append("%s_vs_%s" % (a, b))
if missing:
    print("MISSING:", len(missing), missing[:10])
    print("STOP")
    raise SystemExit(1)
print("all present")
print("")

def besthits(a, b):
    d = {}
    for line in open("%s/%s_vs_%s.tsv" % (H, a, b)):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        q, s = f[0], f[1]
        try:
            p = float(f[2])
        except ValueError:
            continue
        if q not in d:
            d[q] = (s, p)
    return d

aai = {}
rbh = {}
for a, b in itertools.combinations(CL, 2):
    fwd = besthits(a, b)
    rev = besthits(b, a)
    ps = []
    for q, (s, p) in fwd.items():
        r = rev.get(s)
        if r and r[0] == q:
            ps.append((p + r[1]) / 2.0)
    if ps:
        aai[(a, b)] = sum(ps) / len(ps)
        rbh[(a, b)] = len(ps)
    else:
        aai[(a, b)] = float("nan")
        rbh[(a, b)] = 0

def get(a, b):
    return aai.get((a, b), aai.get((b, a), float("nan")))

print("=== AAI MATRIX (reciprocal best hits, mean identity) ===")
print("        " + "".join("%7s" % c for c in CL))
for a in CL:
    row = "%-8s" % a
    for b in CL:
        row += "   ----" if a == b else "%7.2f" % get(a, b)
    print(row)
print("")

with open(RES, "w") as fh:
    fh.write("clusterA\tclusterB\tAAI_pct\tn_RBH\tgroupA\tgroupB\tcomparison\n")
    for a, b in itertools.combinations(CL, 2):
        ga = "LOW" if a in LOW else "HIGH"
        gb = "LOW" if b in LOW else "HIGH"
        comp = "within_LOW" if ga == gb == "LOW" else ("within_HIGH" if ga == gb == "HIGH" else "between")
        fh.write("%s\t%s\t%.2f\t%d\t%s\t%s\t%s\n" % (a, b, get(a, b), rbh[(a, b)], ga, gb, comp))
print("WROTE:", RES)
print("")

groups = collections.defaultdict(list)
for a, b in itertools.combinations(CL, 2):
    ga = "LOW" if a in LOW else "HIGH"
    gb = "LOW" if b in LOW else "HIGH"
    comp = "within_LOW" if ga == gb == "LOW" else ("within_HIGH" if ga == gb == "HIGH" else "between")
    groups[comp].append(get(a, b))

print("=== DISTRIBUTIONS ===")
print("LOW group = C286 C294 C296 C298 C300 (the five below 65 to the type genome)")
for k in ("within_LOW", "within_HIGH", "between"):
    v = sorted(x for x in groups[k] if x == x)
    if not v:
        continue
    print("%-12s n=%3d  min %.2f  median %.2f  max %.2f" % (k, len(v), v[0], v[len(v)//2], v[-1]))
print("")

print("=== IS THERE A GAP? all 136 pairwise values, sorted ===")
allv = sorted(x for x in (get(a, b) for a, b in itertools.combinations(CL, 2)) if x == x)
prev = None
biggest = (0.0, None)
for v in allv:
    if prev is not None and v - prev > biggest[0]:
        biggest = (v - prev, (prev, v))
    prev = v
print("range %.2f to %.2f" % (allv[0], allv[-1]))
print("largest gap: %.2f, between %.2f and %.2f" % (biggest[0], biggest[1][0], biggest[1][1]))
print("")
print("values below 65:", sum(1 for v in allv if v < 65), "of", len(allv))
print("all values:")
print("  " + "  ".join("%.1f" % v for v in allv))
