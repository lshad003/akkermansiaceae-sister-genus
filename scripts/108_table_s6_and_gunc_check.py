#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# 1. Does removing the GUNC failures change the candidate CAZyme median or IQR?
# 2. Build TableS6 from cazy_by_genus.tsv.
import csv, os, shutil, statistics as st
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
def q(v, p):
    v = sorted(v); k = (len(v) - 1) * p; f = int(k)
    return v[f] if f + 1 >= len(v) else v[f] + (k - f) * (v[f + 1] - v[f])

cnt = {r["genome"]: float(r["n_families"]) for r in
       csv.DictReader(open(B + "/results/novel_akk_tree/cazy_per_genome_counts.tsv"), delimiter="\t")
       if r["group"] == "novel"}
fail = set()
for r in csv.DictReader(open(B + "/results/gunc_199/gunc_audit_199.tsv"), delimiter="\t"):
    if r["arm"] == "candidate" and r["pass_gunc"].lower() not in ("true", "yes", "1", "pass"):
        fail.add(r["genome"])
print("candidate genomes: %d   GUNC failures: %d  %s" % (len(cnt), len(fail), sorted(fail)))
allv = list(cnt.values())
keep = [v for g, v in cnt.items() if g not in fail]
for lab, v in (("all 105", allv), ("GUNC pass only", keep)):
    print("  %-16s n=%3d  median %.1f  IQR %.1f to %.1f" % (lab, len(v), st.median(v), q(v, .25), q(v, .75)))

src = B + "/results/novel_akk_tree/cazy_by_genus.tsv"
dst = B + "/results/tables/TableS6_cazy_by_genus.tsv"
os.makedirs(os.path.dirname(dst), exist_ok=True)
shutil.copy(src, dst)
n = sum(1 for _ in open(dst)) - 1
print("\nwrote %s  (%d families)" % (dst, n))
