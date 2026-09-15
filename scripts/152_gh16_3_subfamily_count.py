#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# GH16_3 subfamily prevalence, groups and dbCAN inputs mirrored from repo script 42.
import csv, os, sys
BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM = "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
CEN = BASE + "/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
DIRS = [BASE + "/results/dbcan_verru", BASE + "/results/dbcan_bact_refs",
        BASE + "/results/dbcan_ehi_amphibian", BASE + "/results/dbcan_ehi_nonamph",
        BASE + "/results/dbcan_endo_allphyla", BASE + "/results/dbcan_flavo_refs",
        BASE + "/results/dbcan_scaffold", RUM]
OUT = BASE + "/results/novel_akk_tree/gh16_3_prevalence.tsv"

rows = [r for r in csv.DictReader(open(CEN), delimiter="\t")
        if r["family"] == "Akkermansiaceae" and r["annotated"] == "1"]
def gen(r): return (r["genus"] or "").strip()
def isnov(r): return gen(r) in ("", "unknown", "NO_GENUS")
G = {}
G["NOVEL"]      = [r for r in rows if r["host_class"] == "amphibian" and isnov(r)]
G["AKK_AMPH"]   = [r for r in rows if r["host_class"] == "amphibian" and gen(r) == "Akkermansia"]
G["AKK_MAMMAL"] = [r for r in rows if r["host_class"] == "mammal"    and gen(r) == "Akkermansia"]
G["AKK_GTDB"]   = [r for r in rows if r["host_class"] == "unknown"   and gen(r) == "Akkermansia"]
G["AKK_REPT"]   = [r for r in rows if r["host_class"] == "reptile"   and gen(r) == "Akkermansia"]
sizes = {k: len(v) for k, v in G.items()}
print("group sizes:", sizes)
exp = {"NOVEL": 105, "AKK_AMPH": 94, "AKK_MAMMAL": 71, "AKK_GTDB": 42, "AKK_REPT": 137}
if sizes != exp:
    sys.exit("REFUSED: group sizes do not match script 42 (expected %s). isnov() may differ; paste lines 58 to 68 of 42_novel_genus_function.py." % exp)

path = {}
for d in DIRS:
    if not os.path.isdir(d): print("  MISSING DIR", d); continue
    for fn in os.listdir(d):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            path.setdefault(fn[:-4], os.path.join(d, fn))
print("dbCAN tsv files indexed:", len(path))

def fams(acc):
    p = path.get(acc)
    if p is None and acc[:3] in ("GB_", "RS_"): p = path.get(acc[3:])
    if p is None: p = path.get("GB_" + acc) or path.get("RS_" + acc)
    if p is None and "." in acc: p = path.get(acc.rsplit(".", 1)[0])
    if p is None: return None
    out = set()
    for line in open(p):
        h = line.split("\t", 1)[0]
        if h.endswith(".hmm"): out.add(h[:-4])
    return out

res, missing = {}, []
for k, members in G.items():
    n16_3 = n16any = n = 0
    for r in members:
        acc = r.get("accession") or r.get("genome") or r.get("assembly")
        f = fams(acc)
        if f is None: missing.append((k, acc)); continue
        n += 1
        if "GH16_3" in f: n16_3 += 1
        if any(x == "GH16" or x.startswith("GH16_") for x in f): n16any += 1
    res[k] = (n, n16_3, n16any)

print("\n%-11s %5s %8s %7s %9s %7s" % ("group", "n", "GH16_3", "pct", "GH16_any", "pct"))
for k in ("NOVEL", "AKK_AMPH", "AKK_REPT", "AKK_MAMMAL", "AKK_GTDB"):
    n, a, b = res[k]
    print("%-11s %5d %8d %6.1f%% %9d %6.1f%%" % (k, n, a, 100.0*a/n if n else 0, b, 100.0*b/n if n else 0))
if missing:
    print("\ngenomes with no dbCAN tsv: %d (first 5: %s)" % (len(missing), missing[:5]))
with open(OUT, "w") as o:
    o.write("group\tn_with_dbcan\tn_GH16_3\tpct_GH16_3\tn_GH16_any\tpct_GH16_any\n")
    for k in ("NOVEL", "AKK_AMPH", "AKK_REPT", "AKK_MAMMAL", "AKK_GTDB"):
        n, a, b = res[k]
        o.write("%s\t%d\t%d\t%.1f\t%d\t%.1f\n" % (k, n, a, 100.0*a/n if n else 0, b, 100.0*b/n if n else 0))
print("\nwrote", OUT)
