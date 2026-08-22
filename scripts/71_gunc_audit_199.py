#!/usr/bin/env python3
# Chimerism screening results audited by arm
# Source: ch3-chitin-evolution/scripts/gunc_audit_199.py
# Output: results/gunc_199/gunc_audit_199.tsv
# GUNC chimerism audit across the 199 proteomes (105 candidate genus + 94 amphibian Akkermansia).
# Output: results/gunc_199/gunc_audit_199.tsv
#
# Two criteria, as in the Bacillota audit: GUNC's own pass threshold (CSS 0.45), and a strict
# call requiring CSS >= 0.85, reference representation >= 0.5, and max contamination at family
# level or above. Both are reported; the strict set is what would actually be excluded.

import csv, os, statistics, sys
from collections import Counter

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
GDIR = f"{BASE}/results/gunc_199"
CEN  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
CNT  = f"{BASE}/results/novel_akk_tree/cazy_per_genome_counts.tsv"
OUT  = f"{GDIR}/gunc_audit_199.tsv"

RRS_MIN  = 0.5
CSS_HARD = 0.85

sums = [f for f in os.listdir(GDIR) if f.startswith("GUNC.") and f.endswith(".maxCSS_level.tsv")]
if not sums:
    print("REFUSED: no GUNC.*.maxCSS_level.tsv in", GDIR)
    print("contents:", sorted(os.listdir(GDIR))[:10])
    sys.exit(1)
GUNC = f"{GDIR}/{sums[0]}"
print("reading", GUNC)

g = {}
with open(GUNC) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    I = {k: i for i, k in enumerate(h)}
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) <= max(I.values()):
            continue
        try:
            g[p[I["genome"]]] = dict(
                lvl=p[I["taxonomic_level"]],
                css=float(p[I["clade_separation_score"]]),
                cont=float(p[I["contamination_portion"]]),
                rrs=float(p[I["reference_representation_score"]]),
                ident=float(p[I["mean_hit_identity"]]),
                ncontig=int(p[I["n_contigs"]]),
                passed=p[I["pass.GUNC"]].strip().lower() == "true")
        except (ValueError, KeyError):
            continue
print("GUNC rows: %d" % len(g))

meta = {}
for r in csv.DictReader(open(CEN), delimiter="\t"):
    if r["family"] != "Akkermansiaceae":
        continue
    gen = (r["genus"] or "").strip()
    if (not gen or gen in ("unknown", "NO_GENUS")) and r["host_class"] == "amphibian":
        arm = "candidate"
    elif gen == "Akkermansia" and r["host_class"] == "amphibian":
        arm = "Akkermansia_amph"
    else:
        continue
    try:
        comp = float(r["completeness"])
    except ValueError:
        comp = float("nan")
    meta[r["accession"]] = dict(arm=arm, comp=comp, ds=r["from_dataset"])

nfam = {}
if os.path.exists(CNT):
    for r in csv.DictReader(open(CNT), delimiter="\t"):
        nfam[r["genome"]] = int(r["n_families"])

both = [k for k in g if k in meta]
print("joined to census: %d of %d GUNC rows" % (len(both), len(g)))
if not both:
    print("REFUSED: no join. GUNC genome IDs may differ from census accessions.")
    print("GUNC sample:", sorted(g)[:3])
    sys.exit(1)

def strict(k):
    return (g[k]["css"] >= CSS_HARD and g[k]["rrs"] >= RRS_MIN
            and g[k]["lvl"] in ("kingdom", "phylum", "class", "order", "family"))

print()
print("=== PASS RATE BY ARM (the comparison that matters) ===")
print("%-20s %-6s %-14s %-14s %-14s" % ("arm", "n", "pass.GUNC", "strict chimera", "RRS<0.5"))
for arm in ("candidate", "Akkermansia_amph"):
    sel = [k for k in both if meta[k]["arm"] == arm]
    if not sel:
        continue
    np_ = sum(1 for k in sel if g[k]["passed"])
    print("%-20s %-6d %-14s %-14s %-14s" % (
        arm, len(sel),
        "%d (%.1f%%)" % (np_, 100.0 * np_ / len(sel)),
        "%d (%.1f%%)" % (sum(1 for k in sel if strict(k)),
                         100.0 * sum(1 for k in sel if strict(k)) / len(sel)),
        "%d" % sum(1 for k in sel if g[k]["rrs"] < RRS_MIN)))
print()
print("If the two arms have similar pass rates, the retained free-living-like repertoire in the")
print("candidate genus is NOT explained by chimeric binning. That is the point of this run.")

print()
print("=== TAXONOMIC LEVEL OF MAX CSS ===")
for lvl, n in Counter(g[k]["lvl"] for k in both).most_common():
    fl = sum(1 for k in both if g[k]["lvl"] == lvl and not g[k]["passed"])
    print("  %-12s %4d | %3d fail" % (lvl, n, fl))

rrs = [g[k]["rrs"] for k in both]
print()
print("=== REFERENCE REPRESENTATION ===")
print("  RRS min %.3f  median %.3f  max %.3f" % (min(rrs), statistics.median(rrs), max(rrs)))
print("  RRS < %.1f: %d (%.1f%%). A low RRS means the lineage is undersampled in proGenomes,"
      % (RRS_MIN, sum(1 for x in rrs if x < RRS_MIN),
         100.0 * sum(1 for x in rrs if x < RRS_MIN) / len(rrs)))
print("  so a bad score there is weak evidence of real chimerism. Expect this to be high for a")
print("  genus with no named representative.")

fails = [k for k in both if not g[k]["passed"]]
if fails and nfam:
    fv = [nfam[k] for k in fails if k in nfam]
    pv = [nfam[k] for k in both if g[k]["passed"] and k in nfam]
    if fv and pv:
        print()
        print("=== DO GUNC FAILURES CARRY LARGER REPERTOIRES? ===")
        print("  median CAZy families: fail %.1f | pass %.1f" % (
            statistics.median(fv), statistics.median(pv)))
        print("  If failures carry MORE families, chimerism may be inflating the repertoire.")

print()
print("=== CheckM vs GUNC ===")
if fails:
    cf = [meta[k]["comp"] for k in fails if meta[k]["comp"] == meta[k]["comp"]]
    cp = [meta[k]["comp"] for k in both if g[k]["passed"] and meta[k]["comp"] == meta[k]["comp"]]
    if cf and cp:
        print("  median completeness: fail %.2f | pass %.2f" % (
            statistics.median(cf), statistics.median(cp)))
print("  GUNC detects chimerism CheckM misses, so low overlap is expected.")

with open(OUT, "w") as f:
    f.write("genome\tarm\tfrom_dataset\tcheckm_completeness\tn_cazy_families\t"
            "taxonomic_level\tclade_separation_score\tcontamination_portion\t"
            "reference_representation_score\tmean_hit_identity\tn_contigs\t"
            "pass_gunc\tstrict_chimera\n")
    for k in sorted(both):
        m, v = meta[k], g[k]
        f.write("%s\t%s\t%s\t%.2f\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%s\t%s\n"
                % (k, m["arm"], m["ds"], m["comp"], nfam.get(k, "NA"), v["lvl"],
                   v["css"], v["cont"], v["rrs"], v["ident"], v["ncontig"],
                   "pass" if v["passed"] else "FAIL", "yes" if strict(k) else "no"))
print()
print("wrote", OUT)
print("REPORT BOTH NUMBERS: the default pass.GUNC rate and the strict call with its criterion.")
