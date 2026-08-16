#!/usr/bin/env python3
# Host animals recounted under the corrected animal unit
# Source: ch3-chitin-evolution/scripts/recount_host_animals_corrected_unit.py
# Output: results/akk_composition/host_animal_recount_corrected_unit.tsv
"""
Recompute EVERY host-animal number in CHATINDEX under BOTH units, side by side.

WHY (found 2026-07-14 by guard 3 of rarefy_gh75_amph_reptile.py):
sample_id is NOT the host animal for herptile_MAG. Accession UHM1228.41102_R.bin.103
has sample_id UHM1228.41102, and the SAME animal UHM1228 also appears as sample_id
UHM1228.23072. One animal, two sequencing runs, counted twice.
For EHI_2025 the semantics are inverted: sample_id (AAP82) IS the animal, and the
accession prefix (EHM034581) is a per-genome id.

CORRECTED UNIT: herptile_MAG -> token before the first '.' in accession.
                everything else -> catalog sample_id.

Every number printed as "sample_id" is what CHATINDEX currently says.
Every number printed as "corrected" is what it should say.
"""
import csv, collections, sys

csv.field_size_limit(10_000_000)
BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
V2   = f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
CAT  = f"{BASE}/data/unified_annotation/unified_genome_annotation.tsv"
OUT  = f"{BASE}/results/akk_composition/host_animal_recount_corrected_unit.tsv"

try:
    from scipy.stats import fisher_exact
except ImportError:
    sys.exit("FATAL: scipy unavailable.")

def hr(t=""):
    print("=" * 94)
    if t:
        print(t)
        print("=" * 94)

cat  = {r["accession"]: r for r in csv.DictReader(open(CAT), delimiter="\t")}
rows = list(csv.DictReader(open(V2), delimiter="\t"))
ann  = [r for r in rows if r["annotated"] == "1"]

sid = lambda r: (cat.get(r["accession"], {}).get("sample_id") or "").strip()
def animal(r):
    if r["from_dataset"] == "herptile_MAG":
        return r["accession"].split(".")[0]
    return sid(r)

nu = lambda rs, key: len(set(k for k in (key(r) for r in rs) if k))
def anyrule(rs, key):
    """host animal is GH75+ if ANY of its MAGs carries GH75."""
    d = collections.defaultdict(int)
    for r in rs:
        k = key(r)
        if k:
            d[k] |= (r["gh75_present"] == "1")
    return sum(d.values()), len(d)

out = []
def rec(section, cell, metric, sample_id_val, corrected_val):
    out.append((section, cell, metric, sample_id_val, corrected_val))

# ---------------------------------------------------------------- sanity
hr("SANITY: the genome-level census is UNAFFECTED by the unit bug")
akk_all = [r for r in ann if r["family"] == "Akkermansiaceae"]
print(f"  annotated Verru {len(ann)} (expect 3864), GH75+ "
      f"{sum(1 for r in ann if r['gh75_present']=='1')} (expect 1319)")
print(f"  Akkermansiaceae annotated {len(akk_all)} (expect 742), GH75+ "
      f"{sum(1 for r in akk_all if r['gh75_present']=='1')} (expect 646)")
print("  Result 4 (census), Result 3 (Q3) and the derep species counts are GENOME/SPECIES level.")
print("  They do not use host animals and are NOT touched by this correction.")

# ---------------------------------------------------------------- Result 7, Akkermansiaceae
hr("RESULT 7: AKKERMANSIACEAE COMPOSITION, animals per host_class")
print(f"  {'host_class':<14}{'genomes':>8}{'annot':>7}{'sample_id':>11}{'corrected':>11}"
      f"{'lost':>6}{'gen/animal':>12}")
print("  " + "-" * 72)
akk_raw = [r for r in rows if r["family"] == "Akkermansiaceae"]
by = collections.defaultdict(list)
for r in akk_raw:
    by[r["host_class"]].append(r)
for k in sorted(by, key=lambda k: -len(by[k])):
    rs = by[k]
    a  = [r for r in rs if r["annotated"] == "1"]
    ns, nc = nu(rs, sid), nu(rs, animal)
    gpa = f"{len(rs)/nc:.2f}" if nc else "-"
    print(f"  {k:<14}{len(rs):>8}{len(a):>7}{ns:>11}{nc:>11}{ns-nc:>6}{gpa:>12}")
    if ns or nc:
        rec("Result7_Akkermansiaceae", k, "host_animals", ns, nc)
ts, tc = nu(akk_raw, sid), nu(akk_raw, animal)
print(f"  {'TOTAL':<14}{len(akk_raw):>8}{len(akk_all):>7}{ts:>11}{tc:>11}{ts-tc:>6}")
rec("Result7_Akkermansiaceae", "TOTAL", "host_animals", ts, tc)

# ---------------------------------------------------------------- Result 7, all Verru
hr("RESULT 7: ALL VERRUCOMICROBIOTA, animals per host_class  (annot column is now v2)")
print(f"  {'host_class':<14}{'genomes':>8}{'annot':>7}{'sample_id':>11}{'corrected':>11}{'lost':>6}")
print("  " + "-" * 60)
byv = collections.defaultdict(list)
for r in rows:
    byv[r["host_class"]].append(r)
for k in sorted(byv, key=lambda k: -len(byv[k])):
    rs = byv[k]
    a  = [r for r in rs if r["annotated"] == "1"]
    ns, nc = nu(rs, sid), nu(rs, animal)
    print(f"  {k:<14}{len(rs):>8}{len(a):>7}{ns:>11}{nc:>11}{ns-nc:>6}")
    if ns or nc:
        rec("Result7_allVerru", k, "host_animals", ns, nc)
ts, tc = nu(rows, sid), nu(rows, animal)
print(f"  {'TOTAL':<14}{len(rows):>8}{len(ann):>7}{ts:>11}{tc:>11}{ts-tc:>6}")
rec("Result7_allVerru", "TOTAL", "host_animals", ts, tc)
print("\n  NOTE: CHATINDEX Result 7 still quotes the v1 annotated counts for ALL VERRU")
print("  (3,632 annotated; reptile 167/27). The annot column above is v2 and supersedes them.")

# ---------------------------------------------------------------- Result 5(A)
hr("RESULT 5(A): amphibian Verru GH75 at host-animal level (any-member rule)")
amv  = [r for r in ann if r["host_class"] == "amphibian"]
amvw = [r for r in amv if r["captivity_status"] == "wild"]
for tag, rs in (("all animals", amv), ("wild-only", amvw)):
    p1, n1 = anyrule(rs, sid)
    p2, n2 = anyrule(rs, animal)
    print(f"  {tag:<12} sample_id {p1:>3}/{n1:<3} = {100*p1/n1:>5.1f}%"
          f"      corrected {p2:>3}/{n2:<3} = {100*p2/n2:>5.1f}%")
    rec("Result5A_amph_Verru", tag, "gh75_pos_animals", f"{p1}/{n1}", f"{p2}/{n2}")
    rec("Result5A_amph_Verru", tag, "pct", f"{100*p1/n1:.1f}", f"{100*p2/n2:.1f}")
print("\n  CHATINDEX currently says 95/113 = 84.1% and wild-only 94/111 = 84.7%.")

akkp = [r for r in amv if r["gh75_present"] == "1" and r["family"] == "Akkermansiaceae"]
allp = [r for r in amv if r["gh75_present"] == "1"]
print(f"\n  'EVERY amphibian GH75 is Akkermansiaceae' check: {len(akkp)}/{len(allp)} GH75+ amphibian")
print(f"  Verru MAGs are Akkermansiaceae. Genome-level, unaffected by the unit. (CHATINDEX: 159/159)")
pa, na = anyrule([r for r in amv if r["family"] == "Akkermansiaceae"], animal)
print(f"  animals whose GH75 comes from Akkermansiaceae, corrected unit: {pa}/{na}")

# ---------------------------------------------------------------- Result 8
hr("RESULT 8: the depth-biased any-member test, both units (documents the bias, does not fix it)")
akk_ann = [r for r in ann if r["family"] == "Akkermansiaceae"]
amph = [r for r in akk_ann if r["host_class"] == "amphibian"]
rept = [r for r in akk_ann if r["host_class"] == "reptile"]
for label, key in (("sample_id", sid), ("corrected", animal)):
    ap, an = anyrule(amph, key)
    rp, rn = anyrule(rept, key)
    _, p = fisher_exact([[ap, an - ap], [rp, rn - rp]])
    print(f"  {label:<10} amphibian {ap:>3}/{an:<3} = {100*ap/an:>5.1f}%   "
          f"reptile {rp:>3}/{rn:<3} = {100*rp/rn:>5.1f}%   Fisher p = {p:.4g}")
    rec("Result8_anymember", "amphibian", "gh75_pos_animals", f"{ap}/{an}", "")
    rec("Result8_anymember", "reptile", "gh75_pos_animals", f"{rp}/{rn}", "")
    rec("Result8_anymember", "fisher_p", label, f"{p:.4g}", "")
print("\n  BOTH of these are STILL depth-biased and NEITHER is quotable. The any-member rule gives")
print("  amphibians ~2 draws per animal and reptiles ~1 regardless of which unit is used.")
print("  The valid test is the rarefaction: rarefaction_amph_vs_reptile_gh75.tsv. It is NULL.")

# ---------------------------------------------------------------- MAGs per animal
hr("MAGs PER ANIMAL under the corrected unit (the depth asymmetry that motivated rarefaction)")
for k, rs in (("amphibian", amph), ("reptile", rept)):
    per = collections.Counter(animal(r) for r in rs if animal(r))
    dist = collections.Counter(per.values())
    print(f"  {k:<10} {len(rs)} MAGs / {len(per)} animals = {len(rs)/len(per):.2f} MAGs per animal")
    print(f"  {'':<10} " + "  ".join(f"{n} MAG(s): {c}" for n, c in sorted(dist.items())))
    rec("depth", k, "mags_per_animal_corrected", f"{len(rs)/len(per):.2f}", "")

with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["section", "cell", "metric", "sample_id_unit", "corrected_unit"])
    for r in out:
        w.writerow(r)
print(f"\n  -> {OUT}")
