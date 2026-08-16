#!/usr/bin/env python3
# Family assignments joined onto the census
# Source: ch3-chitin-evolution/scripts/join_herptile_verru_family_v3.py
# Output: results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv
"""
Fill the missing family/genus calls for the 147 herptile Verru, then answer the question
the GH75 census could not:

  Is amphibian GH75 prevalence (62.1%) an AMPHIBIAN signal, or is it just that amphibians
  happen to carry GH75-bearing clades (mostly Akkermansiaceae, 94.1% GH75+ phylum-wide)?

NO GTDB-Tk RUN IS NEEDED. data/herptile_gtdbtk_compiled.tsv (May 27) already classifies
all 147. The unified catalog simply never joined family/genus back in, which is why the
census reported them as 'unclassified'. This script does that join.

THE TEST (indirect standardization)
  For each family f, take the GTDB_r226_rep GH75 prevalence p_f as the clade baseline
  (species-dereplicated, so genome-level is defensible there).
  Expected amphibian carriers = sum_f ( n_amphibian_f * p_f ).
  Compare to observed.
    expected ~= observed  -> amphibian "enrichment" is ENTIRELY clade composition. Dead.
    observed >> expected  -> amphibians carry GH75 above what their clades predict. Real.

WHAT THIS SCRIPT REFUSES TO DO
  1. Refuses to run a MAG-level significance test on the amphibian slice and call it
     significant. Effective replication is host animals (corrected unit, NOT sample_id
     which is the sequencing run for herptile_MAG; see host_animal()), not MAGs. The
     host-animal aggregation is printed instead, and the MAG-level numbers are labelled
     DESCRIPTIVE.
  2. Refuses to compare amphibian vs GTDB prevalence without holding completeness in view.
     Census showed prevalence climbs 20.2% -> 45.2% across completeness bins. If the two
     groups differ in completeness, a raw prevalence gap is partly an assembly artifact.
     A >=90%-complete-only rerun of the headline is printed.
  3. Refuses to report a family baseline from a handful of genomes. Families with fewer
     than MIN_BASELINE GTDB reps are pooled into an 'undercharacterized' bucket and the
     genomes they carry are reported separately, not silently given a fake baseline.

Run:
  sbatch /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/run_join_herptile_family.sh
"""
import csv, math
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_v3.tsv"
GTK  = f"{BASE}/data/herptile_gtdbtk_compiled.tsv"
OUT  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"

MIN_BASELINE = 10   # min GTDB reps in a family before its prevalence is used as a baseline

csv.field_size_limit(10_000_000)

def host_animal(r):
    # sample_id is the SEQUENCING RUN, not the animal, for herptile_MAG (UHM/STP): accession
    # UHM1228.41102_R.bin.103 has sample_id UHM1228.41102, but animal UHM1228 also appears as
    # sample_id UHM1228.23072. For EHI, sample_id already IS the animal. Corrected id matches
    # animal_metadata.tsv 97% (old sample_id unit: 0%). See CHATINDEX "THE UNIT OF REPLICATION".
    if r.get("from_dataset") == "herptile_MAG":
        return r["accession"].split(".")[0]
    return (r.get("sample_id") or "").strip()

def pct(n, d):
    return "  n/a" if d == 0 else f"{100.0*n/d:5.1f}"

def fisher_two_sided(a, b, c, d):
    n = a + b + c + d
    if n == 0 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
        return float("nan")
    def p_of(x):
        return (math.comb(a + b, x) * math.comb(c + d, a + c - x)) / math.comb(n, a + c)
    lo = max(0, a + c - (c + d))
    hi = min(a + b, a + c)
    p_obs = p_of(a)
    tot = 0.0
    for x in range(lo, hi + 1):
        p = p_of(x)
        if p <= p_obs * (1 + 1e-9):
            tot += p
    return min(1.0, tot)

# ----------------------------------------------------------------------------------
# 1. GTDB-Tk calls for herptile MAGs
# ----------------------------------------------------------------------------------
tk = {}
for r in csv.DictReader(open(GTK), delimiter="\t"):
    tk[r["user_genome"]] = (
        (r.get("family") or "").strip(),
        (r.get("genus") or "").strip(),
    )

def clean(rank):
    # 'f__Akkermansiaceae' -> 'Akkermansiaceae'; bare 'f__' -> '' (a real GTDB non-call)
    if not rank:
        return ""
    v = rank.split("__", 1)[-1] if "__" in rank else rank
    return v.strip()

# ----------------------------------------------------------------------------------
# 2. census + fill
# ----------------------------------------------------------------------------------
rows = list(csv.DictReader(open(CEN), delimiter="\t"))
filled = 0
nocall = 0
for r in rows:
    r["family_source"] = "catalog" if (r["family"] or "").strip() else ""
    if (r["family"] or "").strip():
        continue
    acc = r["accession"]
    if acc in tk:
        f, g = tk[acc]
        f, g = clean(f), clean(g)
        if f:
            r["family"] = f
            r["family_source"] = "gtdbtk_compiled"
            filled += 1
        else:
            r["family_source"] = "gtdbtk_no_family_call"
            nocall += 1
        if g and not (r["genus"] or "").strip():
            r["genus"] = g
    else:
        r["family_source"] = "unmatched"

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

ann = [r for r in rows if int(r["annotated"])]

print("=" * 84)
print("HERPTILE VERRU FAMILY JOIN -- from data/herptile_gtdbtk_compiled.tsv")
print("=" * 84)
print(f"family filled from GTDB-Tk:        {filled}")
print(f"GTDB-Tk returned NO family call:   {nocall}   (these stay unclassified, honestly)")
print(f"per-genome table with family ->    {OUT}")

herp = [r for r in ann if r["from_dataset"] == "herptile_MAG"]
print(f"\nherptile_MAG Verru family composition (n={len(herp)}):")
for k, n in Counter(r["family"] or "UNCLASSIFIED" for r in herp).most_common():
    p = sum(int(x["gh75_present"]) for x in herp if (x["family"] or "UNCLASSIFIED") == k)
    print(f"  {n:>4}  {k:<28} GH75+ {p:>3} = {pct(p,n)}%")

# ----------------------------------------------------------------------------------
# 3. baselines from GTDB reps
# ----------------------------------------------------------------------------------
gtdb = [r for r in ann if r["from_dataset"] == "GTDB_r226_rep"]
base = {}
for r in gtdb:
    f = r["family"] or "UNCLASSIFIED"
    base.setdefault(f, [0, 0])
    base[f][0] += 1
    base[f][1] += int(r["gh75_present"])

amph = [r for r in ann if r["host_class"] == "amphibian"]
obs = sum(int(r["gh75_present"]) for r in amph)

print("\n" + "=" * 84)
print("THE TEST: DOES CLADE COMPOSITION ALREADY EXPLAIN THE AMPHIBIAN NUMBER?")
print("=" * 84)
print("Baseline p_f = GH75 prevalence of family f among GTDB_r226_rep (species-dereplicated).")
print(f"Families with < {MIN_BASELINE} GTDB reps have NO usable baseline and are held out.\n")

print(f"{'family':<28}{'amph n':>8}{'amph+':>7}{'amph %':>8}{'GTDB n':>8}{'GTDB %':>8}{'exp+':>7}")
exp_total = 0.0
cover_n = 0
held = []
comp = Counter(r["family"] or "UNCLASSIFIED" for r in amph)
for f, n in comp.most_common():
    a_pos = sum(int(r["gh75_present"]) for r in amph if (r["family"] or "UNCLASSIFIED") == f)
    bn, bp = base.get(f, [0, 0])
    if bn < MIN_BASELINE:
        held.append((f, n, a_pos, bn))
        continue
    p_f = bp / bn
    e = n * p_f
    exp_total += e
    cover_n += n
    print(f"{f:<28}{n:>8}{a_pos:>7}{pct(a_pos,n):>8}{bn:>8}{pct(bp,bn):>8}{e:>7.1f}")

held_n = sum(h[1] for h in held)
held_pos = sum(h[2] for h in held)
print(f"\n{'HELD OUT (no usable baseline)':<28}{held_n:>8}{held_pos:>7}{pct(held_pos,held_n):>8}")
for f, n, a_pos, bn in held:
    print(f"    {f or 'UNCLASSIFIED':<26} amph n={n:<4} GH75+ {a_pos:<4} (only {bn} GTDB reps)")

cov_pos = sum(int(r["gh75_present"]) for r in amph
              if (r["family"] or "UNCLASSIFIED") in base
              and base[(r["family"] or "UNCLASSIFIED")][0] >= MIN_BASELINE)
print("\n" + "-" * 84)
print(f"  amphibian Verru with a usable clade baseline: {cover_n} of {len(amph)}")
print(f"  OBSERVED  GH75+ : {cov_pos:>5}  = {pct(cov_pos,cover_n)}%")
print(f"  EXPECTED  GH75+ : {exp_total:>5.1f}  = {pct(int(round(exp_total)),cover_n)}%   "
      f"(if amphibian clades behaved exactly like their GTDB relatives)")
if exp_total > 0:
    print(f"  observed / expected ratio: {cov_pos/exp_total:.2f}")
print()
print("  READ IT THIS WAY:")
print("    ratio ~= 1.0  -> the amphibian GH75 number is ENTIRELY clade composition.")
print("                     Amphibians are not enriched for GH75, they are enriched for")
print("                     Akkermansiaceae. The amphibian framing is dead.")
print("    ratio >> 1.0  -> amphibian genomes carry GH75 ABOVE their clade baseline.")
print("                     That is a real host-associated signal worth chasing.")
print("    ratio << 1.0  -> amphibian clades have LOST GH75 relative to relatives.")

# ----------------------------------------------------------------------------------
# 4. within-Akkermansiaceae: the single cleanest contrast
# ----------------------------------------------------------------------------------
print("\n" + "=" * 84)
print("WITHIN-AKKERMANSIACEAE: amphibian vs GTDB reference. Clade held fixed.")
print("=" * 84)
aa = [r for r in amph if r["family"] == "Akkermansiaceae"]
ga = [r for r in gtdb if r["family"] == "Akkermansiaceae"]
if aa and ga:
    pa = sum(int(r["gh75_present"]) for r in aa)
    pg = sum(int(r["gh75_present"]) for r in ga)
    p = fisher_two_sided(pa, len(aa) - pa, pg, len(ga) - pg)
    print(f"  amphibian Akkermansiaceae : {pa:>4}/{len(aa):<4} = {pct(pa,len(aa))}%")
    print(f"  GTDB rep Akkermansiaceae  : {pg:>4}/{len(ga):<4} = {pct(pg,len(ga))}%")
    print(f"  Fisher p = {p:.3g}   (DESCRIPTIVE: amphibian MAGs are not independent)")
    print("  If these are the same, amphibian Akkermansia are ordinary Akkermansia.")

# ----------------------------------------------------------------------------------
# 5. completeness control
# ----------------------------------------------------------------------------------
print("\n" + "=" * 84)
print("COMPLETENESS CONTROL -- is any gap just assembly quality?")
print("=" * 84)
def mean_comp(v):
    c = [float(r["completeness"]) for r in v if r["completeness"] not in ("", "NA")]
    return sum(c) / len(c) if c else float("nan")
print(f"  mean completeness, amphibian Verru : {mean_comp(amph):.1f}%")
print(f"  mean completeness, GTDB rep Verru  : {mean_comp(gtdb):.1f}%")
hi_a = [r for r in amph if r["completeness"] not in ("", "NA") and float(r["completeness"]) >= 90]
hi_g = [r for r in gtdb if r["completeness"] not in ("", "NA") and float(r["completeness"]) >= 90]
if hi_a and hi_g:
    pa = sum(int(r["gh75_present"]) for r in hi_a)
    pg = sum(int(r["gh75_present"]) for r in hi_g)
    print(f"\n  RESTRICTED TO >=90% COMPLETE:")
    print(f"    amphibian : {pa:>4}/{len(hi_a):<4} = {pct(pa,len(hi_a))}%")
    print(f"    GTDB rep  : {pg:>4}/{len(hi_g):<4} = {pct(pg,len(hi_g))}%")

# ----------------------------------------------------------------------------------
# 6. host-animal aggregation (the ONLY defensible unit for the amphibian slice)
# ----------------------------------------------------------------------------------
print("\n" + "=" * 84)
print("HOST-ANIMAL AGGREGATION -- effective replication, not MAG count")
print("=" * 84)
by_animal = defaultdict(list)
for r in amph:
    a = host_animal(r)
    if a:
        by_animal[a].append(r)
print(f"  amphibian Verru MAGs: {len(amph)}  behind  {len(by_animal)} distinct host animals")
if by_animal:
    carr = sum(1 for v in by_animal.values() if any(int(r["gh75_present"]) for r in v))
    print(f"  animals with >=1 GH75+ Verru MAG: {carr}/{len(by_animal)} = {pct(carr,len(by_animal))}%")
    akk_only = sum(1 for v in by_animal.values()
                   if any(int(r["gh75_present"]) and r["family"] == "Akkermansiaceae" for r in v))
    print(f"  animals whose GH75 comes from Akkermansiaceae: {akk_only}/{len(by_animal)} = "
          f"{pct(akk_only,len(by_animal))}%")
    print("  This is the number to quote for amphibians. Not the MAG-level percentage.")

wild = [r for r in amph if r["captivity_status"] == "wild"]
if wild:
    pw = sum(int(r["gh75_present"]) for r in wild)
    print(f"\n  WILD-ONLY amphibian Verru: {pw}/{len(wild)} = {pct(pw,len(wild))}%")
    wa = defaultdict(list)
    for r in wild:
        a = host_animal(r)
        if a:
            wa[a].append(r)
    if wa:
        c = sum(1 for v in wa.values() if any(int(r["gh75_present"]) for r in v))
        print(f"  wild animals with >=1 GH75+ Verru: {c}/{len(wa)} = {pct(c,len(wa))}%")

print("\n" + "=" * 84)
