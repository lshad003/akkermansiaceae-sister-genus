#!/usr/bin/env python3
# Genus-level novelty confirmed from classifier RED values
# Source: ch3-chitin-evolution/scripts/step3_red_values.py
# Output: results/novel_akk_tree/red_values_54.txt
import csv
import collections

GT = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/herptile_gtdbtk_compiled.tsv"
CEN = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"

rows = list(csv.DictReader(open(CEN), delimiter="\t"))
nov = [r for r in rows if r["family"] == "Akkermansiaceae"
       and (not r["genus"].strip() or r["genus"] == "unknown")
       and r["host_class"].lower() in ("amphibia", "amphibian")]
acc = set(r["accession"] for r in nov)
print("candidate set from census:", len(acc))
print("  by dataset:", dict(collections.Counter(r["from_dataset"] for r in nov)))
print("")

gt = {}
for r in csv.DictReader(open(GT), delimiter="\t"):
    gt[r["user_genome"]] = r

hit = {a: gt[a] for a in acc if a in gt}
print("candidates found in herptile_gtdbtk_compiled.tsv:", len(hit), "of", len(acc))
print("")

if not hit:
    print("STOP: no overlap. Check whether user_genome uses a different accession form.")
    print("first 3 user_genome values:", list(gt)[:3])
    print("first 3 candidate accessions:", list(acc)[:3])
    raise SystemExit(0)

red = []
for a, r in hit.items():
    try:
        red.append(float(r["red_value"]))
    except (ValueError, TypeError):
        pass

red.sort()
print("RED VALUES, n =", len(red))
if red:
    print("  min %.4f  median %.4f  max %.4f" % (red[0], red[len(red)//2], red[-1]))
print("")

print("classification_method:", dict(collections.Counter(r["classification_method"] for r in hit.values())))
print("")
cgr = collections.Counter((r["closest_genome_reference"] or "<blank>") for r in hit.values())
print("closest_genome_reference, top 5:", dict(cgr.most_common(5)))
print("")

def gslot(s):
    for t in s.split(";"):
        t = t.strip()
        if t.startswith("g__"):
            return t[3:] or "<EMPTY g__>"
    return "<NO g__ FIELD>"

print("g__ slot:", dict(collections.Counter(gslot(r["classification"]) for r in hit.values())))
print("")
print("family slot:", dict(collections.Counter(r["family"] for r in hit.values())))
print("")
print("EXAMPLE ROW")
k = sorted(hit)[0]
for f in ("user_genome", "classification", "classification_method",
          "closest_genome_reference", "closest_genome_ani", "red_value", "warnings"):
    print("  %-26s %s" % (f, hit[k].get(f)))
