#!/usr/bin/env python3
# Host range tested against the full census and the independent EHI collection
# Source: ch3-chitin-evolution/scripts/host_range_test.py
# Output: results/host_range/host_range_test.txt
import csv
import collections
import os

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
CEN = CH3 + "results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
EHI = CH3 + "data/ehi_2025_annotated/ehi_mags_annotated_v2.tsv"
OUT = CH3 + "results/host_range/host_range_test.txt"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

L = []
def say(s):
    print(s)
    L.append(s)

say("HOST RANGE TEST for the candidate genus")
say("")

rows = list(csv.DictReader(open(CEN), delimiter="\t"))
say("=== 1. FULL CENSUS: %s ===" % CEN)
say("rows: %d" % len(rows))
akk = [r for r in rows if r["family"] == "Akkermansiaceae"]
say("Akkermansiaceae genomes: %d" % len(akk))
say("")

def blank(v):
    return (v or "").strip() in ("", "unknown", "NO_GENUS")

cand = [r for r in akk if blank(r["genus"])]
say("unassigned-genus Akkermansiaceae: %d" % len(cand))
c = collections.Counter(r["host_class"] for r in cand)
for k, v in c.most_common():
    say("   %-12s %d" % (k, v))
say("")
say("non-amphibian candidates, listed:")
for r in cand:
    if r["host_class"].lower() not in ("amphibia", "amphibian"):
        say("   %-34s %-10s %-14s %s" % (r["accession"], r["host_class"],
            r.get("from_dataset", ""), r.get("host_animal_type", "")))
say("")

say("Akkermansia itself, by host class, same census:")
named = [r for r in akk if (r["genus"] or "").strip() == "Akkermansia"]
for k, v in collections.Counter(r["host_class"] for r in named).most_common():
    say("   %-12s %d" % (k, v))
say("   total %d" % len(named))
say("")

say("candidates by source dataset:")
for k, v in collections.Counter(r.get("from_dataset", "?") for r in cand).most_common():
    say("   %-16s %d" % (k, v))
say("")

say("=== 2. EHI INDEPENDENT COLLECTION: %s ===" % EHI)
e = list(csv.DictReader(open(EHI), delimiter="\t"))
say("EHI MAGs: %d" % len(e))
for k, v in collections.Counter(r["host_class"] for r in e).most_common():
    say("   %-12s %d" % (k, v))
say("")
ea = [r for r in e if "Akkermansiaceae" in (r["family"] or "")]
say("EHI Akkermansiaceae: %d" % len(ea))
ec = [r for r in ea if blank(r["genus"])]
say("EHI unassigned-genus Akkermansiaceae: %d" % len(ec))
say("   non-amphibian among them: %d" % sum(1 for r in ec
    if r["host_class"].lower() not in ("amphibia", "amphibian")))
say("")
say("EHI candidate host species:")
for k, v in collections.Counter(r.get("host_animal_type", "?") for r in ec).most_common():
    say("   %-28s %d" % (k, v))
say("")
say("EHI Akkermansia by host class, showing the family is sampled broadly:")
en = [r for r in ea if (r["genus"] or "").strip() == "Akkermansia"]
for k, v in collections.Counter(r["host_class"] for r in en).most_common():
    say("   %-12s %d" % (k, v))

open(OUT, "w").write("\n".join(L) + "\n")
print("")
print("WROTE:", OUT)
