#!/usr/bin/env python3
# Candidate set counted before and after the non-amphibian exclusion
# Source: ch3-chitin-evolution/scripts/check_novel_107_vs_105.py
# Output: results/novel_akk_tree/novel_107_vs_105.txt
import csv

P = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
rows = list(csv.DictReader(open(P), delimiter="\t"))

nov = [r for r in rows if r["family"] == "Akkermansiaceae"
       and (not r["genus"].strip() or r["genus"] == "unknown")]

print("novel-candidate rows in census:", len(nov))
print()
print("%-34s %-14s %-16s %-12s %s" % ("accession", "from_dataset", "host_class", "genome_type", "host_animal_type"))
for r in sorted(nov, key=lambda x: x["host_class"]):
    if r["host_class"].lower() not in ("amphibia", "amphibian"):
        print("%-34s %-14s %-16s %-12s %s" % (
            r["accession"], r["from_dataset"], r["host_class"] or "(blank)",
            r["genome_type"], r["host_animal_type"] or "(blank)"))
print()
import collections
print("host_class counts:", dict(collections.Counter(r["host_class"] for r in nov)))
