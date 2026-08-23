#!/usr/bin/env python3
# Functional categories of the ancestral losses summarized
# Source: ch3-chitin-evolution/scripts/summarize_polarity_B.py
# Output: results/pangenome/polarity_B_summary.tsv
import csv, collections, os, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
EGG  = f"{BASE}/results/pangenome/eggnog_polarity"
OUT  = f"{BASE}/results/pangenome/polarity_B_summary.tsv"

COG = {
 "J":"translation, ribosome","A":"RNA processing","K":"transcription","L":"replication and repair",
 "B":"chromatin","D":"cell cycle and division","Y":"nuclear structure","V":"defence",
 "T":"signal transduction","M":"cell wall, membrane, envelope","N":"cell motility",
 "Z":"cytoskeleton","W":"extracellular structures","U":"trafficking and secretion",
 "O":"protein turnover and chaperones","X":"mobilome, prophages, transposons",
 "C":"energy production","G":"carbohydrate transport and metabolism",
 "E":"amino acid transport and metabolism","F":"nucleotide transport and metabolism",
 "H":"coenzyme transport and metabolism","I":"lipid transport and metabolism",
 "P":"inorganic ion transport","Q":"secondary metabolites","R":"general function only",
 "S":"function unknown"}

def load(cat):
    p = f"{EGG}/{cat}.emapper.annotations"
    if not os.path.exists(p):
        return None
    rows = []
    hdr = None
    for line in open(p):
        if line.startswith("#query"):
            hdr = line.lstrip("#").rstrip("\n").split("\t")
            continue
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if hdr and len(f) >= len(hdr):
            rows.append(dict(zip(hdr, f)))
    return rows

TOT = {"B_lost_at_gut_ancestor": 1047, "A_lost_in_Akkermansia": 57,
       "C_Akkermansia_enriched": 34, "D_novel_specific": 55}

print("%-26s %8s %10s %8s" % ("category", "total", "annotated", "rate"))
data = {}
for cat in ("B_lost_at_gut_ancestor", "A_lost_in_Akkermansia",
            "C_Akkermansia_enriched", "D_novel_specific"):
    r = load(cat)
    if r is None:
        print("%-26s %8s %10s" % (cat, TOT.get(cat, "?"), "NOT RUN")); continue
    data[cat] = r
    n = len(set(x["query"] for x in r))
    print("%-26s %8d %10d %7.1f%%" % (cat, TOT[cat], n, 100.0*n/TOT[cat]))

B = data.get("B_lost_at_gut_ancestor", [])
if not B:
    print("REFUSED: no category B rows"); sys.exit(1)

cats = collections.Counter()
for r in B:
    c = (r.get("COG_category") or "-").strip()
    if c in ("", "-"):
        cats["(none)"] += 1; continue
    for ch in c:                      # multi-letter assignments count once each
        cats[ch] += 1

print()
print("=== COG categories among the %d annotated category-B orthogroups ===" % len(B))
print("(orthogroups with several letters are counted under each)")
print("%-5s %-38s %6s %7s" % ("code", "category", "n", "pct"))
for k, v in cats.most_common(24):
    print("%-5s %-38s %6d %6.1f%%" % (k, COG.get(k, "unclassified" if k != "(none)" else "no COG"),
                                      v, 100.0*v/len(B)))

known = sum(v for k, v in cats.items() if k not in ("S", "R", "(none)"))
print()
print("assigned to a specific function (excluding S, R and none): %d of %d (%.1f%%)"
      % (known, len(B), 100.0*known/len(B)))

named = [r for r in B if (r.get("Preferred_name") or "-").strip() not in ("", "-")]
print("carrying a preferred gene name: %d" % len(named))
print()
print("=== 25 most frequent gene names among the category-B losses ===")
for k, v in collections.Counter((r["Preferred_name"] or "").strip()
                                for r in named).most_common(25):
    print("  %-14s %d" % (k, v))

caz = [r for r in B if (r.get("CAZy") or "-").strip() not in ("", "-")]
print()
print("category-B orthogroups with a CAZy assignment: %d" % len(caz))
for k, v in collections.Counter((r["CAZy"] or "").strip() for r in caz).most_common(15):
    print("  %-16s %d" % (k, v))

with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["orthogroup", "COG_category", "Preferred_name", "Description",
                "EC", "KEGG_ko", "CAZy", "PFAMs"])
    for r in sorted(B, key=lambda x: x["query"]):
        w.writerow([r["query"], r.get("COG_category", ""), r.get("Preferred_name", ""),
                    (r.get("Description", "") or "")[:160], r.get("EC", ""),
                    r.get("KEGG_ko", ""), r.get("CAZy", ""), r.get("PFAMs", "")])
print()
print("wrote", OUT)
