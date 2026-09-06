#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# TableS4: genus-boundary metrics for the 13 neighbouring Akkermansiaceae genera.
# AAI anchor  = UHM1073.23039_R.bin.125 (current type genome)
# POCP anchor = EHM058340 (2,415 proteins), the anchor used when POCP was computed
import csv, sys
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
AAI  = B + "/results/within_genus_aai/aai_vs_type_v2_neighbours.tsv"
POCP = B + "/results/aai_pocp/pocp_fixed.tsv"
OUT  = B + "/results/tables/TableS4_genus_boundary_metrics.tsv"

a = {r["genome"]: r for r in csv.DictReader(open(AAI), delimiter="\t")}
p = {r["rep"]: r for r in csv.DictReader(open(POCP), delimiter="\t")}
if len(a) != 13 or len(p) != 13:
    print("REFUSED: expected 13 rows in each file, got AAI %d POCP %d" % (len(a), len(p))); sys.exit(1)
if set(a) != set(p):
    print("REFUSED: representative genomes differ between the two files")
    print("  AAI only :", sorted(set(a) - set(p)))
    print("  POCP only:", sorted(set(p) - set(a))); sys.exit(1)
anch = set(r["reference_type_genome"] for r in a.values())
if anch != {"UHM1073.23039_R.bin.125"}:
    print("REFUSED: unexpected AAI anchor %s" % anch); sys.exit(1)
t1 = set(r["T1"] for r in p.values())
if t1 != {"2415"}:
    print("REFUSED: POCP anchor proteome size is not constant: %s" % t1); sys.exit(1)

rows = []
for g in sorted(a, key=lambda k: -float(a[k]["AAI_pct"])):
    ar, pr = a[g], p[g]
    if ar["label"] != pr["genus"]:
        print("REFUSED: genus label mismatch for %s: %s vs %s" % (g, ar["label"], pr["genus"])); sys.exit(1)
    rows.append([pr["genus"], g, ar["AAI_pct"], ar["n_RBH"], ar["AAI_pct_old_EHM058340"],
                 pr["POCP_pct"], pr["T2"]])
with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["genus", "representative_genome", "AAI_pct_to_type_genome", "n_reciprocal_best_hits",
                "AAI_pct_to_EHM058340", "POCP_pct_to_EHM058340", "n_proteins_representative"])
    w.writerows(rows)
print("wrote %s  (%d genera)" % (OUT, len(rows)))
for r in rows:
    print("%-18s AAI %6s  POCP %6s  proteins %5s" % (r[0], r[2], r[5], r[6]))
