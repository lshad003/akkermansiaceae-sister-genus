#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Rename the type-genome column and add two derived columns to TableS1.
# MIMAG HQ = >90% complete, <5% contamination, 5S/16S/23S present, >=18 tRNA types,
# and the 16S confirmed to belong to the genome (results/rrna16s/hq_16s_verdict.txt).
import csv, sys
T = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus/tables/TableS1_genome_quality.tsv"
VALID   = {"EHM058340", "EHM060011", "UHM973.23044_R.bin.53"}
INVALID = {"EHM059280", "EHM059433"}

with open(T, newline="") as f:
    rd = csv.DictReader(f, delimiter="\t")
    cols = list(rd.fieldnames); rows = list(rd)
if len(rows) != 105:
    print("REFUSED: expected 105 rows, got %d" % len(rows)); sys.exit(1)
if cols.count("meets_all_criteria") != 1:
    print("REFUSED: meets_all_criteria found %d times" % cols.count("meets_all_criteria")); sys.exit(1)
for c in ("rRNA_16S_validated", "MIMAG_high_quality"):
    if c in cols:
        print("REFUSED: %s already exists, nothing written" % c); sys.exit(1)

cols[cols.index("meets_all_criteria")] = "meets_type_genome_criteria"
for r in rows:
    r["meets_type_genome_criteria"] = r.pop("meets_all_criteria")

def base_ok(r):
    return (float(r["completeness"]) > 90 and float(r["contamination"]) < 5
            and r["rRNA_5S"].strip() == "+" and r["rRNA_16S"].strip() == "+"
            and r["rRNA_23S"].strip() == "+" and int(r["tRNA_types"]) >= 18)

need = [r["genome"] for r in rows if base_ok(r)]
nover = [g for g in need if g not in VALID and g not in INVALID]
if nover:
    print("REFUSED: these pass all other MIMAG criteria but have no 16S containment verdict: %s" % nover)
    sys.exit(1)

for r in rows:
    g = r["genome"]
    if g in VALID:      r["rRNA_16S_validated"] = "yes"
    elif g in INVALID:  r["rRNA_16S_validated"] = "no"
    elif r["rRNA_16S"].strip() == "+": r["rRNA_16S_validated"] = "not_assessed"
    else:               r["rRNA_16S_validated"] = "not_applicable"
    r["MIMAG_high_quality"] = "yes" if (base_ok(r) and r["rRNA_16S_validated"] == "yes") else "no"

n = sum(1 for r in rows if r["MIMAG_high_quality"] == "yes")
if n != 3:
    print("REFUSED: MIMAG_high_quality = %d, expected 3. Nothing written." % n); sys.exit(1)

cols += ["rRNA_16S_validated", "MIMAG_high_quality"]
with open(T, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows)
print("passed all other MIMAG criteria: %d  ->  %s" % (len(need), need))
print("MIMAG_high_quality yes:", [r["genome"] for r in rows if r["MIMAG_high_quality"] == "yes"])
print("rewrote", T)
