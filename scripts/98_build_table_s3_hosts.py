#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# TableS3: one row per animal for the 105 candidate genomes.
# Animal unit: EHI_2025 -> ncbi_biosample ; herptile_MAG -> sample_id prefix before first dot
# 46 EHI biosamples + 20 herptile animals = 66 (verify_animal_count.py rule)
import csv, os, sys
from collections import defaultdict
BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN  = BASE + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
EHI  = BASE + "/data/ehi_2025_annotated/ehi_mags_annotated_v2.tsv"
HERP = BASE + "/data/herptile_gtdbtk_compiled.tsv"
OUT  = BASE + "/results/tables/TableS3_host_metadata.tsv"
os.makedirs(os.path.dirname(OUT), exist_ok=True)

cen = list(csv.DictReader(open(CEN), delimiter="\t"))
cand = [r for r in cen if r["family"] == "Akkermansiaceae"
        and (not r["genus"].strip() or r["genus"] in ("unknown", "NO_GENUS"))
        and r["host_class"] == "amphibian"]
if len(cand) != 105:
    print("REFUSED: expected 105 candidate genomes, got %d" % len(cand)); sys.exit(1)

ehi = {r["accession"]: r for r in csv.DictReader(open(EHI), delimiter="\t")}
herp = {}
for r in csv.DictReader(open(HERP), delimiter="\t"):
    herp[r["user_genome"]] = r

animals = defaultdict(lambda: defaultdict(set))
counts = defaultdict(int)
for r in cand:
    acc, ds = r["accession"], r["from_dataset"]
    if ds == "EHI_2025":
        e = ehi.get(acc)
        if e is None:
            print("REFUSED: %s not in EHI table" % acc); sys.exit(1)
        aid = e["ncbi_biosample"].strip()
        rec = {"host_species": e["host_animal_type"], "host_class": e["host_class"],
               "captivity_status": e["captivity_status"], "host_animal_type": e["host_label"],
               "sample_id": e["sample_id"]}
    else:
        aid = r["sample_id"].strip().split(".")[0]
        rec = {"host_species": "", "host_class": r["host_class"],
               "captivity_status": r["captivity_status"], "host_animal_type": r["host_animal_type"],
               "sample_id": r["sample_id"]}
    if not aid:
        print("REFUSED: empty animal id for %s (%s)" % (acc, ds)); sys.exit(1)
    key = (ds, aid)
    counts[key] += 1
    for k, v in rec.items():
        if v: animals[key][k].add(v)
    animals[key]["genomes"].add(acc)

if len(counts) != 66:
    print("REFUSED: expected 66 animals, got %d" % len(counts)); sys.exit(1)

def one(key, field):
    v = sorted(animals[key][field])
    return v[0] if len(v) == 1 else ("|".join(v) if v else "NA")

with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["animal_id", "dataset", "host_species", "host_animal_type", "host_class",
                "captivity_status", "n_candidate_genomes", "sample_id", "genomes"])
    for (ds, aid) in sorted(counts, key=lambda k: (k[0], k[1])):
        key = (ds, aid)
        w.writerow([aid, ds, one(key, "host_species"), one(key, "host_animal_type"),
                    one(key, "host_class"), one(key, "captivity_status"), counts[key],
                    one(key, "sample_id"), ";".join(sorted(animals[key]["genomes"]))])

from collections import Counter
print("animals: %d  genomes: %d" % (len(counts), sum(counts.values())))
print("by dataset:", Counter(ds for ds, a in counts))
print("captivity:", Counter(one((ds, a), "captivity_status") for ds, a in counts))
print("genomes per animal, top 8:", sorted(counts.items(), key=lambda kv: -kv[1])[:8])
print("wrote", OUT)

# ---- enrich herptile rows from the 16S metadata sheet, keyed on animal_number ----
import pandas as pd
XL = BASE + "/data/251027_16S_metadata_final.xlsx"
d = pd.read_excel(XL, sheet_name="251027_16S_metadata_final")
meta = {}
for _, r in d.iterrows():
    a = str(r["animal_number"]).strip()
    if a and a != "nan" and a not in meta:
        meta[a] = r
rows = list(csv.DictReader(open(OUT), delimiter="\t"))
cols = ["animal_id", "dataset", "host_species", "host_order", "diet", "host_animal_type",
        "host_class", "captivity_status", "n_candidate_genomes", "sample_id", "genomes"]
nmiss = 0
for x in rows:
    if x["dataset"] == "herptile_MAG":
        m = meta.get(x["animal_id"])
        if m is None:
            nmiss += 1; x["host_order"] = "NA"; x["diet"] = "NA"; continue
        x["host_species"] = "%s %s" % (m["host_genus"], m["host_species"])
        x["host_order"] = str(m["Clade_Order"])
        x["diet"] = str(m["Diet"])
    else:
        x["host_order"] = "NA"; x["diet"] = "NA"
if nmiss:
    print("REFUSED: %d herptile animals absent from the metadata sheet" % nmiss); sys.exit(1)
with open(OUT, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    w.writeheader(); w.writerows(rows)
print("enriched %d herptile rows from %s" % (sum(1 for x in rows if x["dataset"] == "herptile_MAG"), XL))
print("EHI host_order and diet left NA: not in that sheet")
