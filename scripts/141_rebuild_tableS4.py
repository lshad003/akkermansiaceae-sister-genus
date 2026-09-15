import csv, statistics

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/delim58"
OUT = W + "/TableS4_genus_boundary_metrics.tsv"
TYPE = "CAND_C288"

mem = {r["label"]: r["genome"] for r in csv.DictReader(open(W + "/members.tsv"), delimiter="\t")}
nprot = {}
for r in csv.DictReader(open(W + "/members.tsv"), delimiter="\t"):
    nprot[r["label"]] = sum(1 for l in open(r["path"], errors="replace") if l[:1] == ">")

def mat(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    hdr = rows[0][1:]
    return {(r[0], hdr[j]): float(v) for r in rows[1:] for j, v in enumerate(r[1:])}

AAI = mat(W + "/aai_matrix.tsv")
POCP = mat(W + "/pocp_matrix.tsv")

labs = list(mem)
akk = [x for x in labs if x.startswith("AKK_clade")]
free = [x for x in labs if not x.startswith(("CAND_", "AKK_clade"))]

rows = []
a = [AAI[(TYPE, b)] for b in akk]
p = [POCP[(TYPE, b)] for b in akk]
rows.append(("Akkermansia (24 clade representatives)", "see members.tsv",
             "%.2f" % statistics.median(a), "%.2f-%.2f" % (min(a), max(a)),
             "%.2f" % statistics.median(p), "%.2f-%.2f" % (min(p), max(p)), ""))
for g in sorted(free, key=lambda x: -AAI[(TYPE, x)]):
    rows.append((g, mem[g], "%.2f" % AAI[(TYPE, g)], "",
                 "%.2f" % POCP[(TYPE, g)], "", str(nprot[g])))

with open(OUT, "w") as o:
    o.write("\t".join(["genus", "representative_genome", "AAI_pct_to_type_genome",
                       "AAI_range", "POCP_pct_to_type_genome", "POCP_range",
                       "n_proteins_representative"]) + "\n")
    for r in rows:
        o.write("\t".join(r) + "\n")

print("type genome:", mem[TYPE], "  proteins:", nprot[TYPE])
print("rows:", len(rows))
print()
for r in rows:
    print("  %-40s AAI %-6s %-14s POCP %-6s %s" % (r[0][:40], r[2], r[3], r[4], r[5]))
print()
print("FREE-LIVING AAI RANGE: %.2f to %.2f across %d genera"
      % (min(AAI[(TYPE, g)] for g in free), max(AAI[(TYPE, g)] for g in free), len(free)))
print("WROTE", OUT)

AKKREP = "AKK_clade01"
for lab in labs:
    if lab.startswith("AKK_clade") and mem[lab].replace("GB_","").replace("RS_","") == "GCF_040616545.1":
        AKKREP = lab
with open(OUT) as f: body = f.read().splitlines()
newrow = "\t".join(["Akkermansia (representative RS_GCF_040616545.1)", mem[AKKREP],
                    "%.2f" % AAI[(TYPE, AKKREP)], "", "%.2f" % POCP[(TYPE, AKKREP)], "",
                    str(nprot[AKKREP])])
body.insert(1, newrow)
open(OUT, "w").write("\n".join(body) + "\n")
print()
print("ADDED:", newrow)
