#!/usr/bin/env python3
# Per-sample abundance of the genus assembled from the amplicon survey
# Source: ch3-chitin-evolution/scripts/rurik_per_sample_table.py
# Output: results/rurik_16s/per_sample_abundance.tsv
import csv, sys, openpyxl
CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
OUT = CH3 + "results/rurik_16s/per_sample_abundance.tsv"
csv.field_size_limit(sys.maxsize)

wb = openpyxl.load_workbook(CH3+"data/251027_16S_metadata_final.xlsx", read_only=True, data_only=True)
ws = wb.worksheets[0]
it = ws.iter_rows(values_only=True)
hdr = ['' if c is None else str(c) for c in next(it)]
I = {h:i for i,h in enumerate(hdr)}
meta = {}
for row in it:
    c = ['' if x is None else str(x).strip() for x in row]
    if I["sample_name"] < len(c):
        meta[c[I["sample_name"]]] = (
            c[I["Clade_Order"]] if I["Clade_Order"] < len(c) else "",
            c[I["env_broad_scale"]] if I["env_broad_scale"] < len(c) else "",
            (c[I["host_genus"]]+" "+c[I["host_species"]]).strip() if I["host_species"] < len(c) else "")
wb.close()

want = set(l.strip() for l in open(CH3+"results/rurik_16s/matched_ids.txt") if l.strip())
with open(CH3+"data/251001_16S_count_table.csv", newline="") as fh:
    r = csv.reader(fh)
    samples = next(r)[1:]
    depth = [0.0]*len(samples)
    ours  = [0.0]*len(samples)
    for row in r:
        mine = row[0].strip().strip('"') in want
        for i, v in enumerate(row[1:]):
            try: x = float(v)
            except ValueError: continue
            depth[i] += x
            if mine: ours[i] += x

n = 0
with open(OUT, "w") as o:
    o.write("sample\thost_order\tenv\thost_species\tdepth\tgenus_reads\tpercent\n")
    for i, s in enumerate(samples):
        m = meta.get(s)
        if not m or m[0] not in ("Caudata","Anura","Squamata","Testudines") or depth[i] <= 0:
            continue
        env = "captive" if m[1].lower() in ("zoo","captive","lab population") else "wild"
        o.write("%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.6f\n" % (s, m[0], env, m[2], depth[i], ours[i], 100.0*ours[i]/depth[i]))
        n += 1
print("rows written:", n)
print("WROTE:", OUT)
