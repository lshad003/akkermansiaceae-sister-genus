import csv, os, glob

REPO = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus"
CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
T1 = REPO + "/tables/TableS1_genome_quality.tsv"
FA = CH3 + "/data/amphibia_gtdbtk_input"

rows = list(csv.DictReader(open(T1), delimiter="\t"))
if "N50_bp" in rows[0]:
    print("REFUSED: N50_bp already present"); raise SystemExit

def n50(p):
    L = []
    n = 0
    for line in open(p, errors="replace"):
        if line[:1] == ">":
            if n: L.append(n)
            n = 0
        else:
            n += len(line.strip())
    if n: L.append(n)
    L.sort(reverse=True)
    tot = sum(L); run = 0
    for x in L:
        run += x
        if run >= tot / 2: return x, len(L), tot
    return 0, len(L), tot

miss = []
for r in rows:
    p = None
    for e in (".fa", ".fna", ".fasta"):
        c = os.path.join(FA, r["genome"] + e)
        if os.path.exists(c): p = c; break
    if not p:
        miss.append(r["genome"]); r["N50_bp"] = "NA"; continue
    v, nc, tot = n50(p)
    r["N50_bp"] = str(v)
    if r["genome"] == "UHM1073.23039_R.bin.125":
        print("TYPE GENOME: N50 %d, contigs %d (table says %s), total %.2f Mb (table %s)"
              % (v, nc, r["contigs"], tot / 1e6, r["size_Mb"]))

if miss:
    print("no assembly for %d genomes, e.g. %s" % (len(miss), miss[:3]))

cols = list(rows[0].keys())
with open(T1, "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
print("WROTE", T1, "rows:", len(rows))
