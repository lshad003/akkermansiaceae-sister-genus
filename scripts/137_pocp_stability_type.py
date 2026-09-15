import os, csv, subprocess, shutil

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
OUT = CH3 + "/results/aai_pocp_type"
W = OUT + "/stability_tmp"
os.makedirs(W, exist_ok=True)
REF = "RS_GCF_040616545.1"
SIX = ["EHM060011", "UHM979.41089_R.bin.103", "UHM967.23060_R.bin.185",
       "EHM034674", "EHM059669", "UHM1073.23039_R.bin.125"]

D = shutil.which("diamond")
if not D:
    print("REFUSED: diamond not on PATH, run: module load diamond/2.2.6"); raise SystemExit(1)
v = subprocess.run([D, "version"], capture_output=True, text=True).stdout.strip()
print(v)
if "2.2.6" not in v:
    print("REFUSED: need DIAMOND 2.2.6"); raise SystemExit(1)

DIRS = [CH3 + "/results/pangenome/" + d for d in ("gff199", "akkfam", "outgroups", "gtdb_akk")]
DIRS += [CH3 + "/results/" + d for d in ("dbcan_scaffold", "dbcan_ehi_amphibian", "dbcan_ehi_nonamph")]
idx = {}
for d in DIRS:
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d, f))

def find(a):
    return idx.get(a) or idx.get(a.replace("GB_", "").replace("RS_", ""))

def nprot(p):
    return sum(1 for l in open(p, errors="replace") if l[:1] == ">")

def conserved(q, s, tag):
    db = "%s/db_%s" % (W, tag)
    o = "%s/h_%s.tsv" % (W, tag)
    subprocess.run([D, "makedb", "--in", s, "-d", db, "--quiet"], check=True)
    subprocess.run([D, "blastp", "-q", q, "-d", db, "-o", o, "--outfmt", "6",
                    "qseqid", "sseqid", "pident", "length", "qlen", "evalue",
                    "--max-target-seqs", "1", "--more-sensitive", "--evalue", "1e-5",
                    "--quiet", "--threads", "8"], check=True)
    n = 0
    for line in open(o):
        f = line.split("\t")
        if float(f[2]) > 40.0 and int(f[4]) > 0 and (int(f[3]) / int(f[4])) > 0.5:
            n += 1
    os.remove(o); os.remove(db + ".dmnd")
    return n

rp = find(REF)
if not rp:
    print("REFUSED: no proteome for", REF); raise SystemExit(1)
rows = []
for g in SIX:
    gp = find(g)
    if not gp:
        print("REFUSED: no proteome for", g); raise SystemExit(1)
    c1 = conserved(gp, rp, "fwd"); c2 = conserved(rp, gp, "rev")
    t1 = nprot(gp); t2 = nprot(rp)
    p = 100.0 * (c1 + c2) / (t1 + t2)
    rows.append((g, c1, c2, t1, t2, p))
    print("  %-26s POCP %6.2f" % (g, p))

with open(OUT + "/pocp_stability_type.tsv", "w") as o:
    o.write("candidate_genome\tC1\tC2\tT1\tT2\tPOCP_pct\n")
    for r in rows:
        o.write("%s\t%d\t%d\t%d\t%d\t%.2f\n" % r)
v = [r[5] for r in rows]
print()
print("RANGE %.2f to %.2f   (old EHM058340-anchored set: 46.60 to 54.53)" % (min(v), max(v)))
print("WROTE", OUT + "/pocp_stability_type.tsv")
