#!/usr/bin/env python3
# Amino acid identity of every candidate genome to the type genome
# Source: ch3-chitin-evolution/scripts/within_genus_aai.py
# Output: results/within_genus_aai/within_genus_aai_to_type.tsv
import csv, os, statistics, subprocess, sys, shutil

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
P = f"{BASE}/results/pangenome"
FAAD = [f"{P}/gff199", f"{P}/gtdb_akk", f"{P}/outgroups", f"{P}/akkfam",
        f"{BASE}/results/dbcan_ehi_nonamph", f"{BASE}/results/dbcan_ehi_amphibian"]
OUT  = f"{BASE}/results/within_genus_aai"
TYPE = "EHM058340"

D = shutil.which("diamond")
if D is None:
    print("REFUSED: diamond not on PATH. Load the module in the job script.")
    sys.exit(1)

idx = {}
for d in FAAD:
    if not os.path.isdir(d):
        continue
    for f in os.listdir(d):
        if f.endswith(".faa"):
            idx[f[:-4]] = os.path.join(d, f)

rows = [r for r in csv.DictReader(open(CEN), delimiter="\t")
        if r["family"] == "Akkermansiaceae"
        and (not r["genus"].strip() or r["genus"] in ("unknown", "NO_GENUS"))
        and r["host_class"] == "amphibian"]

genomes = sorted(set(r["accession"] for r in rows))
if len(genomes) != 105:
    print("REFUSED: expected 105 amphibian novel genomes, got %d" % len(genomes))
    sys.exit(1)

missing = [g for g in genomes if g not in idx]
if TYPE not in idx:
    print("REFUSED: type genome %s has no proteome in the dbCAN directories" % TYPE)
    sys.exit(1)
print("genomes: %d, with proteome: %d, missing: %d" % (
    len(genomes), len(genomes) - len(missing), len(missing)))
if missing:
    print("REFUSED: %d of 105 novel genomes have no proteome. A within-genus" % len(missing))
    print("distribution computed on a subset is not the denominator the exclusion needs.")
    for m in missing:
        print("   ", m)
    sys.exit(1)

def nbest(path):
    best = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        pid = float(f[2])
        if f[0] not in best or pid > best[f[0]][1]:
            best[f[0]] = (f[1], pid)
    return best

def rbh_aai(qa, sa):
    q, s = idx[qa], idx[sa]
    db1 = f"{OUT}/tmp_s"; db2 = f"{OUT}/tmp_q"
    f1  = f"{OUT}/tmp_qs.tsv"; f2 = f"{OUT}/tmp_sq.tsv"
    subprocess.run([D, "makedb", "--in", s, "-d", db1, "--quiet"], check=True)
    subprocess.run([D, "blastp", "-q", q, "-d", db1, "-o", f1, "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "evalue", "--max-target-seqs", "1",
        "--more-sensitive", "--evalue", "1e-5", "--quiet", "--threads", "8"], check=True)
    subprocess.run([D, "makedb", "--in", q, "-d", db2, "--quiet"], check=True)
    subprocess.run([D, "blastp", "-q", s, "-d", db2, "-o", f2, "--outfmt", "6",
        "qseqid", "sseqid", "pident", "length", "evalue", "--max-target-seqs", "1",
        "--more-sensitive", "--evalue", "1e-5", "--quiet", "--threads", "8"], check=True)
    fwd = nbest(f1); rev = nbest(f2)
    ids = [pid for q0, (s0, pid) in fwd.items() if rev.get(s0, (None,))[0] == q0]
    for x in (f1, f2, db1 + ".dmnd", db2 + ".dmnd"):
        try: os.remove(x)
        except OSError: pass
    if not ids:
        return float("nan"), 0
    return statistics.mean(ids), len(ids)

res = []
for i, g in enumerate(genomes, 1):
    if g == TYPE or g not in idx:
        continue
    a, n = rbh_aai(g, TYPE)
    res.append((g, a, n))
    print("%4d/%d  %-34s AAI %.2f  RBH %d" % (i, len(genomes), g, a, n), flush=True)

with open(f"{OUT}/within_genus_aai_to_type.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome", "reference_type_genome", "AAI_pct", "n_RBH"])
    for g, a, n in res:
        w.writerow([g, TYPE, "%.2f" % a, n])

vals = sorted(a for _, a, _ in res)
print()
print("WITHIN-GENUS AAI to type genome %s, n=%d" % (TYPE, len(vals)))
print("  min %.2f  Q1 %.2f  median %.2f  Q3 %.2f  max %.2f" % (
    vals[0], vals[len(vals)//4], statistics.median(vals),
    vals[3*len(vals)//4], vals[-1]))
print()
print("TORTOISE for comparison: 53.26 and 53.28")
print("(results/tortoise/tortoise_exclusion_aai.tsv)")
print("wrote", f"{OUT}/within_genus_aai_to_type.tsv")
