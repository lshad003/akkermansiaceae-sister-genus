#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# AAI of amphibian-associated Akkermansia to the three Dali et al. clade references.
# Same RBH DIAMOND method as within_genus_aai_v2.py. 30 percent identity floor at reporting.
import csv, os, statistics as st, subprocess, sys, shutil
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
OUT = B + "/results/amph_akk_clade"
os.makedirs(OUT, exist_ok=True)
REFS = {"RS_GCF_000020225.1": ("A. muciniphila", "mammal"),
        "RS_GCF_023516715.1": ("A. massiliensis", "mammal"),
        "RS_GCF_026072915.1": ("A. biwaensis", "mammal"),
        "RS_GCF_001683795.1": ("A. glycaniphila", "reptile"),
        "GB_GCA_019115065.1": ("A. intestinavium", "avian"),
        "GB_GCA_019114365.1": ("A. intestinigallinarum", "avian")}
FAAD = [B + "/results/pangenome/gff199", B + "/results/pangenome/gtdb_akk",
        B + "/results/pangenome/outgroups", B + "/results/pangenome/akkfam",
        B + "/results/dbcan_ehi_nonamph", B + "/results/dbcan_ehi_amphibian"]
idx = {}
for d in FAAD:
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d, f))
D = shutil.which("diamond")
if D is None:
    print("REFUSED: diamond not on PATH, load the module in the job script"); sys.exit(1)
cen = list(csv.DictReader(open(B + "/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"), delimiter="\t"))
amph = sorted(set(r["accession"] for r in cen
                  if r["genus"] == "Akkermansia" and r["host_class"] == "amphibian" and r["accession"] in idx))
if len(amph) != 94:
    print("REFUSED: expected 94 amphibian Akkermansia with proteomes, got %d" % len(amph)); sys.exit(1)
for r in REFS:
    if r not in idx:
        print("REFUSED: no proteome for reference %s" % r); sys.exit(1)
# one representative per amphibian species cluster would be ideal; use all 94, it is cheap enough
def nbest(p):
    b = {}
    for line in open(p):
        f = line.rstrip("\n").split("\t"); pid = float(f[2])
        if f[0] not in b or pid > b[f[0]][1]: b[f[0]] = (f[1], pid)
    return b
DBC = {}
def mkdb(fa, tag):
    if tag not in DBC:
        d = OUT + "/db_" + tag
        subprocess.run([D, "makedb", "--in", fa, "-d", d, "--quiet"], check=True)
        DBC[tag] = d
    return DBC[tag]
def rbh(q, s, qtag, stag):
    db1 = mkdb(s, stag); db2 = mkdb(q, qtag)
    f1, f2 = OUT + "/t_qs.tsv", OUT + "/t_sq.tsv"
    subprocess.run([D, "blastp", "-q", q, "-d", db1, "-o", f1, "--outfmt", "6", "qseqid", "sseqid",
        "pident", "length", "evalue", "--max-target-seqs", "1", "--more-sensitive",
        "--evalue", "1e-5", "--quiet", "--threads", "8"], check=True)
    subprocess.run([D, "blastp", "-q", s, "-d", db2, "-o", f2, "--outfmt", "6", "qseqid", "sseqid",
        "pident", "length", "evalue", "--max-target-seqs", "1", "--more-sensitive",
        "--evalue", "1e-5", "--quiet", "--threads", "8"], check=True)
    fwd, rev = nbest(f1), nbest(f2)
    ids = [p for q0, (s0, p) in fwd.items() if rev.get(s0, (None,))[0] == q0]
    for x in (f1, f2):
        try: os.remove(x)
        except OSError: pass
    return (st.mean(ids), len(ids)) if ids else (float("nan"), 0)
rows = []
for i, g in enumerate(amph, 1):
    for r, (name, clade) in REFS.items():
        a, n = rbh(idx[g], idx[r], g, r)
        rows.append([g, r, name, clade, "%.2f" % a, n])
    if i % 10 == 0: print("  %d of %d done" % (i, len(amph)), flush=True)
with open(OUT + "/amph_akk_vs_dali.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["amphibian_genome", "ref_accession", "ref_species", "dali_clade", "AAI_pct", "n_RBH"])
    w.writerows(rows)
print("\nmean AAI of the 94 amphibian genomes to each reference:")
for r, (name, clade) in REFS.items():
    v = [float(x[4]) for x in rows if x[1] == r and x[4] != "nan"]
    print("  %-26s %-8s n=%2d  min %.2f  median %.2f  max %.2f" % (name, clade, len(v), min(v), st.median(v), max(v)))
for d in DBC.values():
    try: os.remove(d + ".dmnd")
    except OSError: pass
print("\nwrote", OUT + "/amph_akk_vs_dali.tsv")
