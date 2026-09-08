#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Stage 6 ingroup and 8 outgroup proteomes for the NovInvenio run.
import os, shutil, sys, csv
B   = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
IN  = B + "/results/pangenome/orthofinder_in"
OUT = B + "/results/novinvenio_set"
ING = {"NOVEL_sp03_n16.faa": ("C288", "UHM1073.23039_R.bin.125", "proposed type genome"),
       "NOVEL_sp01_n18.faa": ("C286", "EHM034720", "reduced size and GC cluster"),
       "NOVEL_sp02_n16.faa": ("C287", "EHM034692", "large cluster"),
       "NOVEL_sp04_n15.faa": ("C289", "EHM034674", "large cluster"),
       "NOVEL_sp09_n4.faa":  ("C294", "UHM1210.23070_R.bin.101", "below 65 percent AAI to type"),
       "NOVEL_sp13_n1.faa":  ("C298", "STP248.12601_R.bin.118", "deepest cluster, 95.9 percent complete")}
akk = sorted(f for f in os.listdir(IN) if f.startswith("AKK_clade"))[:8]
if len(akk) != 8:
    print("REFUSED: found %d AKK proteomes, expected at least 8" % len(akk)); sys.exit(1)
miss = [f for f in ING if not os.path.exists(os.path.join(IN, f))]
if miss:
    print("REFUSED: missing %s" % miss); sys.exit(1)
for d in ("ingroup", "outgroup"):
    os.makedirs(os.path.join(OUT, d), exist_ok=True)
rows = []
for f, (c, g, why) in sorted(ING.items()):
    shutil.copy(os.path.join(IN, f), os.path.join(OUT, "ingroup", f))
    n = sum(1 for l in open(os.path.join(IN, f)) if l.startswith(">"))
    rows.append([f, "ingroup", "candidate_genus", c, g, n, why])
for f in akk:
    shutil.copy(os.path.join(IN, f), os.path.join(OUT, "outgroup", f))
    n = sum(1 for l in open(os.path.join(IN, f)) if l.startswith(">"))
    rows.append([f, "outgroup", "Akkermansia", f.split("_")[1], "clade representative", n, "GTDB and MAG Akkermansia"])
with open(OUT + "/groups.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["proteome", "role", "group", "cluster_or_clade", "source_genome", "n_proteins", "note"])
    w.writerows(rows)
print("staged %d ingroup and %d outgroup proteomes in %s" % (len(ING), len(akk), OUT))
os.system("column -t " + OUT + "/groups.tsv")
