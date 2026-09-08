#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Map every novinvenio_set proteome to its source genome and stage the nucleotide FASTA.
import csv, os, shutil, sys
B    = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SET  = B + "/results/novinvenio_set"
FAAD = [B + "/results/pangenome/gff199", B + "/results/pangenome/gtdb_akk",
        B + "/results/pangenome/outgroups", B + "/results/pangenome/akkfam",
        B + "/results/dbcan_ehi_nonamph", B + "/results/dbcan_ehi_amphibian"]
FNAD = [B + "/data/amphibia_gtdbtk_input", B + "/results/novel_akk_ani/ref_fna",
        B + "/data/ehi_genomes", B + "/data/gtdb_akk_fna"]

# index: first protein header -> source genome name
idx = {}
for d in FAAD:
    if not os.path.isdir(d): continue
    for f in os.listdir(d):
        if not f.endswith(".faa"): continue
        with open(os.path.join(d, f)) as fh:
            for line in fh:
                if line.startswith(">"):
                    idx.setdefault(line[1:].split()[0], []).append(f[:-4]); break

# index: genome name -> nucleotide fasta
fna = {}
for d in FNAD:
    if not os.path.isdir(d):
        print("note: no such directory, skipped: " + d); continue
    for f in os.listdir(d):
        for e in (".fa", ".fna", ".fasta"):
            if f.endswith(e):
                fna.setdefault(f[:-len(e)], os.path.join(d, f)); break

os.makedirs(SET + "/genomes", exist_ok=True)
rows, missing = [], []
for role in ("ingroup", "outgroup"):
    for f in sorted(os.listdir(SET + "/" + role)):
        if not f.endswith(".faa"): continue
        h = open(os.path.join(SET, role, f)).readline()[1:].split()[0]
        g = sorted(set(idx.get(h, [])))
        if len(g) != 1:
            print("REFUSED: header %s from %s matched %d distinct genomes: %s" % (h, f, len(g), g))
            sys.exit(1)
        g = g[0]
        p = fna.get(g)
        if p:
            shutil.copy(p, SET + "/genomes/" + g + ".fna")
        else:
            missing.append((f, g))
        rows.append([f, role, g, h, p or "MISSING"])

with open(SET + "/genome_map.tsv", "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["proteome", "role", "source_genome", "first_protein_header", "nucleotide_fasta"])
    w.writerows(rows)
for r in rows:
    print("%-42s %-9s %-26s %s" % (r[0], r[1], r[2], "OK" if r[4] != "MISSING" else "MISSING FASTA"))
print()
print("staged %d genomes in %s/genomes" % (len(rows) - len(missing), SET))
if missing:
    print("MISSING nucleotide FASTA for %d:" % len(missing))
    for f, g in missing: print("   %s -> %s" % (f, g))
