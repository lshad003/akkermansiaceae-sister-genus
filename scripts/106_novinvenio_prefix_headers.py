#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Re-emit every novinvenio_set proteome with the source genome prefixed onto each header.
import csv, os, sys
SET = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novinvenio_set"
m = {r["proteome"]: r["source_genome"] for r in csv.DictReader(open(SET + "/genome_map.tsv"), delimiter="\t")}
if len(m) != 14:
    print("REFUSED: genome_map.tsv has %d rows, expected 14" % len(m)); sys.exit(1)
tot = 0
for role in ("ingroup", "outgroup"):
    d = SET + "/" + role + "_prefixed"
    os.makedirs(d, exist_ok=True)
    for f in sorted(os.listdir(SET + "/" + role)):
        if not f.endswith(".faa"): continue
        g = m.get(f)
        if g is None:
            print("REFUSED: %s not in genome_map.tsv" % f); sys.exit(1)
        n = 0
        with open(os.path.join(SET, role, f)) as fi, open(os.path.join(d, f), "w") as fo:
            for line in fi:
                if line.startswith(">"):
                    n += 1
                    fo.write(">" + g + "|" + line[1:].split()[0] + "\n")
                else:
                    fo.write(line)
        tot += n
        print("%-42s %-26s %5d proteins" % (f, g, n))
print("total %d proteins written to ingroup_prefixed and outgroup_prefixed" % tot)
