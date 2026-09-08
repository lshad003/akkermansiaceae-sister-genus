#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Pick one GTDB genome per Dali et al. Akkermansia reference species; report proteome availability.
import os, csv
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
WANT = ["Akkermansia muciniphila", "Akkermansia massiliensis", "Akkermansia biwaensis",
        "Akkermansia glycaniphila", "Akkermansia intestinavium", "Akkermansia intestinigallinarum"]
CLADE = {"Akkermansia muciniphila": "mammal", "Akkermansia massiliensis": "mammal",
         "Akkermansia biwaensis": "mammal", "Akkermansia glycaniphila": "reptile",
         "Akkermansia intestinavium": "avian", "Akkermansia intestinigallinarum": "avian"}
hits = {}
for line in open(B + "/data/bac120_taxonomy_r226.tsv"):
    acc, tax = line.rstrip("\n").split("\t")[:2]
    if ";g__Akkermansia;" not in tax: continue
    sp = tax.split(";s__")[-1].strip()
    if sp in WANT:
        hits.setdefault(sp, []).append(acc)
FAAD = [B + "/results/pangenome/gff199", B + "/results/pangenome/gtdb_akk",
        B + "/results/pangenome/outgroups", B + "/results/pangenome/akkfam",
        B + "/results/dbcan_ehi_nonamph", B + "/results/dbcan_ehi_amphibian"]
have = {}
for d in FAAD:
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): have[f[:-4]] = os.path.join(d, f)
print("%-34s %-8s %6s  %s" % ("species", "clade", "n_GTDB", "proteome we already have"))
need = []
for sp in WANT:
    accs = hits.get(sp, [])
    got = [a for a in accs if a in have]
    print("%-34s %-8s %6d  %s" % (sp, CLADE[sp], len(accs), got[0] if got else "NONE -> need " + (accs[0] if accs else "NO GTDB ENTRY")))
    if not got and accs: need.append((sp, accs[0]))
if need:
    print("\nto download:")
    for sp, a in need: print("   %-34s %s" % (sp, a))
