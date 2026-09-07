#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Do the five clusters below 65% AAI to the type genome form a clade in the bac120 tree?
# LOW set defined in scripts/compute_cluster_aai.py line 12.
import csv, sys
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
NWK = B + "/results/novel_akk_tree/akk_placement.nwk"
S1  = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus/tables/TableS1_genome_quality.tsv"
LOW = {"C286", "C294", "C296", "C298", "C300"}

cl = {}
for r in csv.DictReader(open(S1), delimiter="\t"):
    cl[r["genome"]] = r["species_cluster"]
if len(cl) != 105:
    print("REFUSED: TableS1 has %d genomes, expected 105" % len(cl)); sys.exit(1)

t = open(NWK).read().strip().rstrip(";")
pos = [0]
def node():
    if t[pos[0]] == "(":
        pos[0] += 1
        kids = [node()]
        while t[pos[0]] == ",":
            pos[0] += 1; kids.append(node())
        if t[pos[0]] != ")":
            print("REFUSED: malformed newick at %d" % pos[0]); sys.exit(1)
        pos[0] += 1
        while pos[0] < len(t) and t[pos[0]] not in "(),": pos[0] += 1
        return kids
    st = pos[0]
    while pos[0] < len(t) and t[pos[0]] not in "(),": pos[0] += 1
    return t[st:pos[0]].split(":")[0]
root = node()

clades = []
def walk(n):
    if isinstance(n, str): return {n}
    s = set()
    for k in n: s |= walk(k)
    clades.append(s)
    return s
alltips = walk(root)
if len(alltips) != 784:
    print("REFUSED: parsed %d tips, expected 784" % len(alltips)); sys.exit(1)

nov = {x for x in alltips if x.split("|")[0] == "NOVEL"}
if len(nov) != 105:
    print("REFUSED: %d NOVEL tips, expected 105" % len(nov)); sys.exit(1)
def gname(tip): return tip.split("|")[-1]
miss = [x for x in nov if gname(x) not in cl]
if miss:
    print("REFUSED: %d NOVEL tips not in TableS1, e.g. %s" % (len(miss), miss[:3])); sys.exit(1)

lowtips = {x for x in nov if cl[gname(x)] in LOW}
print("LOW clusters %s -> %d genomes" % (sorted(LOW), len(lowtips)))
print("LOW set is a clade:", lowtips in clades)
for c in sorted(LOW):
    s = {x for x in nov if cl[gname(x)] == c}
    print("  %-6s n=%2d  monophyletic: %s" % (c, len(s), s in clades))
sm = [c for c in clades if lowtips <= c]
if sm:
    small = min(sm, key=len)
    extra = sorted({cl[gname(x)] for x in small if x in nov and cl[gname(x)] not in LOW})
    print("smallest clade containing all LOW tips: %d tips, also contains clusters %s" % (len(small), extra))
