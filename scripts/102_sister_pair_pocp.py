#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Which of the 28 calibrated genus pairs are sisters in the bac120 tree?
import csv, re, sys
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
NWK = B + "/results/novel_akk_tree/akk_placement.nwk"
CAL = B + "/results/aai_pocp/pocp_calibration.tsv"
s = open(NWK).read().strip().rstrip(";")

def parse(t):
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
            st = pos[0]
            while pos[0] < len(t) and t[pos[0]] not in "(),": pos[0] += 1
            return kids
        st = pos[0]
        while pos[0] < len(t) and t[pos[0]] not in "(),": pos[0] += 1
        return t[st:pos[0]].split(":")[0]
    return node()

root = parse(s)
def tips(n):
    if isinstance(n, str): return [n]
    out = []
    for k in n: out += tips(k)
    return out
def genus(tip): return tip.split("|")[0]

all_tips = tips(root)
if len(all_tips) != 784:
    print("REFUSED: parsed %d tips, expected 784" % len(all_tips)); sys.exit(1)

# for every internal node record the set of genera below it
clades = []
def walk(n):
    if isinstance(n, str): return {genus(n)}
    g = set()
    for k in n: g |= walk(k)
    clades.append(g)
    return g
walk(root)

def is_sister(a, b):
    # sisters if some clade contains exactly {a, b} and nothing else
    return any(c == {a, b} for c in clades)

cal = list(csv.DictReader(open(CAL), delimiter="\t"))
print("%-18s %-18s %8s  %s" % ("genusA", "genusB", "POCP", "sister_in_tree"))
sis, non = [], []
for r in cal:
    a, b, p = r["genusA"], r["genusB"], float(r["POCP_pct"])
    ok = is_sister(a, b)
    (sis if ok else non).append(p)
    print("%-18s %-18s %8.2f  %s" % (a, b, p, "YES" if ok else "no"))
print()
print("sister pairs     n=%d  %s" % (len(sis), ["%.2f" % x for x in sorted(sis)]))
print("non-sister pairs n=%d  min %.2f  max %.2f" % (len(non), min(non), max(non)) if non else "")
print()
print("NOVEL and Akkermansia sisters in this tree:", is_sister("NOVEL", "Akkermansia"))
