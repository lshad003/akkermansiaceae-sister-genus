from ete3 import Tree
import csv

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TREE = CH3 + "/results/novel_akk_tree/akk_placement.nwk"
W = CH3 + "/results/delim58"

rows = list(csv.DictReader(open(W + "/members.tsv"), delimiter="\t"))
def norm(x): return x.replace("GB_", "").replace("RS_", "").strip()
wanted = {norm(r["genome"]): r["label"] for r in rows}
print("wanted:", len(wanted))

t = Tree(TREE, format=1)
og = [l for l in t.iter_leaves() if l.name.split("|")[0].startswith("OUTGROUP")]
anc = t.get_common_ancestor(og)
if anc is t:
    print("REFUSED: outgroup MRCA is root"); raise SystemExit(1)
t.set_outgroup(anc)
nov = [l for l in t.iter_leaves() if l.name.split("|")[0] == "NOVEL"]
akk = [l for l in t.iter_leaves() if l.name.split("|")[0] == "Akkermansia"]
print("NOVEL+Akk clade after rooting:", len(t.get_common_ancestor(nov+akk).get_leaves()), "(expect 292)")

sup = {}
for tag, tips in (("CAND", nov), ("AKK", akk), ("SHARED", nov + akk)):
    nd = t.get_common_ancestor(tips)
    sup[tag] = (nd.name or "").strip() or "NA"
with open(W + "/node_support.tsv", "w") as o:
    o.write("node\tsupport\n")
    for k, v in sup.items():
        o.write("%s\t%s\n" % (k, v))
print("support:", sup)

matched = {}
for leaf in t.iter_leaves():
    p = leaf.name.split("|")
    acc = norm(p[4]) if len(p) >= 5 else norm(p[-1])
    if acc in wanted and acc not in matched:
        matched[acc] = leaf.name
missing = sorted(set(wanted) - set(matched))
if missing:
    print("REFUSED: %d missing" % len(missing))
    for x in missing[:8]: print("  ", x)
    raise SystemExit(1)

t.prune(list(matched.values()), preserve_branch_length=True)
rev = {v: wanted[k] for k, v in matched.items()}
for leaf in t.iter_leaves():
    leaf.name = rev[leaf.name]
t.ladderize()
t.write(outfile=W + "/delim58_tree.nwk", format=1)

order = [l.name for l in t.iter_leaves()]
with open(W + "/delim58_tree_order.tsv", "w") as o:
    o.write("order\tlabel\n")
    for i, lab in enumerate(order, 1):
        o.write("%d\t%s\n" % (i, lab))

def block(pref):
    ix = [i for i, x in enumerate(order) if x.startswith(pref)]
    return ix[0]+1, ix[-1]+1, ix == list(range(ix[0], ix[-1]+1))
print("tips:", len(order))
print("CAND block %d-%d contiguous %s" % block("CAND_"))
print("AKK  block %d-%d contiguous %s" % block("AKK_clade"))
print()
for lab in order: print("  ", lab)
