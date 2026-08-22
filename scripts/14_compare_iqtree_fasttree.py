#!/usr/bin/env python3
# Topology compared between the two inference methods
# Source: ch3-chitin-evolution/scripts/compare_iqtree_fasttree.py
# Output: stdout
from ete3 import Tree
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_tree"
IQ=f"{B}/iqtree_placement/akk_placement_iqtree.treefile"
FT=f"{B}/akk_placement.nwk"

def load(p):
    for fmt in (1,0,5,3):
        try: return Tree(p, format=fmt)
        except Exception: continue
    return None

def grp(n): return n.split("|")[0]

for lab,path in (("IQ-TREE",IQ),("FastTree",FT)):
    t=load(path)
    if t is None: print(lab,"PARSE FAILED",path); continue
    leaves=t.get_leaves()
    outg=[l.name for l in leaves if grp(l.name).startswith("OUTGROUP")]
    t.set_outgroup(t.get_common_ancestor(outg))
    nov=[l.name for l in leaves if grp(l.name)=="NOVEL"]
    akk=[l.name for l in leaves if grp(l.name)=="Akkermansia"]
    print(f"\n=== {lab} ({len(leaves)} tips) ===")
    for name,tips,exp in (("NOVEL clade",nov,105),("Akkermansia clade",akk,187),
                          ("NOVEL+Akkermansia parent",nov+akk,292)):
        node=t.get_common_ancestor(tips)
        got=[l.name for l in node.get_leaves()]
        intr=[x for x in got if x not in tips]
        print(f"  {name:26s} subtends {len(got):3d} (expect {exp})  support {node.support}  intruders {len(intr)}")
        if intr: print("      e.g.", intr[:3])
    sis=t.get_common_ancestor(nov)
    up=sis.up
    sibs=[c for c in up.children if c is not sis] if up else []
    if sibs:
        st=[l.name for l in sibs[0].get_leaves()]
        pure=all(grp(x)=="Akkermansia" for x in st)
        print(f"  sister of NOVEL = {len(st)} tips, all Akkermansia: {pure}")
