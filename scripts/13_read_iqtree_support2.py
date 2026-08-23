#!/usr/bin/env python3
# Support values decoded at the three key nodes
# Source: ch3-chitin-evolution/scripts/read_iqtree_support2.py
# Output: stdout
import re
from ete3 import Tree
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_tree"
IQ=f"{B}/iqtree_placement/akk_placement_iqtree.treefile"

raw=open(IQ).read().strip()
print("labels found:", len(re.findall(r"\)[0-9.]+/[0-9.]+:", raw)))

# turn "SH/UF" into a single numeric label we can recover: encode as SH*1000+UF
enc=re.sub(r"\)([0-9.]+)/([0-9.]+):",
           lambda m: ")%f:" % (float(m.group(1))*1000.0+float(m.group(2))), raw)
t=Tree(enc, format=0)
def grp(n): return n.split("|")[0]
leaves=t.get_leaves()
outg=[l.name for l in leaves if grp(l.name).startswith("OUTGROUP")]
t.set_outgroup(t.get_common_ancestor(outg))
nov=[l.name for l in leaves if grp(l.name)=="NOVEL"]
akk=[l.name for l in leaves if grp(l.name)=="Akkermansia"]

def dec(v):
    sh=int(v)//1000; uf=round(v-sh*1000,1)
    return sh,uf
for name,tips,exp in (("NOVEL clade",nov,105),("Akkermansia clade",akk,187),
                      ("NOVEL+Akk parent",nov+akk,292)):
    n=t.get_common_ancestor(tips)
    got=len(n.get_leaves())
    sh,uf=dec(n.support)
    print(f"{name:20s} {got:3d} tips (expect {exp})   SH-aLRT {sh}   UFBoot {uf}   raw {n.support}")
