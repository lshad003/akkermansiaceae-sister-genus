#!/usr/bin/env python3
# Figure 1 tree and clade colour annotations
# Source: ch3-chitin-evolution/scripts/fig1_akk_clades_iqtree.py
# Output: results/figures/fig1_collapsed/fig1_iqtree.nwk, results/figures/fig1_collapsed/fig1_iqtree_colors.txt
import os, collections, re
from ete3 import Tree

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
TREE=f"{BASE}/results/novel_akk_tree/iqtree_placement/akk_placement_iqtree.treefile"
TRI=f"{BASE}/results/novel_akk_ani/novel_triangle.tsv"
OUTD=f"{BASE}/results/figures/fig1_collapsed"
os.makedirs(OUTD, exist_ok=True)

raw=open(TREE).read().strip()
print("support labels:", len(re.findall(r"\)[0-9.]+/[0-9.]+:", raw)))
t=Tree(raw, format=1)   # keeps "SH/UF" as internal node names

def fld(n,i):
    p=n.split("|"); return p[i] if i<len(p) else ""
def grp(n): return fld(n,0)
def acc(n): return n.split("|")[-1]
DSN={"herptile_MAG":"UHM","EHI_2025":"EHI","GTDB_r226_rep":"GTDB"}

leaves=t.get_leaves()
outg=[l.name for l in leaves if grp(l.name).startswith("OUTGROUP")]
t.set_outgroup(t.get_common_ancestor(outg))
leaves=t.get_leaves()
print("tips:",len(leaves))

def base(x):
    x=os.path.basename(x.strip())
    for e in (".fa",".fna",".fasta"):
        if x.endswith(e): x=x[:-len(e)]
    return x
rows=[l.rstrip("\n") for l in open(TRI) if l.strip()]
st=1 if rows[0].strip().isdigit() else 0
names=[];vals=[]
for l in rows[st:]:
    p=l.split("\t"); names.append(base(p[0])); vals.append(p[1:])
par={}
def f_(a):
    par.setdefault(a,a)
    while par[a]!=a: par[a]=par[par[a]]; a=par[a]
    return a
def u_(a,b):
    ra,rb=f_(a),f_(b)
    if ra!=rb: par[ra]=rb
for n in names: f_(n)
for i,row in enumerate(vals):
    for j,x in enumerate(row):
        try: v=float(x)
        except: continue
        if j<len(names) and i!=j and v>=95.0: u_(names[i],names[j])
cl=collections.defaultdict(list)
for n in names: cl[f_(n)].append(n)

acc2tip={acc(l.name):l.name for l in leaves}
keep=[]; label={}
for i,(k,mem) in enumerate(sorted(cl.items(), key=lambda kv:-len(kv[1])),1):
    tip=next((acc2tip[m] for m in mem if m in acc2tip), None)
    if tip is None: continue
    ehi=sum(1 for m in mem if m.startswith("EHM"))
    ds="EHI" if ehi==len(mem) else ("UHM" if ehi==0 else "mixed")
    keep.append(tip)
    g="genome" if len(mem)==1 else "genomes"
    label[tip]=f"NOVEL_sp{i:02d}_{ds}_{len(mem)}{g}"

akkset=set(l.name for l in leaves if grp(l.name)=="Akkermansia")
akknode=t.get_common_ancestor(list(akkset))
def collapse(thr):
    groups=[]; done=set()
    for node in akknode.traverse("preorder"):
        if id(node) in done: continue
        lv=[l.name for l in node.get_leaves()]
        if not lv or not set(lv)<=akkset: continue
        d=max([node.get_distance(l) for l in node.get_leaves()] or [0])
        if d<=thr:
            groups.append(lv)
            for x in node.traverse(): done.add(id(x))
    return groups
best=None
for thr in [x/2000.0 for x in range(1,800)]:
    g=collapse(thr)
    if sum(len(x) for x in g)==len(akkset) and 14<=len(g)<=24:
        best=(thr,g); break
if best is None: best=(0.05, collapse(0.05))
thr,groups=best
print("Akkermansia clades:",len(groups),"threshold",thr)
for i,g in enumerate(sorted(groups,key=len,reverse=True),1):
    hosts=collections.Counter(fld(x,2) for x in g)
    dss=collections.Counter(DSN.get(fld(x,1),fld(x,1)) for x in g)
    htxt="-".join(f"{h}{n}" for h,n in hosts.most_common()[:2])
    dtxt="-".join(f"{d}{n}" for d,n in dss.most_common()[:2])
    keep.append(g[0]); label[g[0]]=f"Akkermansia_clade{i:02d}_n{len(g)}_{htxt}_{dtxt}"

byg=collections.defaultdict(list)
for l in leaves:
    if grp(l.name) not in ("NOVEL","Akkermansia"): byg[grp(l.name)].append(l.name)
for g,mem in sorted(byg.items(), key=lambda kv:-len(kv[1])):
    if len(mem)<3: continue
    keep.append(mem[0]); label[mem[0]]=f"{g}_{len(mem)}genomes"

t.prune(keep, preserve_branch_length=True)
for l in t.get_leaves():
    l.name=re.sub(r"[^A-Za-z0-9_.-]","",label.get(l.name,l.name))

nwk=f"{OUTD}/fig1_iqtree.nwk"
t.write(outfile=nwk, format=1)
out=["TREE_COLORS","SEPARATOR TAB","DATA"]
for l in t.get_leaves():
    n=l.name
    if n.startswith("NOVEL"):
        c="#c96a3f" if "_EHI_" in n else ("#7a5ea8" if "_mixed_" in n else "#4a7ebb")
    elif n.startswith("Akkermansia"):
        c="#d9a441" if "amphibian" in n else "#b8b6ac"
    else: continue
    out.append(f"{n}\trange\t{c}\t{n}")
open(f"{OUTD}/fig1_iqtree_colors.txt","w").write("\n".join(out)+"\n")
print("wrote",nwk)
print("wrote",f"{OUTD}/fig1_iqtree_colors.txt")
print("tips kept:",len(t.get_leaves()))
