#!/usr/bin/env python3
# Representative proteome set assembled
# Source: ch3-chitin-evolution/scripts/build_orthofinder_set.py
# Output: results/pangenome/orthofinder_in
import os, csv, shutil, collections
from ete3 import Tree

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"
TREE=f"{BASE}/results/novel_akk_tree/akk_placement.nwk"
TRI=f"{BASE}/results/novel_akk_ani/novel_triangle.tsv"
OF=f"{P}/orthofinder_in"
os.makedirs(OF, exist_ok=True)

SRC=[f"{P}/gff199", f"{P}/outgroups", f"{P}/akkfam", f"{P}/gtdb_akk"]
faa={}
for d in SRC:
    if not os.path.isdir(d): continue
    for fn in os.listdir(d):
        if fn.endswith(".faa"): faa.setdefault(fn[:-4], os.path.join(d,fn))
print("proteomes indexed:", len(faa))

def find(a):
    if a in faa: return faa[a]
    if a[:3] in ("GB_","RS_") and a[3:] in faa: return faa[a[3:]]
    for p in ("GB_","RS_"):
        if p+a in faa: return faa[p+a]
    return None

t=None
for fmt in (1,0,5,3):
    try: t=Tree(TREE, format=fmt); break
    except Exception: continue
def fld(n,i):
    p=n.split("|"); return p[i] if i<len(p) else ""
def acc(n): return n.split("|")[-1]
leaves=t.get_leaves()
outg=[l.name for l in leaves if fld(l.name,0).startswith("OUTGROUP")]
t.set_outgroup(t.get_common_ancestor(outg))
leaves=t.get_leaves()

# ---- NOVEL: one rep per 95% ANI cluster ----
def base(x):
    x=os.path.basename(x.strip())
    for e in (".fa",".fna",".fasta"):
        if x.endswith(e): x=x[:-len(e)]
    return x
raw=[l.rstrip("\n") for l in open(TRI) if l.strip()]
st=1 if raw[0].strip().isdigit() else 0
names=[];vals=[]
for l in raw[st:]:
    q=l.split("\t"); names.append(base(q[0])); vals.append(q[1:])
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

sel=[]; miss=[]
for i,(k,mem) in enumerate(sorted(cl.items(), key=lambda kv:-len(kv[1])),1):
    rep=next((m for m in mem if find(m)), None)
    if rep: sel.append(("NOVEL", f"NOVEL_sp{i:02d}_n{len(mem)}", rep))
    else: miss.append(("NOVEL cluster %d"%i, mem[0]))

# ---- Akkermansia: one rep per monophyletic clade ----
akkset=set(l.name for l in leaves if fld(l.name,0)=="Akkermansia")
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
for thr in [x/2000.0 for x in range(1,600)]:
    g=collapse(thr)
    if sum(len(x) for x in g)==len(akkset) and 14<=len(g)<=24:
        best=(thr,g); break
if best is None: best=(0.05, collapse(0.05))
thr,groups=best
print("Akkermansia clades:", len(groups), "threshold", thr)
for i,g in enumerate(sorted(groups,key=len,reverse=True),1):
    hosts=collections.Counter(fld(x,2) for x in g)
    rep=next((acc(x) for x in g if find(acc(x))), None)
    tag="-".join(f"{h}{n}" for h,n in hosts.most_common()[:2])
    if rep: sel.append(("AKK", f"AKK_clade{i:02d}_n{len(g)}_{tag}", rep))
    else: miss.append(("AKK clade %d"%i, acc(g[0])))

# ---- free-living: reps per genus (2 for genera >=30 genomes, else 1) ----
byg=collections.defaultdict(list)
for l in leaves:
    g=fld(l.name,0)
    if g in ("NOVEL","Akkermansia"): continue
    byg[g].append(acc(l.name))
for g,mem in sorted(byg.items(), key=lambda kv:-len(kv[1])):
    k = 2 if len(mem)>=30 else 1
    got=0
    for a in mem:
        if got>=k: break
        if find(a):
            sel.append(("FREE", f"{g}_rep{got+1}_of{len(mem)}", a)); got+=1
    if got==0: miss.append((g, mem[0]))

print("\n%-6s %-38s %s" % ("set","label","accession"))
for s,lab,a in sel: print("%-6s %-38s %s" % (s,lab,a))
c=collections.Counter(s for s,_,_ in sel)
print("\nSELECTED:", dict(c), " total", len(sel))
if miss:
    print("\nNO PROTEOME FOUND for:")
    for m in miss: print("  ",m)

n=0
for s,lab,a in sel:
    p=find(a)
    if not p: continue
    shutil.copy(p, os.path.join(OF, lab+".faa")); n+=1
print("\ncopied %d proteomes to %s" % (n, OF))
