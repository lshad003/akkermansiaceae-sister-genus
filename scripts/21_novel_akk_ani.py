#!/usr/bin/env python3
# Nucleotide identity triangle and species clustering
# Source: ch3-chitin-evolution/scripts/novel_akk_ani.py
# Output: results/novel_akk_ani/novel_triangle.tsv
import csv, os, subprocess, sys, gzip, shutil, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SKANI="/bigdata/stajichlab/lshad003/condaenvs/drep/bin/skani"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
PATHS=f"{BASE}/data/verru_annotation_paths.txt"
FA=f"{BASE}/data/amphibia_gtdbtk_input"
OUT=f"{BASE}/results/novel_akk_ani"; os.makedirs(OUT,exist_ok=True)
TMP=f"{OUT}/ref_fna"; os.makedirs(TMP,exist_ok=True)
THREADS=8

cen={r["accession"]:r for r in csv.DictReader(open(CEN),delimiter="\t")}
def g(r): return (r["genus"] or "").strip()
novel=[a for a,r in cen.items() if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian"
       and r["annotated"]=="1" and g(r) in ("","unknown","NO_GENUS")]
refs=[a for a,r in cen.items() if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"
      and g(r)=="Akkermansia" and r["from_dataset"]=="GTDB_r226_rep"]
print(f"novel (unnamed amphibian Akkermansiaceae): {len(novel)}")
print(f"named Akkermansia GTDB refs:               {len(refs)}")

pmap={}
for l in open(PATHS):
    l=l.strip()
    if not l: continue
    b=os.path.basename(l).replace("_genomic.fna.gz","")
    pmap[b]=l
def refpath(acc):
    for k in (acc, acc[3:] if acc[:3] in ("GB_","RS_") else acc):
        if k in pmap: return pmap[k]
    return None

nlist=f"{OUT}/novel.txt"
with open(nlist,"w") as fh:
    for a in novel: fh.write(f"{FA}/{a}.fa\n")

rlist=f"{OUT}/refs.txt"; got=0; missing=[]
with open(rlist,"w") as fh:
    for a in refs:
        p=refpath(a)
        if not p or not os.path.exists(p): missing.append(a); continue
        d=f"{TMP}/{a}.fna"
        if not os.path.exists(d):
            try:
                with gzip.open(p,"rb") as i, open(d,"wb") as o: shutil.copyfileobj(i,o)
            except Exception as e: missing.append(a); continue
        fh.write(d+"\n"); got+=1
print(f"ref FASTAs resolved: {got}/{len(refs)}   missing {len(missing)}")
if missing[:5]: print("  missing e.g.",missing[:5])
if got==0: print("NO REFS RESOLVED, stopping"); sys.exit(1)

m1=f"{OUT}/novel_triangle.tsv"
print(f"\n[1/2] skani triangle on {len(novel)} novel genomes")
r=subprocess.run([SKANI,"triangle","-l",nlist,"-o",m1,"--full-matrix","-t",str(THREADS)],
                 capture_output=True,text=True)
if r.returncode!=0: print("FAILED\n",r.stderr[-2000:]); sys.exit(1)

m2=f"{OUT}/novel_vs_akkermansia.tsv"
print(f"[2/2] skani dist: {len(novel)} novel vs {got} named Akkermansia")
r=subprocess.run([SKANI,"dist","--ql",nlist,"--rl",rlist,"-o",m2,"-t",str(THREADS)],
                 capture_output=True,text=True)
if r.returncode!=0: print("FAILED\n",r.stderr[-2000:]); sys.exit(1)

lines=[l.rstrip("\n") for l in open(m1) if l.strip()]
names=[];vals=[]
for l in lines[1:]:
    f=l.split("\t"); names.append(os.path.basename(f[0]).replace(".fa",""))
    vals.append([float(x) if x not in ("","NA") else 0.0 for x in f[1:]])
def cluster(th):
    par={g:g for g in names}
    def find(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            v=vals[i][j] if j<len(vals[i]) else 0.0
            if v>=th:
                a,b=find(names[i]),find(names[j])
                if a!=b: par[b]=a
    c=collections.defaultdict(list)
    for gname in names: c[find(gname)].append(gname)
    return c
print("\n=== ARE THE 105 ONE LINEAGE? ===")
for th in (95,90,85,80,75):
    c=cluster(th)
    big=max((len(v) for v in c.values()), default=0)
    print(f"  >= {th}% ANI: {len(c)} clusters, largest {big}/{len(names)}")
nz=[vals[i][j] for i in range(len(names)) for j in range(i+1,len(names))
    if j<len(vals[i]) and vals[i][j]>0]
if nz:
    nz.sort()
    print(f"  pairwise ANI (nonzero, n={len(nz)}): min {nz[0]:.1f}  median {nz[len(nz)//2]:.1f}  max {nz[-1]:.1f}")
print(f"  pairs with NO ANI computable: {len(names)*(len(names)-1)//2 - len(nz)}")

print("\n=== ARE THEY OUTSIDE AKKERMANSIA? ===")
hits=collections.defaultdict(list)
for r_ in csv.DictReader(open(m2),delimiter="\t"):
    q=os.path.basename(r_["Query_file"]).replace(".fa","")
    try: hits[q].append(float(r_["ANI"]))
    except: pass
print(f"  novel genomes with ANY ANI hit to a named Akkermansia: {len(hits)}/{len(novel)}")
print(f"  novel genomes with NO hit at all:                      {len(novel)-len(hits)}/{len(novel)}")
if hits:
    best=sorted((max(v),k) for k,v in hits.items())
    b=[x[0] for x in best]
    print(f"  best-hit ANI: min {b[0]:.1f}  median {b[len(b)//2]:.1f}  max {b[-1]:.1f}")
    print("  5 highest:", [f"{v:.1f}" for v,_ in best[-5:]])
print(f"\nwrote {m1}\nwrote {m2}")
