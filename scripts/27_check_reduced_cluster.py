#!/usr/bin/env python3
# Reduced-genome cluster characterized
# Source: ch3-chitin-evolution/scripts/check_reduced_cluster.py
# Output: stdout
import csv, os, statistics, collections

TRI="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_ani/novel_triangle.tsv"
SGC="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novel_akk_tree/novel_size_gc.tsv"
THR=95.0

raw=[l.rstrip("\n") for l in open(TRI) if l.strip()]
print("first 3 lines of triangle file:")
for l in raw[:3]: print("   ", l[:140])

def base(x):
    x=os.path.basename(x.strip())
    for ext in (".fa",".fna",".fasta"): 
        if x.endswith(ext): x=x[:-len(ext)]
    return x

pairs=[]
hdr=raw[0].split("\t")
if any("ANI" in c for c in hdr) and len(hdr)>=3:
    # pairwise format
    r=csv.DictReader(open(TRI),delimiter="\t")
    kq=[c for c in r.fieldnames if "Query" in c or "query" in c][0]
    kr=[c for c in r.fieldnames if "Ref" in c or "ref" in c][0]
    ka=[c for c in r.fieldnames if c.strip()=="ANI" or "ANI" in c][0]
    for row in r:
        try: v=float(row[ka])
        except: continue
        pairs.append((base(row[kq]),base(row[kr]),v))
    print("\nparsed as PAIRWISE, rows:",len(pairs))
else:
    # skani triangle: line0 = n, then name + values
    names=[]; vals=[]
    start=1 if raw[0].strip().isdigit() else 0
    for l in raw[start:]:
        f=l.split("\t")
        names.append(base(f[0]))
        vals.append([x for x in f[1:]])
    for i,row in enumerate(vals):
        for j,x in enumerate(row):
            try: v=float(x)
            except: continue
            if j<len(names): pairs.append((names[i],names[j],v))
    print("\nparsed as TRIANGLE, genomes:",len(names)," pairs:",len(pairs))

# single-linkage at 95
parent={}
def find(a):
    parent.setdefault(a,a)
    while parent[a]!=a:
        parent[a]=parent[parent[a]]; a=parent[a]
    return a
def union(a,b):
    ra,rb=find(a),find(b)
    if ra!=rb: parent[ra]=rb

sz={}; gc={}
for row in csv.DictReader(open(SGC),delimiter="\t"):
    sz[row["genome"]]=float(row["genome_size_bp"]); gc[row["genome"]]=float(row["gc_percent"])
for g in sz: find(g)
for a,b,v in pairs:
    if a!=b and v>=THR: union(a,b)

cl=collections.defaultdict(list)
for g in sz: cl[find(g)].append(g)
print("\nclusters at %.0f%% ANI: %d  (genomes with size/GC: %d)" % (THR,len(cl),len(sz)))

out=[]
for k,mem in cl.items():
    s=[sz[m] for m in mem]; c=[gc[m] for m in mem]
    out.append((len(mem), statistics.median(s)/1e6, statistics.median(c), mem))
out.sort(key=lambda x:-x[0])
print("\n%-6s %-12s %-8s" % ("n","med_Mb","med_GC"))
for n,ms,mg,mem in out:
    flag=""
    if mg<45.5 or ms<2.6: flag="   <-- REDUCED-like"
    print("%-6d %-12.3f %-8.2f%s" % (n,ms,mg,flag))

red=[o for o in out if o[2]<45.5 or o[1]<2.6]
print("\nreduced-like clusters: %d, total genomes: %d" % (len(red), sum(o[0] for o in red)))
if len(red)==1 and red[0][0]==18:
    print("VERDICT: the reduced set IS exactly one 95%% ANI cluster (n=18). Sentence OK.")
else:
    print("VERDICT: reduced set is NOT a single n=18 cluster. Rewrite the sentence.")
    for n,ms,mg,mem in red:
        print("   cluster n=%d med %.2f Mb / %.2f GC" % (n,ms,mg))
