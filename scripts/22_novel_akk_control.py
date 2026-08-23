#!/usr/bin/env python3
# Screening floor measured on named reference genomes
# Source: ch3-chitin-evolution/scripts/novel_akk_control.py
# Output: results/novel_akk_ani/named_vs_akkermansia.tsv
import csv, os, subprocess, collections, sys
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SKANI="/bigdata/stajichlab/lshad003/condaenvs/drep/bin/skani"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
FA=f"{BASE}/data/amphibia_gtdbtk_input"
OUT=f"{BASE}/results/novel_akk_ani"
cen={r["accession"]:r for r in csv.DictReader(open(CEN),delimiter="\t")}
def g(r): return (r["genus"] or "").strip()

named=[a for a,r in cen.items() if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian"
       and r["annotated"]=="1" and g(r)=="Akkermansia" and os.path.exists(f"{FA}/{a}.fa")]
print(f"POSITIVE CONTROL: {len(named)} NAMED amphibian Akkermansia vs the same 60 refs")
cl=f"{OUT}/named_ctrl.txt"
with open(cl,"w") as fh:
    for a in named: fh.write(f"{FA}/{a}.fa\n")
mc=f"{OUT}/named_vs_akkermansia.tsv"
r=subprocess.run([SKANI,"dist","--ql",cl,"--rl",f"{OUT}/refs.txt","-o",mc,"-t","8"],
                 capture_output=True,text=True)
if r.returncode!=0: print("FAILED\n",r.stderr[-1500:]); sys.exit(1)
h=collections.defaultdict(list)
for x in csv.DictReader(open(mc),delimiter="\t"):
    try: h[os.path.basename(x["Query_file"]).replace(".fa","")].append(float(x["ANI"]))
    except: pass
print(f"  named genomes WITH a hit: {len(h)}/{len(named)}")
if h:
    b=sorted(max(v) for v in h.values())
    print(f"  best-hit ANI: min {b[0]:.1f}  median {b[len(b)//2]:.1f}  max {b[-1]:.1f}")
    print("  VERDICT: skani works. The 0/105 for unnamed genomes is REAL.")
else:
    print("  VERDICT: skani returns 0 for KNOWN Akkermansia too. SCRIPT IS BROKEN, ignore the 0/105.")

print("\n=== CLUSTER COMPOSITION at 80% ANI: do EHI and UHM converge? ===")
lines=[l.rstrip("\n") for l in open(f"{OUT}/novel_triangle.tsv") if l.strip()]
names=[];vals=[]
for l in lines[1:]:
    f=l.split("\t"); names.append(os.path.basename(f[0]).replace(".fa",""))
    vals.append([float(x) if x not in ("","NA") else 0.0 for x in f[1:]])
par={x:x for x in names}
def find(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
for i in range(len(names)):
    for j in range(i+1,len(names)):
        v=vals[i][j] if j<len(vals[i]) else 0.0
        if v>=80.0:
            a,b=find(names[i]),find(names[j])
            if a!=b: par[b]=a
cl2=collections.defaultdict(list)
for x in names: cl2[find(x)].append(x)
def animal(a):
    r=cen[a]
    return a.split(".")[0] if r["from_dataset"]=="herptile_MAG" else (r.get("sample_id") or "?")
for k,v in sorted(cl2.items(), key=lambda x:-len(x[1])):
    ds=collections.Counter(cen[a]["from_dataset"] for a in v)
    hs=collections.Counter(cen[a].get("host_animal_type","?") for a in v)
    an=len(set(animal(a) for a in v))
    print(f"\n  cluster n={len(v)}  animals={an}")
    print(f"    datasets: {dict(ds)}")
    print(f"    hosts:    {dict(hs.most_common(4))}")
