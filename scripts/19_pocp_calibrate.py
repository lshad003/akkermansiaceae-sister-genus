#!/usr/bin/env python3
# POCP calibrated on known genus pairs
# Source: ch3-chitin-evolution/scripts/pocp_calibrate.py
# Output: results/aai_pocp/pocp_calibration.tsv
import csv, os, subprocess, collections, itertools
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"; OUT=f"{BASE}/results/aai_pocp"
D="/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
FAA=[f"{P}/gff199",f"{P}/gtdb_akk",f"{P}/outgroups",f"{P}/akkfam",
     f"{BASE}/results/dbcan_ehi_nonamph",f"{BASE}/results/dbcan_ehi_amphibian"]
idx={}
for d in FAA:
    if not os.path.isdir(d): continue
    for f in os.listdir(d):
        if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d,f))
rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"]
def gn(r):
    g=(r["genus"] or "").strip()
    return "NOVEL" if g in ("","unknown","NO_GENUS") else g
byg=collections.defaultdict(list)
for r in rows:
    if r["accession"] in idx: byg[gn(r)].append(r["accession"])
comp={r["accession"]:r for r in rows}
def c(x):
    try: return float(comp[x].get("completeness","0"))
    except: return 0.0
reps={g:sorted(a,key=lambda x:-c(x))[0] for g,a in byg.items()}
reps["NOVEL"]="EHM058340"

def nprot(p): return sum(1 for l in open(p) if l.startswith(">"))
def cons(qa,sa):
    q,s=idx[qa],idx[sa]; db=f"{OUT}/c_{sa}"
    subprocess.run([D,"makedb","--in",s,"-d",db,"--quiet"],check=True)
    o=f"{OUT}/c_{qa}_{sa}.tsv"
    subprocess.run([D,"blastp","-q",q,"-d",db,"-o",o,"--outfmt","6","qseqid","sseqid",
        "pident","length","qlen","evalue","--max-target-seqs","1","--more-sensitive",
        "--evalue","1e-5","--quiet","--threads","8"],check=True)
    n=0
    for line in open(o):
        f=line.split("\t")
        if float(f[2])>40.0 and int(f[4])>0 and int(f[3])/int(f[4])>0.5: n+=1
    for x in (o,db+".dmnd"):
        try: os.remove(x)
        except: pass
    return n,nprot(q)
def pocp(a,b):
    c1,t1=cons(a,b); c2,t2=cons(b,a)
    return 100.0*(c1+c2)/(t1+t2)

BIG=[g for g in reps if len(byg.get(g,[]))>=5 and g!="NOVEL"]
print("CALIBRATION: POCP between ACCEPTED genera of Akkermansiaceae")
print("If accepted genus pairs also exceed 50%, the threshold does not")
print("discriminate in this family and the NOVEL 51.8% is uninformative.\n")
print(f"{'pair':<42}{'POCP%':>8}  verdict")
print("-"*62)
out=[]
for a,b in itertools.combinations(sorted(BIG),2):
    v=pocp(reps[a],reps[b])
    out.append((a,b,v))
    print(f"{a+' vs '+b:<42}{v:>8.1f}  {'SAME (threshold fails)' if v>=50 else 'different'}")
v=[x[2] for x in out]
v.sort()
print(f"\naccepted-genus pairs: n={len(v)}  min {v[0]:.1f}  median {v[len(v)//2]:.1f}  max {v[-1]:.1f}")
print(f"pairs at or above 50%: {sum(1 for x in v if x>=50)}/{len(v)}")
print(f"\nNOVEL vs Akkermansia = 51.8  <- compare to the max above")
with open(f"{OUT}/pocp_calibration.tsv","w") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["genusA","genusB","POCP_pct"])
    for a,b,x in out: w.writerow([a,b,f"{x:.2f}"])
print(f"wrote {OUT}/pocp_calibration.tsv")
