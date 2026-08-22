#!/usr/bin/env python3
# POCP recomputed against a fixed reference
# Source: ch3-chitin-evolution/scripts/pocp_fixed.py
# Output: results/aai_pocp/pocp_fixed.tsv
import csv, os, subprocess, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"; OUT=f"{BASE}/results/aai_pocp"
D="/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
TYPE="EHM058340"
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
reps={}
for g,a in byg.items():
    def c(x):
        try: return float(comp[x].get("completeness","0"))
        except: return 0.0
    reps[g]=sorted(a,key=lambda x:-c(x))[0]
reps["NOVEL"]=TYPE

def nprot(p): return sum(1 for l in open(p) if l.startswith(">"))
def conserved(qa,sa):
    """Qin 2014: evalue<1e-5, identity>40%, alignment >50% of QUERY length."""
    q,s=idx[qa],idx[sa]
    db=f"{OUT}/d_{sa}"
    subprocess.run([D,"makedb","--in",s,"-d",db,"--quiet"],check=True)
    o=f"{OUT}/h_{qa}_{sa}.tsv"
    subprocess.run([D,"blastp","-q",q,"-d",db,"-o",o,"--outfmt","6","qseqid","sseqid",
        "pident","length","qlen","evalue","--max-target-seqs","1","--more-sensitive",
        "--evalue","1e-5","--quiet","--threads","8"],check=True)
    n=0; raw=0
    for line in open(o):
        f=line.split("\t")
        pid=float(f[2]); alen=int(f[3]); qlen=int(f[4])
        raw+=1
        if pid>40.0 and qlen>0 and (alen/qlen)>0.5: n+=1
    for x in (o,db+".dmnd"):
        try: os.remove(x)
        except: pass
    return n, raw, nprot(q)

print(f"POCP, Qin 2014 criteria (e<1e-5, id>40%, aln>50% of query)")
print(f"type genome: {TYPE}\n")
print(f"{'vs genus':<18}{'C1':>6}{'C2':>6}{'T1':>6}{'T2':>6}{'POCP%':>8}  verdict")
print("-"*62)
res=[]
for g in sorted(reps):
    if g=="NOVEL": continue
    c1,r1,t1=conserved(TYPE,reps[g])
    c2,r2,t2=conserved(reps[g],TYPE)
    pocp=100.0*(c1+c2)/(t1+t2)
    v="SAME genus" if pocp>=50 else "different genus"
    print(f"{g:<18}{c1:>6}{c2:>6}{t1:>6}{t2:>6}{pocp:>8.1f}  {v}")
    res.append((g,reps[g],c1,c2,t1,t2,pocp))
with open(f"{OUT}/pocp_fixed.tsv","w") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["genus","rep","C1","C2","T1","T2","POCP_pct"])
    for r in res: w.writerow(list(r[:6])+[f"{r[6]:.2f}"])
print(f"\nwrote {OUT}/pocp_fixed.tsv")
