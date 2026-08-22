#!/usr/bin/env python3
# POCP stability across alternative representatives
# Source: ch3-chitin-evolution/scripts/pocp_stability.py
# Output: results/aai_pocp/pocp_stability.tsv
import csv, os, subprocess, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"; OUT=f"{BASE}/results/aai_pocp"
D="/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
AKK_REP="RS_GCF_040616545.1"          # same Akkermansia rep as the 51.8 run

faa={}
for d in (f"{P}/gff199",f"{P}/gtdb_akk",f"{P}/outgroups",f"{P}/akkfam",
          f"{BASE}/results/dbcan_ehi_nonamph",f"{BASE}/results/dbcan_ehi_amphibian"):
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): faa.setdefault(f[:-4], os.path.join(d,f))

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1" and r["accession"] in faa]
def gen(r): return (r["genus"] or "").strip()
NOV=[r for r in rows if r["host_class"]=="amphibian" and gen(r) in ("","unknown","NO_GENUS")]
print(f"NOVEL genomes with proteins: {len(NOV)}")
def comp(r):
    try: return float(r["completeness"])
    except: return 0.0

# pick reps spanning the size range: use protein count as the proxy
def np_(a): return sum(1 for l in open(faa[a]) if l.startswith(">"))
cand=[(np_(r["accession"]), r["accession"], comp(r)) for r in NOV if comp(r)>=90]
cand.sort()
picks=[cand[0], cand[len(cand)//4], cand[len(cand)//2], cand[3*len(cand)//4], cand[-1]]
if "EHM058340" not in [p[1] for p in picks]:
    picks.append((np_("EHM058340"),"EHM058340",100.0))
print(f"\nreps chosen (spanning the protein-count range, completeness>=90):")
for n,a,c in picks: print(f"  {a:<30}{n:>6} proteins  comp {c:.1f}%")

def cons(qa,sa):
    q,s=faa[qa],faa[sa]; db=f"{OUT}/s_{sa}"
    subprocess.run([D,"makedb","--in",s,"-d",db,"--quiet"],check=True)
    o=f"{OUT}/s_{qa}_{sa}.tsv"
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
    return n, sum(1 for l in open(q) if l.startswith(">"))

print(f"\nPOCP vs Akkermansia ({AKK_REP}), Qin 2014 criteria")
print(f"{'NOVEL rep':<30}{'C1':>6}{'C2':>6}{'T1':>6}{'T2':>6}{'POCP%':>8}  verdict")
print("-"*76)
res=[]
for _,a,_ in picks:
    c1,t1=cons(a,AKK_REP); c2,t2=cons(AKK_REP,a)
    p=100.0*(c1+c2)/(t1+t2)
    v="SAME genus" if p>=50 else "different genus"
    print(f"{a:<30}{c1:>6}{c2:>6}{t1:>6}{t2:>6}{p:>8.1f}  {v}")
    res.append((a,p))
v=sorted(x[1] for x in res)
print(f"\nPOCP range across {len(v)} representatives: {v[0]:.1f} to {v[-1]:.1f}")
print(f"reps giving >=50 (SAME genus): {sum(1 for x in v if x>=50)}/{len(v)}")
print("\nREAD IT: if some reps fall below 50, POCP is unstable at this boundary and the")
print("  51.8 dissent is an artifact of representative choice. If ALL are >=50, the")
print("  dissent is real and must be reported.")
with open(f"{OUT}/pocp_stability.tsv","w") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["novel_rep","POCP_vs_Akkermansia"])
    for a,p in res: w.writerow([a,f"{p:.2f}"])
print(f"wrote {OUT}/pocp_stability.tsv")
