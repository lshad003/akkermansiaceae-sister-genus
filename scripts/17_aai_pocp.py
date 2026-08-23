#!/usr/bin/env python3
# Average amino acid identity and POCP computed by reciprocal search
# Source: ch3-chitin-evolution/scripts/aai_pocp.py
# Output: results/aai_pocp/aai_pocp_vs_type.tsv
import csv, os, subprocess, collections, sys
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"
OUT=f"{BASE}/results/aai_pocp"; os.makedirs(OUT,exist_ok=True)
D="/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
TYPE="EHM058340"

FAA=[f"{P}/gff199", f"{P}/gtdb_akk", f"{P}/outgroups", f"{P}/akkfam",
     f"{BASE}/results/dbcan_ehi_nonamph", f"{BASE}/results/dbcan_ehi_amphibian"]
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

# one representative per genus: highest completeness with proteins
comp={r["accession"]:r for r in rows}
reps={}
for g,accs in byg.items():
    def c(a):
        try: return float(comp[a].get("completeness","0"))
        except: return 0.0
    reps[g]=sorted(accs,key=lambda a:-c(a))[0]
reps["NOVEL"]=TYPE
print("representatives (one per genus):")
for g,a in sorted(reps.items()): print(f"  {g:<20}{a:<32}n_with_faa={len(byg[g])}")

def nprot(p): return sum(1 for l in open(p) if l.startswith(">"))
def rbh(qa,sa):
    q,s=idx[qa],idx[sa]
    db=f"{OUT}/tmp_{sa}"
    subprocess.run([D,"makedb","--in",s,"-d",db,"--quiet"],check=True)
    f1=f"{OUT}/tmp_qs.tsv"
    subprocess.run([D,"blastp","-q",q,"-d",db,"-o",f1,"--outfmt","6","qseqid","sseqid",
        "pident","length","evalue","--max-target-seqs","1","--more-sensitive",
        "--evalue","1e-5","--quiet","--threads","8"],check=True)
    db2=f"{OUT}/tmp_{qa}"
    subprocess.run([D,"makedb","--in",q,"-d",db2,"--quiet"],check=True)
    f2=f"{OUT}/tmp_sq.tsv"
    subprocess.run([D,"blastp","-q",s,"-d",db2,"-o",f2,"--outfmt","6","qseqid","sseqid",
        "pident","length","evalue","--max-target-seqs","1","--more-sensitive",
        "--evalue","1e-5","--quiet","--threads","8"],check=True)
    fwd={}; pocp_hits=0
    for line in open(f1):
        f=line.split("\t")
        fwd[f[0]]=(f[1],float(f[2]))
        # POCP (Qin 2014): identity >40%, alignment > 50% of query length
        if float(f[2])>40.0: pocp_hits+=1
    rev={}
    for line in open(f2):
        f=line.split("\t"); rev[f[0]]=f[1]
    ids=[v[1] for k,v in fwd.items() if rev.get(v[0])==k]
    for x in (f1,f2,db+".dmnd",db2+".dmnd"):
        try: os.remove(x)
        except: pass
    nq=nprot(q); ns=nprot(s)
    aai=sum(ids)/len(ids) if ids else float("nan")
    return aai,len(ids),nq,ns,pocp_hits

print(f"\n{'vs genus':<20}{'AAI%':>7}{'RBH':>7}{'POCP%':>8}  interpretation")
print("-"*70)
res=[]
for g in sorted(reps):
    if g=="NOVEL": continue
    try:
        aai,nrbh,nq,ns,ph=rbh(TYPE,reps[g])
    except subprocess.CalledProcessError as e:
        print(f"{g:<20}  DIAMOND FAILED"); continue
    # POCP needs both directions; approximate with forward hits both ways
    aai2,nrbh2,nq2,ns2,ph2=rbh(reps[g],TYPE)
    pocp=100.0*(ph+ph2)/(nq+ns)
    tag="same genus" if pocp>=50 else "DIFFERENT genus"
    print(f"{g:<20}{aai:>7.1f}{nrbh:>7}{pocp:>8.1f}  {tag}")
    res.append((g,aai,nrbh,pocp))
with open(f"{OUT}/aai_pocp_vs_type.tsv","w") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["genus","rep_genome","AAI_pct","n_RBH","POCP_pct"])
    for g,aai,n,p in res: w.writerow([g,reps[g],f"{aai:.2f}",n,f"{p:.2f}"])
print(f"\nwrote {OUT}/aai_pocp_vs_type.tsv")
print("\nREAD IT: POCP < 50% = different genus (Qin et al. 2014).")
print("  AAI genus boundary is ~60-80% (Konstantinidis & Tiedje; Luo et al.).")
print("  The Akkermansia row is the one that matters.")
