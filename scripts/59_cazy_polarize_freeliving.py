#!/usr/bin/env python3
# Enzyme families polarized against free-living outgroups
# Source: ch3-chitin-evolution/scripts/cazy_polarize_freeliving.py
# Output: results/novel_akk_tree/cazy_polarized.tsv
import csv, os, collections, statistics
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUT=f"{BASE}/results/novel_akk_tree/cazy_polarized.tsv"
DIRS=[f"{BASE}/results/dbcan_verru", f"{BASE}/results/dbcan_bact_refs",
      f"{BASE}/results/dbcan_ehi_amphibian", f"{BASE}/results/dbcan_ehi_nonamph",
      f"{BASE}/results/dbcan_endo_allphyla", f"{BASE}/results/dbcan_flavo_refs",
      f"{BASE}/results/dbcan_scaffold", RUM]
E=1e-10; COV=0.30
FREE_GENERA={"Luteolibacter","Haloferula","Rubritalea","Roseibacillus","Roseibacillus_B",
             "SW10","Oceaniferula","UBA956","Persicirhabdus","WTJZ01"}

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t") if r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()

G=collections.OrderedDict()
G["NOVEL"]=[r for r in rows if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian"
            and gen(r) in ("","unknown","NO_GENUS")]
G["AKK_amph"]=[r for r in rows if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian"
               and gen(r)=="Akkermansia"]
G["AKK_all"]=[r for r in rows if r["family"]=="Akkermansiaceae" and gen(r)=="Akkermansia"]
G["FREE"]=[r for r in rows if gen(r) in FREE_GENERA]
for k,v in G.items(): print(f"{k:10s} n={len(v)}")
fg=collections.Counter(gen(r) for r in G["FREE"])
print("free-living genera:", dict(fg.most_common()))

idx={}
for d in DIRS:
    if not os.path.isdir(d): continue
    for fn in os.listdir(d):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            idx.setdefault(fn[:-4], os.path.join(d,fn))
def res(a):
    if a in idx: return idx[a]
    if a[:3] in ("GB_","RS_") and a[3:] in idx: return idx[a[3:]]
    return None
def fam_of(h): return h[:-4].split("_")[0] if h.endswith(".hmm") else h.split("_")[0]

fams_of_genome={}
for k,v in G.items():
    for r in v:
        a=r["accession"]
        if a in fams_of_genome: continue
        p=res(a)
        if not p: continue
        s=set()
        for line in open(p):
            f=line.rstrip("\n").split("\t")
            if len(f)<10: continue
            try: ev=float(f[4]); cv=float(f[9])
            except: continue
            if cv<COV or ev>E: continue
            s.add(fam_of(f[0]))
        fams_of_genome[a]=s

def tot(k): return [len(fams_of_genome[r["accession"]]) for r in G[k] if r["accession"] in fams_of_genome]
print("\nTOTAL CAZy FAMILIES PER GENOME (repertoire size)")
for k in G:
    t=tot(k)
    if t: print("  %-10s n=%3d  median %5.1f  IQR %.0f-%.0f" % (
        k,len(t),statistics.median(t),
        statistics.quantiles(t,n=4)[0] if len(t)>3 else min(t),
        statistics.quantiles(t,n=4)[2] if len(t)>3 else max(t)))

def prev(k,f):
    v=[r["accession"] for r in G[k] if r["accession"] in fams_of_genome]
    if not v: return float("nan"),0,0
    n=sum(1 for a in v if f in fams_of_genome[a])
    return 100.0*n/len(v), n, len(v)

allf=sorted(set().union(*fams_of_genome.values())) if fams_of_genome else []
res_rows=[]
for f in allf:
    pn,_,_=prev("NOVEL",f); pa,_,_=prev("AKK_amph",f)
    pA,_,_=prev("AKK_all",f); pf,_,_=prev("FREE",f)
    res_rows.append((f,pn,pa,pA,pf))

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["family","pct_NOVEL","pct_AKK_amph","pct_AKK_all","pct_FREE","polarity"])
    for f,pn,pa,pA,pf in res_rows:
        pol=""
        if pn>=50 and pA<=10: pol="RETAINED_by_novel_LOST_by_Akk" if pf>=30 else "GAINED_by_novel"
        elif pA>=50 and pn<=10: pol="AKK_specific" if pf<=10 else "LOST_by_novel"
        w.writerow([f,"%.1f"%pn,"%.1f"%pa,"%.1f"%pA,"%.1f"%pf,pol])
print("\nwrote",OUT)

print("\n=== POLARIZED: novel has it, Akkermansia lacks it, FREE-LIVING ALSO HAS IT")
print("    (= present in the ancestor, LOST on the Akkermansia branch) ===")
print("%-8s %8s %8s %8s %8s" % ("family","NOVEL","AKKamph","AKKall","FREE"))
n_ret=0
for f,pn,pa,pA,pf in sorted(res_rows,key=lambda z:-(z[1]-z[3])):
    if pn>=50 and pA<=10 and pf>=30:
        print("%-8s %7.1f%% %7.1f%% %7.1f%% %7.1f%%" % (f,pn,pa,pA,pf)); n_ret+=1
print("  total:",n_ret)

print("\n=== NOVEL-ONLY (free-living also lack it = gained by the novel lineage) ===")
for f,pn,pa,pA,pf in sorted(res_rows,key=lambda z:-z[1]):
    if pn>=50 and pA<=10 and pf<30:
        print("%-8s %7.1f%% %7.1f%% %7.1f%% %7.1f%%" % (f,pn,pa,pA,pf))

print("\n=== AKKERMANSIA-SPECIFIC (Akk has it, novel and free-living lack it) ===")
for f,pn,pa,pA,pf in sorted(res_rows,key=lambda z:-z[3]):
    if pA>=50 and pn<=10 and pf<=20:
        print("%-8s %7.1f%% %7.1f%% %7.1f%% %7.1f%%" % (f,pn,pa,pA,pf))

# ---- repertoire-matched control ----
print("\n=== REPERTOIRE-MATCHED CHECK (does the pattern survive size normalisation?) ===")
tn=[(r["accession"],len(fams_of_genome[r["accession"]])) for r in G["NOVEL"] if r["accession"] in fams_of_genome]
ta=[(r["accession"],len(fams_of_genome[r["accession"]])) for r in G["AKK_amph"] if r["accession"] in fams_of_genome]
lo=max(min(x[1] for x in tn), min(x[1] for x in ta))
hi=min(max(x[1] for x in tn), max(x[1] for x in ta))
mn=[a for a,c in tn if lo<=c<=hi]; ma=[a for a,c in ta if lo<=c<=hi]
print("overlapping repertoire range: %d-%d families | NOVEL %d genomes, AKK_amph %d genomes"%(lo,hi,len(mn),len(ma)))
if mn and ma:
    print("%-8s %10s %10s %8s" % ("family","NOVELmatch","AKKmatch","diff"))
    for f,pn,pa,pA,pf in sorted(res_rows,key=lambda z:-(z[1]-z[2]))[:15]:
        a1=100.0*sum(1 for a in mn if f in fams_of_genome[a])/len(mn)
        a2=100.0*sum(1 for a in ma if f in fams_of_genome[a])/len(ma)
        print("%-8s %9.1f%% %9.1f%% %+7.1f" % (f,a1,a2,a1-a2))
