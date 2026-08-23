#!/usr/bin/env python3
# Chitin panel prevalence across lineage-matched groups
# Source: ch3-chitin-evolution/scripts/chitin_panel_5groups.py
# Output: results/novel_akk_tree/chitin_panel_5groups.tsv
import csv, os, math, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUT=f"{BASE}/results/novel_akk_tree/chitin_panel_5groups.tsv"
DIRS=[f"{BASE}/results/dbcan_verru",f"{BASE}/results/dbcan_bact_refs",
      f"{BASE}/results/dbcan_ehi_amphibian",f"{BASE}/results/dbcan_ehi_nonamph",
      f"{BASE}/results/dbcan_endo_allphyla",f"{BASE}/results/dbcan_flavo_refs",
      f"{BASE}/results/dbcan_scaffold",RUM]
PANEL=["GH18","GH19","AA10","CBM5","CBM12","CBM14","CE4","CE14","GH8","GH46","GH75","GH3","GH20","GH84","CE11"]
MODULE={"GH18":"A attack","GH19":"A attack","AA10":"A attack","CBM5":"A bind","CBM12":"A bind",
        "CBM14":"A bind","CE4":"B deacetyl","CE14":"B deacetyl","GH8":"C chitosan",
        "GH46":"C chitosan","GH75":"C chitosan","GH3":"D terminal","GH20":"D terminal",
        "GH84":"D terminal","CE11":"control LpxC"}
E1,E2,COV=1e-10,1e-9,0.30

def logf(n): return math.lgamma(n+1)
def hlp(a,b,c,d): return (logf(a+b)+logf(c+d)+logf(a+c)+logf(b+d)-logf(a+b+c+d)-logf(a)-logf(b)-logf(c)-logf(d))
def fisher(a,b,c,d):
    if (a+b)==0 or (c+d)==0 or (a+c)==0 or (b+d)==0: return 1.0
    obs=hlp(a,b,c,d); r1=a+b; c1=a+c; lo=max(0,c1-(c+d)); hi=min(r1,c1); t=0.0
    for x in range(lo,hi+1):
        lp=hlp(x,r1-x,c1-x,r2 if False else (c+d)-(c1-x))
        if lp<=obs+1e-9: t+=math.exp(lp)
    return min(1.0,t)

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()
G=collections.OrderedDict()
G["NOVEL_amph"]=[r for r in rows if r["host_class"]=="amphibian" and gen(r) in ("","unknown","NO_GENUS")]
G["AKK_amph"]=[r for r in rows if r["host_class"]=="amphibian" and gen(r)=="Akkermansia"]
G["AKK_podarcis"]=[r for r in rows if r["host_class"]=="reptile" and gen(r)=="Akkermansia"]
G["AKK_mammal"]=[r for r in rows if r["host_class"]=="mammal" and gen(r)=="Akkermansia"]
G["AKK_gtdb"]=[r for r in rows if r["host_class"]=="unknown" and gen(r)=="Akkermansia"]
print("GROUPS")
for k,v in G.items(): print(f"  {k:<16}{len(v)}")

idx={}
for d in DIRS:
    if not os.path.isdir(d): continue
    for fn in os.listdir(d):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"): idx.setdefault(fn[:-4],os.path.join(d,fn))
def res(a):
    if a in idx: return idx[a]
    if a[:3] in ("GB_","RS_") and a[3:] in idx: return idx[a[3:]]
    return None
allr=[r for v in G.values() for r in v]
un=[r["accession"] for r in allr if res(r["accession"]) is None]
print(f"resolved {len(allr)-len(un)}/{len(allr)}")
if un: print("UNRESOLVED",un[:5]); raise SystemExit(1)

def famof(h): return h[:-4].split("_")[0] if h.endswith(".hmm") else h.split("_")[0]
pres={}
for r in allr:
    a=r["accession"]; s1=set(); s2=set()
    for line in open(res(a)):
        f=line.rstrip("\n").split("\t")
        if len(f)<10: continue
        try: ev=float(f[4]); cv=float(f[9])
        except: continue
        if cv<COV: continue
        fam=famof(f[0])
        if fam not in PANEL: continue
        if ev<=E2: s2.add(fam)
        if ev<=E1: s1.add(fam)
    pres[a]={E1:s1,E2:s2}

def prev(g,f,cut):
    n=sum(1 for r in G[g] if f in pres[r["accession"]][cut]); return n,100.0*n/len(G[g])

print("\n"+"="*100)
print("CHITIN PANEL, 15 FAMILIES x 5 GROUPS, uniform filter E<=1e-10 AND cov>=0.30")
print("="*100)
hdr=f"{'family':<7}{'module':<14}{'NOVEL':>8}{'AKKamph':>9}{'Podarcis':>10}{'mammal':>8}{'GTDB':>8}{'amph-vs-ref':>13}{'q1e-9 flip':>11}"
print(hdr); print("-"*100)
res_rows=[]
for f in PANEL:
    n1,p1=prev("NOVEL_amph",f,E1); n2,p2=prev("AKK_amph",f,E1)
    n3,p3=prev("AKK_podarcis",f,E1); n4,p4=prev("AKK_mammal",f,E1); n5,p5=prev("AKK_gtdb",f,E1)
    ref_pos=n4+n5; ref_n=len(G["AKK_mammal"])+len(G["AKK_gtdb"])
    eff=p2-100.0*ref_pos/ref_n
    p=fisher(n2,len(G["AKK_amph"])-n2,ref_pos,ref_n-ref_pos)
    m2,_=prev("AKK_amph",f,E2)
    p9=fisher(m2,len(G["AKK_amph"])-m2,ref_pos,ref_n-ref_pos)
    flip="YES" if (p<0.05)!=(p9<0.05) else ""
    print(f"{f:<7}{MODULE[f]:<14}{p1:>7.1f}%{p2:>8.1f}%{p3:>9.1f}%{p4:>7.1f}%{p5:>7.1f}%{eff:>+12.1f}{flip:>11}")
    res_rows.append([f,MODULE[f],n1,f"{p1:.1f}",n2,f"{p2:.1f}",n3,f"{p3:.1f}",n4,f"{p4:.1f}",
                     n5,f"{p5:.1f}",f"{eff:+.1f}",f"{p:.3g}",flip])
with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["family","module","n_novel","pct_novel","n_akk_amph","pct_akk_amph",
                "n_podarcis","pct_podarcis","n_mammal","pct_mammal","n_gtdb","pct_gtdb",
                "amph_minus_ref","fisher_p","flips_at_1e9"])
    w.writerows(res_rows)
print(f"\nwrote {OUT}")
print("\nREAD IT: 'amph-vs-ref' is amphibian Akkermansia minus (mammal+GTDB) Akkermansia.")
print("  Only GH75 down      -> specific loss.")
print("  All module C down   -> pathway loss.")
print("  Everything down     -> genome reduction.")
print("  GH75 down, GH18 up  -> SUBSTRATE SWITCH.")
print("\nCAVEATS: MAG vs isolate (amph/Podarcis are MAGs, mammal/GTDB largely isolates).")
print("  n is animals not genomes. 1e-10 rejects real GH75 genes (Result 10); flips flagged.")
print("  'Podarcis' is one wall-lizard genus, not reptiles.")
