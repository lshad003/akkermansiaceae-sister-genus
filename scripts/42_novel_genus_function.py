#!/usr/bin/env python3
# Mucin family prevalence across groups
# Source: ch3-chitin-evolution/scripts/novel_genus_function.py
# Output: results/novel_akk_tree/novel_genus_function.tsv
import csv, os, math, random, collections

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
OUTD=f"{BASE}/results/novel_akk_tree"
OUT=f"{OUTD}/novel_genus_function.tsv"
os.makedirs(OUTD, exist_ok=True)

DIRS=[f"{BASE}/results/dbcan_verru", f"{BASE}/results/dbcan_bact_refs",
      f"{BASE}/results/dbcan_ehi_amphibian", f"{BASE}/results/dbcan_ehi_nonamph",
      f"{BASE}/results/dbcan_endo_allphyla", f"{BASE}/results/dbcan_flavo_refs",
      f"{BASE}/results/dbcan_scaffold", RUM]
MUCIN=set("GH2 GH20 GH27 GH29 GH33 GH35 GH84 GH85 GH89 GH95 GH98 GH101 GH109 GH110 GH123 GH129".split())
E_PRIMARY=1e-10
E_ALT=1e-9
COV=0.30
FLOOR=0.10
DRAWS=499
SEED=20260715
BAR="="*78

def logf(n):
    return math.lgamma(n+1)
def hyp_logp(a,b,c,d):
    n=a+b+c+d
    return (logf(a+b)+logf(c+d)+logf(a+c)+logf(b+d)-logf(n)
            -logf(a)-logf(b)-logf(c)-logf(d))
def fisher2(a,b,c,d):
    obs=hyp_logp(a,b,c,d)
    r1=a+b; r2=c+d; c1=a+c; n=r1+r2
    lo=max(0,c1-r2); hi=min(r1,c1)
    tot=0.0
    for x in range(lo,hi+1):
        lp=hyp_logp(x,r1-x,c1-x,r2-c1+x)
        if lp<=obs+1e-9: tot+=math.exp(lp)
    return min(1.0,tot)
def bh(pairs):
    m=len(pairs)
    s=sorted(pairs,key=lambda t:t[1])
    q={}; prev=1.0
    for i in range(m-1,-1,-1):
        k,p=s[i]
        v=min(prev,p*m/(i+1))
        q[k]=v; prev=v
    return q
def median(v):
    if not v: return float("nan")
    v=sorted(v); n=len(v)
    return v[n//2] if n%2 else (v[n//2-1]+v[n//2])/2.0

print(BAR); print("NOVEL GENUS FUNCTION: genus-matched CAZyme repertoire"); print(BAR)

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()
def isnov(r): return gen(r) in ("","unknown","NO_GENUS")

G=collections.OrderedDict()
G["NOVEL"]=[r for r in rows if r["host_class"]=="amphibian" and isnov(r)]
G["AKK_AMPH"]=[r for r in rows if r["host_class"]=="amphibian" and gen(r)=="Akkermansia"]
G["AKK_MAMMAL"]=[r for r in rows if r["host_class"]=="mammal" and gen(r)=="Akkermansia"]
G["AKK_GTDB"]=[r for r in rows if r["host_class"]=="unknown" and gen(r)=="Akkermansia"]
G["AKK_PODARCIS"]=[r for r in rows if r["host_class"]=="reptile" and gen(r)=="Akkermansia"]
LABEL={"NOVEL":"NOVEL (unnamed genus)","AKK_AMPH":"Akkermansia amphibian",
       "AKK_MAMMAL":"Akkermansia mammal","AKK_GTDB":"Akkermansia GTDB",
       "AKK_PODARCIS":"Akkermansia Podarcis (wall lizard)"}
for k,v in G.items(): print(f"  {LABEL[k]:<38} n={len(v)}")
if len(G["NOVEL"])!=105 or len(G["AKK_AMPH"])!=94:
    print(f"GUARD FAILED: NOVEL={len(G['NOVEL'])} expected 105, AKK_AMPH={len(G['AKK_AMPH'])} expected 94")
    raise SystemExit(1)
print("  guard passed")

idx={}
for d in DIRS:
    if not os.path.isdir(d): print("  MISSING DIR",d); continue
    for fn in os.listdir(d):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            idx.setdefault(fn[:-4], os.path.join(d,fn))
def resolve(a):
    if a in idx: return idx[a]
    if a[:3] in ("GB_","RS_") and a[3:] in idx: return idx[a[3:]]
    return None

allr=[r for v in G.values() for r in v]
unres=[r["accession"] for r in allr if resolve(r["accession"]) is None]
print(f"  resolved {len(allr)-len(unres)}/{len(allr)} tsvs")
if unres:
    print("  UNRESOLVED:",unres[:10]); raise SystemExit(1)

def fam_of(h):
    return h[:-4].split("_")[0] if h.endswith(".hmm") else h.split("_")[0]

pres={}   # accession -> {eval_cut: set(families)}
npass=nfail=0
raw=set()
for r in allr:
    a=r["accession"]; p=resolve(a)
    s1=set(); s2=set()
    for line in open(p):
        f=line.rstrip("\n").split("\t")
        if len(f)<10: continue
        try:
            ev=float(f[4]); cv=float(f[9])
        except: continue
        raw.add(f[0])
        if cv<COV:
            nfail+=1; continue
        fam=fam_of(f[0])
        if ev<=E_ALT: s2.add(fam)
        if ev<=E_PRIMARY: s1.add(fam); npass+=1
        else: nfail+=1
    pres[a]={E_PRIMARY:s1, E_ALT:s2}
print(f"  domain lines: passed {npass}  failed {nfail}")
print(f"  raw HMM names {len(raw)} collapsed to {len(set(fam_of(x) for x in raw))} families")

def animal(r):
    return r["accession"].split(".")[0] if r["from_dataset"]=="herptile_MAG" else (r["sample_id"] or "?")
ANI={g:collections.defaultdict(list) for g in ("NOVEL","AKK_AMPH")}
for g in ("NOVEL","AKK_AMPH"):
    for r in G[g]: ANI[g][animal(r)].append(r["accession"])
print(f"  NOVEL: {len(G['NOVEL'])} genomes / {len(ANI['NOVEL'])} animals")
print(f"  AKK_AMPH: {len(G['AKK_AMPH'])} genomes / {len(ANI['AKK_AMPH'])} animals")

fams=set()
for g,v in G.items():
    cnt=collections.Counter()
    for r in v:
        for f in pres[r["accession"]][E_PRIMARY]: cnt[f]+=1
    for f,n in cnt.items():
        if n/len(v)>=FLOOR: fams.add(f)
fams=sorted(fams)
print(f"  families clearing the {int(FLOOR*100)}% floor in >=1 group: {len(fams)}")

def prev(g,f,cut):
    n=sum(1 for r in G[g] if f in pres[r["accession"]][cut])
    return n, 100.0*n/len(G[g])

def rarefy(g,cut):
    rnd=random.Random(SEED)
    animals=sorted(ANI[g].keys())
    acc={f:[] for f in fams}
    for _ in range(DRAWS):
        pick=[rnd.choice(ANI[g][an]) for an in animals]
        c=collections.Counter()
        for a in pick:
            for f in pres[a][cut]:
                if f in acc: c[f]+=1
        for f in fams: acc[f].append(100.0*c[f]/len(animals))
    return {f:median(v) for f,v in acc.items()}
RAR={g:{cut:rarefy(g,cut) for cut in (E_PRIMARY,E_ALT)} for g in ("NOVEL","AKK_AMPH")}

def contrast(gA,gB_list,cut):
    eff={}; pv={}
    for f in fams:
        a=sum(1 for r in G[gA] if f in pres[r["accession"]][cut]); b=len(G[gA])-a
        c=sum(1 for g in gB_list for r in G[g] if f in pres[r["accession"]][cut])
        d=sum(len(G[g]) for g in gB_list)-c
        eff[f]=100.0*a/(a+b)-100.0*c/(c+d)
        pv[f]=fisher2(a,b,c,d)
    q=bh(list(pv.items()))
    return eff,pv,q

A_eff,A_p,A_q=contrast("NOVEL",["AKK_AMPH"],E_PRIMARY)
B_eff,B_p,B_q=contrast("AKK_AMPH",["AKK_GTDB","AKK_MAMMAL"],E_PRIMARY)
A_eff9,A_p9,A_q9=contrast("NOVEL",["AKK_AMPH"],E_ALT)
B_eff9,B_p9,B_q9=contrast("AKK_AMPH",["AKK_GTDB","AKK_MAMMAL"],E_ALT)
A_flip={f:((A_q[f]<0.05)!=(A_q9[f]<0.05)) for f in fams}
B_flip={f:((B_q[f]<0.05)!=(B_q9[f]<0.05)) for f in fams}
print(f"  Contrast A families flipping BH-significance at 1e-9: {sum(A_flip.values())}/{len(fams)}")
print(f"  Contrast B families flipping BH-significance at 1e-9: {sum(B_flip.values())}/{len(fams)}")

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["family","in_mucin_panel","n_novel","pct_novel","pct_novel_rarefied",
                "n_akk_amph","pct_akk_amph","pct_akk_amph_rarefied","n_akk_mammal","pct_akk_mammal",
                "n_akk_gtdb","pct_akk_gtdb","n_akk_podarcis","pct_akk_podarcis",
                "A_effect","A_fisher_p","A_bh_q","A_flips_at_1e9",
                "B_effect","B_fisher_p","B_bh_q","B_flips_at_1e9"])
    for f in fams:
        nn,pn=prev("NOVEL",f,E_PRIMARY); na,pa=prev("AKK_AMPH",f,E_PRIMARY)
        nm,pm=prev("AKK_MAMMAL",f,E_PRIMARY); ng,pg=prev("AKK_GTDB",f,E_PRIMARY)
        np_,pp=prev("AKK_PODARCIS",f,E_PRIMARY)
        w.writerow([f,int(f in MUCIN),nn,f"{pn:.1f}",f"{RAR['NOVEL'][E_PRIMARY][f]:.1f}",
                    na,f"{pa:.1f}",f"{RAR['AKK_AMPH'][E_PRIMARY][f]:.1f}",nm,f"{pm:.1f}",
                    ng,f"{pg:.1f}",np_,f"{pp:.1f}",
                    f"{A_eff[f]:+.1f}",f"{A_p[f]:.3g}",f"{A_q[f]:.3g}",int(A_flip[f]),
                    f"{B_eff[f]:+.1f}",f"{B_p[f]:.3g}",f"{B_q[f]:.3g}",int(B_flip[f])])

print("\n"+BAR); print("1. MUCIN PANEL (pre-declared). Does NOVEL retain the Akkermansia toolkit?"); print(BAR)
print(f"  {'family':<9}{'NOVEL':>8}{'rar':>7}{'AKK_am':>8}{'rar':>7}{'AKK_mam':>9}{'AKK_gtdb':>10}{'Podarcis':>10}{'A_eff':>8}{'A_q':>10}  flip")
retained=absent=0
for f in sorted(MUCIN):
    if f not in fams:
        nn,pn=prev("NOVEL",f,E_PRIMARY); na,pa=prev("AKK_AMPH",f,E_PRIMARY)
        print(f"  {f:<9}{pn:>7.1f}%{'':>7}{pa:>7.1f}%{'':>7}{'':>9}{'':>10}{'':>10}{'':>8}{'below floor':>10}")
        continue
    nn,pn=prev("NOVEL",f,E_PRIMARY); na,pa=prev("AKK_AMPH",f,E_PRIMARY)
    nm,pm=prev("AKK_MAMMAL",f,E_PRIMARY); ng,pg=prev("AKK_GTDB",f,E_PRIMARY)
    np_,pp=prev("AKK_PODARCIS",f,E_PRIMARY)
    print(f"  {f:<9}{pn:>7.1f}%{RAR['NOVEL'][E_PRIMARY][f]:>7.1f}{pa:>7.1f}%{RAR['AKK_AMPH'][E_PRIMARY][f]:>7.1f}"
          f"{pm:>8.1f}%{pg:>9.1f}%{pp:>9.1f}%{A_eff[f]:>+8.1f}{A_q[f]:>10.2g}  {'YES' if A_flip[f] else ''}")
    if pn>=50: retained+=1
    if pn<10: absent+=1
inpanel=[f for f in sorted(MUCIN) if f in fams]
print(f"\n  mucin families clearing the floor: {len(inpanel)}/{len(MUCIN)}")
print(f"  present in >=50% of NOVEL: {retained}    present in <10% of NOVEL: {absent}")
if inpanel:
    dn=median([prev("NOVEL",f,E_PRIMARY)[1] for f in inpanel])
    da=median([prev("AKK_AMPH",f,E_PRIMARY)[1] for f in inpanel])
    print(f"  median mucin-panel prevalence: NOVEL {dn:.1f}%  vs  amphibian Akkermansia {da:.1f}%")
    verdict = "RETAINED" if dn>=50 else ("PARTIAL" if dn>=20 else "NOT RETAINED")
    print(f"  VERDICT: NOVEL {verdict} the Akkermansia mucin toolkit")
else:
    verdict="UNTESTABLE"
    print("  VERDICT: UNTESTABLE, no mucin family cleared the floor")

def famline(f,eff,p,q,flip):
    nn,pn=prev("NOVEL",f,E_PRIMARY); na,pa=prev("AKK_AMPH",f,E_PRIMARY)
    tag=" [MUCIN]" if f in MUCIN else ""
    return f"  {f:<10}{eff[f]:>+7.1f}  NOVEL {pn:>5.1f}%  AKK_amph {pa:>5.1f}%  q={q[f]:.2g}{tag}"
stableA=[f for f in fams if not A_flip[f]]
sa=sorted(stableA,key=lambda f:A_eff[f])
print("\n"+BAR); print("2. CONTRAST A: NOVEL vs amphibian Akkermansia (genus effect, host held constant)"); print(BAR)
print("\n  TOP 20 NOVEL-ENRICHED:")
for f in sa[::-1][:20]: print(famline(f,A_eff,A_p,A_q,A_flip))
print("\n  TOP 20 NOVEL-DEPLETED:")
for f in sa[:20]: print(famline(f,A_eff,A_p,A_q,A_flip))
print(f"\n  ({sum(A_flip.values())} families excluded for flipping between 1e-10 and 1e-9)")

stableB=[f for f in fams if not B_flip[f]]
sb=sorted(stableB,key=lambda f:B_eff[f])
def famlineB(f):
    na,pa=prev("AKK_AMPH",f,E_PRIMARY); ng,pg=prev("AKK_GTDB",f,E_PRIMARY)
    tag=" [MUCIN]" if f in MUCIN else ""
    return f"  {f:<10}{B_eff[f]:>+7.1f}  AKK_amph {pa:>5.1f}%  ref {pg:>5.1f}%  q={B_q[f]:.2g}{tag}"
print("\n"+BAR); print("3. CONTRAST B: amphibian Akkermansia vs reference Akkermansia (host effect, genus held constant)"); print(BAR)
print("\n  TOP 15 amphibian-ENRICHED:")
for f in sb[::-1][:15]: print(famlineB(f))
print("\n  TOP 15 amphibian-DEPLETED:")
for f in sb[:15]: print(famlineB(f))
print(f"\n  ({sum(B_flip.values())} families excluded for flipping)")

print("\n"+BAR); print("4. Median completeness and median retained-family count per genome"); print(BAR)
for g in G:
    comps=[]
    for r in G[g]:
        try: comps.append(float(r["completeness"]))
        except: pass
    reps=[sum(1 for f in fams if f in pres[r["accession"]][E_PRIMARY]) for r in G[g]]
    print(f"  {LABEL[g]:<38} n={len(G[g]):>4}  completeness {median(comps):>5.1f}%  retained families/genome {median(reps):>4.0f}")

print("\n"+BAR); print("5. GH75 by group (CONTEXT ONLY, not the question here)"); print(BAR)
for g in G:
    n,p=prev(g,"GH75",E_PRIMARY)
    print(f"  {LABEL[g]:<38} GH75 {p:>5.1f}%  ({n}/{len(G[g])})")
if "GH75" in fams:
    print(f"  Contrast A: effect {A_eff['GH75']:+.1f}  p={A_p['GH75']:.2g}  q={A_q['GH75']:.2g}  "
          f"flip@1e-9={'yes' if A_flip['GH75'] else 'no'}")
else:
    print("  GH75 did not clear the prevalence floor")

print("\n"+BAR); print("6. CAVEATS (verbatim, do not soften)"); print(BAR)
print("   - Two different genera differ in many families by definition. Contrast A is a")
print("     CHARACTERIZATION of an unnamed genus, not a hypothesis test, and a long significant")
print("     list is expected, not a finding.")
print("   - NOVEL is amphibian-only, so genus and host are perfectly confounded within NOVEL.")
print("     Contrast A controls host by comparing to amphibian Akkermansia from the same datasets.")
print("     Contrast B is where a host effect could appear, and its reference is largely isolates.")
print(f"   - n is ANIMALS not genomes: NOVEL {len(ANI['NOVEL'])} animals / {len(G['NOVEL'])} genomes; "
      f"AKK_AMPH {len(ANI['AKK_AMPH'])} animals / {len(G['AKK_AMPH'])} genomes.")
print("   - MAG vs isolate: NOVEL and AKK_AMPH are MAGs; AKK_GTDB and AKK_MAMMAL are largely")
print("     isolates. Contrast A is MAG vs MAG and so is not exposed to this; Contrast B is.")
print("   - Threshold: 1e-10 is pre-registered but CHATINDEX Result 10 shows it rejects real genes.")
print("     Families flagged as flipping at 1e-9 are not reportable.")
print("   - This says what genes are PRESENT. It does not show the lineage degrades mucin in vivo.")
print(f"\nwrote {OUT}")
print("\nSuggested PROJECT_LOG line (do NOT append it here):")
print(f"  2026-07-15 | novel genus function (genus-matched, replaces Result 11) | "
      f"NOVEL vs amphibian Akkermansia: mucin toolkit {verdict}, {len(stableA)} stable retained "
      f"families, GH75 context printed; wrote novel_genus_function.tsv")
