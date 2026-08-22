#!/usr/bin/env python3
# Genome-quality explanation tested against assembly type
# Source: ch3-chitin-evolution/scripts/check_mag_vs_isolate.py
# Output: stdout
import csv, collections, statistics
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
FREE=["Luteolibacter","Haloferula","Rubritalea","Roseibacillus_B","SW10","Oceaniferula","UBA956"]

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t") if r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()
G=collections.OrderedDict()
G["NOVEL"]=[r for r in rows if r["family"]=="Akkermansiaceae" and r["host_class"]=="amphibian" and gen(r) in ("","unknown","NO_GENUS")]
G["AKK_all"]=[r for r in rows if r["family"]=="Akkermansiaceae" and gen(r)=="Akkermansia"]
for g in FREE:
    v=[r for r in rows if gen(r)==g]
    if len(v)>=10: G[g]=v

print("%-16s %5s  %-28s  %s" % ("group","n","genome_type","completeness median (IQR)"))
for k,v in G.items():
    t=collections.Counter((r.get("genome_type") or "?") for r in v)
    comps=[]
    for r in v:
        try: comps.append(float(r["completeness"]))
        except: pass
    if comps:
        q=statistics.quantiles(comps,n=4) if len(comps)>3 else [min(comps),statistics.median(comps),max(comps)]
        s="%.1f (%.1f-%.1f)"%(statistics.median(comps),q[0],q[2])
    else: s="n/a"
    print("%-16s %5d  %-28s  %s" % (k,len(v),dict(t.most_common()),s))

print("\nPOOLED:")
gut=G["NOVEL"]+G["AKK_all"]
free=[r for k in FREE if k in G for r in G[k]]
for lab,v in (("gut (novel+Akk)",gut),("free-living",free)):
    t=collections.Counter((r.get("genome_type") or "?") for r in v)
    comps=[float(r["completeness"]) for r in v if r["completeness"]]
    iso=100.0*t.get("isolate",0)/len(v)
    print("  %-18s n=%4d  isolate %.1f%%  MAG %.1f%%  median completeness %.1f" % (
        lab,len(v),iso,100.0*t.get("MAG",0)/len(v),statistics.median(comps)))
