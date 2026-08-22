#!/usr/bin/env python3
# Genus calls resolved for amphibian Akkermansiaceae
# Source: ch3-chitin-evolution/scripts/resolve_akk_genus.py
# Output: stdout
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
CANDS=[f"{BASE}/data/herptile_gtdbtk_compiled.tsv",
       f"{BASE}/data/herptile_cazyme_taxonomy_joined.tsv"]

for p in CANDS:
    print("="*70)
    print(p, "EXISTS" if os.path.exists(p) else "MISSING")
    if not os.path.exists(p): continue
    with open(p) as fh:
        hdr=fh.readline().rstrip("\n").split("\t")
        row1=fh.readline().rstrip("\n").split("\t")
    print("  cols:", len(hdr))
    for i,(h,v) in enumerate(zip(hdr,row1),1):
        v=(v[:60]+"...") if len(v)>60 else v
        print(f"    {i:>3} {h:<26} {v}")

rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"]
need=[r for r in rows if r["host_class"]=="amphibian"
      and (r["genus"] or "").strip() in ("","unknown","NO_GENUS")]
print("\n"+"="*70)
print(f"amphibian Akkermansiaceae needing a genus call: {len(need)}")
print("  examples:", [r["accession"] for r in need[:3]])

def parse_g(s):
    for tok in (s or "").split(";"):
        tok=tok.strip()
        if tok.startswith("g__"): return tok[3:] or "UNCLASSIFIED_AT_GENUS"
    return None

for p in CANDS:
    if not os.path.exists(p): continue
    print("\n"+"-"*70)
    print("TRY:",p)
    with open(p) as fh: hdr=fh.readline().rstrip("\n").split("\t")
    idcol=next((c for c in ("bin_id","user_genome","accession","genome","Name") if c in hdr), None)
    gcol =next((c for c in ("genus","g","Genus") if c in hdr), None)
    ccol =next((c for c in ("gtdb_classification","classification","Classification") if c in hdr), None)
    print(f"  id col={idcol}  genus col={gcol}  classification col={ccol}")
    if not idcol: print("  NO ID COLUMN, skipping"); continue
    lut={}
    for r in csv.DictReader(open(p),delimiter="\t"):
        k=(r[idcol] or "").strip()
        g=None
        if gcol and (r.get(gcol) or "").strip(): g=r[gcol].strip()
        if not g and ccol: g=parse_g(r.get(ccol))
        if k and g: lut[k]=g
    print(f"  usable genus calls in table: {len(lut)}")
    hit=collections.Counter(); missn=0
    for r in need:
        g=lut.get(r["accession"])
        if g is None: missn+=1
        else: hit[g]+=1
    print(f"  RESOLVED {len(need)-missn}/{len(need)}   still missing {missn}")
    for g,n in hit.most_common(12): print(f"    {g:<30} {n}")
