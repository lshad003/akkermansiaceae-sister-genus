#!/usr/bin/env python3
# Mucin presence matrix assembled
# Source: ch3-chitin-evolution/scripts/build_mucin_matrix.py
# Output: results/novel_akk_tree/mucin_asr_matrix.tsv
import csv, os, collections

BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"
TIPS=f"{BASE}/results/novel_akk_tree/tree_tip_labels.txt"
OUT=f"{BASE}/results/novel_akk_tree/mucin_asr_matrix.tsv"

DIRS=[f"{BASE}/results/dbcan_verru", f"{BASE}/results/dbcan_bact_refs",
      f"{BASE}/results/dbcan_ehi_amphibian", f"{BASE}/results/dbcan_ehi_nonamph",
      f"{BASE}/results/dbcan_endo_allphyla", f"{BASE}/results/dbcan_flavo_refs",
      f"{BASE}/results/dbcan_scaffold", RUM]

MUCIN=sorted("GH2 GH20 GH27 GH29 GH33 GH35 GH84 GH85 GH89 GH95 GH98 GH101 GH109 GH110 GH123 GH129".split())
E_PRIMARY=1e-10
COV=0.30

def fam_of(h):
    return h[:-4].split("_")[0] if h.endswith(".hmm") else h.split("_")[0]

idx={}
for d in DIRS:
    if not os.path.isdir(d):
        print("  MISSING DIR", d); continue
    for fn in os.listdir(d):
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            idx.setdefault(fn[:-4], os.path.join(d,fn))

def resolve(a):
    if a in idx: return idx[a]
    if a[:3] in ("GB_","RS_") and a[3:] in idx: return idx[a[3:]]
    return None

def parse(p):
    fams=set()
    for line in open(p):
        f=line.rstrip("\n").split("\t")
        if len(f)<10: continue
        try:
            ev=float(f[4]); cv=float(f[9])
        except: continue
        if cv<COV: continue
        if ev<=E_PRIMARY: fams.add(fam_of(f[0]))
    return fams

labels=[l.strip() for l in open(TIPS) if l.strip()]
print("n_tips", len(labels))

rows=[]
per_taxon=collections.OrderedDict()
for lab in labels:
    parts=lab.split("|")
    taxon=parts[0]
    acc=parts[-1]
    p=resolve(acc)
    st=per_taxon.setdefault(taxon, {"n":0,"resolved":0,"noann":0,"mucin_sum":0})
    st["n"]+=1
    if p is None:
        rows.append([lab]+["?"]*len(MUCIN))
        st["noann"]+=1
        continue
    fams=parse(p)
    st["resolved"]+=1
    nmuc=sum(1 for m in MUCIN if m in fams)
    st["mucin_sum"]+=nmuc
    rows.append([lab]+["1" if m in fams else "0" for m in MUCIN])

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["tip_label"]+MUCIN)
    for r in rows: w.writerow(r)

print("wrote", OUT)
print()
print("GATE TABLE  (does each taxon actually have annotation?)")
print(f"{'taxon':<22}{'n_tips':>8}{'resolved':>10}{'no_ann':>8}{'mean_mucin_fams':>18}")
for tx,st in sorted(per_taxon.items(), key=lambda kv: -kv[1]['n']):
    mm = st["mucin_sum"]/st["resolved"] if st["resolved"] else float("nan")
    print(f"{tx:<22}{st['n']:>8}{st['resolved']:>10}{st['noann']:>8}{mm:>18.2f}")
