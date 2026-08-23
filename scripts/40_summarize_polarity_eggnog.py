#!/usr/bin/env python3
# Annotated losses summarized by branch
# Source: ch3-chitin-evolution/scripts/summarize_polarity_eggnog.py
# Output: results/pangenome/polarity_annotated.tsv
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
D=f"{BASE}/results/pangenome/eggnog_polarity"
OUT=f"{BASE}/results/pangenome/polarity_annotated.tsv"
COG={"J":"translation","A":"RNA processing","K":"transcription","L":"replication/repair",
"B":"chromatin","D":"cell cycle","Y":"nuclear","V":"defense","T":"signal transduction",
"M":"cell wall/membrane","N":"motility","Z":"cytoskeleton","W":"extracellular",
"U":"trafficking/secretion","O":"post-translational/chaperone","X":"mobilome",
"C":"energy production","G":"carbohydrate metabolism","E":"amino acid metabolism",
"F":"nucleotide metabolism","H":"coenzyme metabolism","I":"lipid metabolism",
"P":"inorganic ion transport","Q":"secondary metabolites","R":"general function",
"S":"unknown function"}

rows=[]
for C in ("A_lost_in_Akkermansia","C_Akkermansia_enriched","D_novel_specific"):
    p=f"{D}/{C}.emapper.annotations"
    if not os.path.isfile(p):
        print("MISSING",p); continue
    hdr=None
    n=0
    cogs=collections.Counter(); named=[]
    for line in open(p):
        if line.startswith("##"): continue
        if line.startswith("#"):
            hdr=line[1:].rstrip("\n").split("\t"); continue
        f=line.rstrip("\n").split("\t")
        if not hdr or len(f)<len(hdr)-2: continue
        d=dict(zip(hdr,f)); n+=1
        cat=d.get("COG_category","-") or "-"
        for ch in cat:
            if ch in COG: cogs[ch]+=1
        gname=d.get("Preferred_name","-")
        desc=d.get("Description","-")
        kegg=d.get("KEGG_ko","-")
        rows.append([C,d.get("query","?"),cat,gname,kegg,desc[:90]])
        if gname and gname!="-": named.append((gname,desc[:70]))
    print(f"\n=== {C}: {n} annotated ===")
    print("  top COG categories:")
    for ch,c in cogs.most_common(8):
        print(f"    {ch} {COG[ch]:<32} {c}")
    print("  named genes:")
    for g,ds in named[:25]:
        print(f"    {g:<12} {ds}")

with open(OUT,"w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["category","orthogroup","COG","gene","KEGG_ko","description"])
    w.writerows(rows)
print("\nwrote",OUT)
