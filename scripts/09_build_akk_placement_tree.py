#!/usr/bin/env python3
# Placement alignment extracted
# Source: ch3-chitin-evolution/scripts/build_akk_placement_tree.py
# Output: results/novel_akk_tree/akk_placement.faa
import gzip, csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
MSA=f"{BASE}/results/gtdbtk_amphibia/align/align/gtdbtk.bac120.msa.fasta.gz"
CEN=f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
OUT=f"{BASE}/results/novel_akk_tree"; os.makedirs(OUT,exist_ok=True)
cen={r["accession"]:r for r in csv.DictReader(open(CEN),delimiter="\t")}
def g(r): return (r["genus"] or "").strip()

want={}
for a,r in cen.items():
    fam=r["family"]; gen=g(r); hc=r["host_class"]; ds=r["from_dataset"]
    if fam=="Akkermansiaceae":
        if hc=="amphibian" and gen in ("","unknown","NO_GENUS"):
            want[a]=f"NOVEL|{ds}|{r.get('host_animal_type','NA')}|GH75_{r['gh75_present']}"
        elif gen=="Akkermansia":
            want[a]=f"Akkermansia|{ds}|{hc}|GH75_{r['gh75_present']}"
        elif gen:
            want[a]=f"{gen}|{ds}|{hc}|GH75_{r['gh75_present']}"
    elif fam in ("Verrucomicrobiaceae","Chthoniobacteraceae") and ds=="GTDB_r226_rep":
        want[a]=f"OUTGROUP_{fam}|{gen or 'NA'}|{hc}|GH75_{r['gh75_present']}"

alt={}
for a in want:
    alt[a]=a
    if a[:3] in ("GB_","RS_"): alt[a[3:]]=a
print("targets:",len(want))

recs={}
with gzip.open(MSA,"rt") as fh:
    keep=None; buf=[]
    for line in fh:
        if line.startswith(">"):
            if keep: recs[keep]="".join(buf)
            i=line[1:].split()[0]
            keep=alt.get(i); buf=[]
        elif keep: buf.append(line.strip())
    if keep: recs[keep]="".join(buf)
print("found in MSA:",len(recs),"/",len(want))

c=collections.Counter(want[a].split("|")[0] for a in recs)
for k,v in c.most_common(): print(f"  {k:<32}{v}")

fa=f"{OUT}/akk_placement.faa"
with open(fa,"w") as o:
    for a,s in recs.items():
        o.write(f">{want[a]}|{a}\n")
        for i in range(0,len(s),80): o.write(s[i:i+80]+"\n")
print("\nwrote",fa,len(recs),"seqs, aln",len(next(iter(recs.values()))))
