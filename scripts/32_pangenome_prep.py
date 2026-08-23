#!/usr/bin/env python3
# Group assignments and pooled proteins assembled
# Source: ch3-chitin-evolution/scripts/pangenome_prep.py
# Output: results/pangenome/genome_groups.tsv, results/pangenome/all.faa
import csv, os, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN=f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
OUT=f"{BASE}/results/pangenome"; os.makedirs(OUT,exist_ok=True)
FAA_DIRS=[f"{BASE}/results/dbcan_ehi_amphibian",
          f"{BASE}/results/dbcan_ehi_nonamph",
          "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/bin_proteins",
          "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/prodigal"]
idx={}
for d in FAA_DIRS:
    if not os.path.isdir(d): print("MISSING",d); continue
    n=0
    for f in os.listdir(d):
        if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d,f)); n+=1
    print(f"  {d}  {n} faa")
rows=[r for r in csv.DictReader(open(CEN),delimiter="\t")
      if r["family"]=="Akkermansiaceae" and r["annotated"]=="1"]
def gen(r): return (r["genus"] or "").strip()
G={"NOVEL":[r for r in rows if r["host_class"]=="amphibian" and gen(r) in ("","unknown","NO_GENUS")],
   "AKK_AMPH":[r for r in rows if r["host_class"]=="amphibian" and gen(r)=="Akkermansia"]}
print()
for k,v in G.items(): print(f"{k}: {len(v)} genomes")
tot=miss=0
with open(f"{OUT}/genome_groups.tsv","w") as gt, open(f"{OUT}/all.faa","w") as fa:
    gt.write("genome\tgroup\thost_animal\tdataset\n")
    for grp,v in G.items():
        for r in v:
            a=r["accession"]; p=idx.get(a)
            if not p: miss+=1; print("  NO FAA:",grp,a); continue
            an=a.split(".")[0] if r["from_dataset"]=="herptile_MAG" else (r["sample_id"] or "?")
            gt.write(f"{a}\t{grp}\t{an}\t{r['from_dataset']}\n")
            n=0
            for line in open(p):
                if line.startswith(">"):
                    n+=1; fa.write(f">{a}|prot{n}\n")
                else: fa.write(line)
            tot+=1
print(f"\ngenomes written {tot}, missing faa {miss}")
print(f"wrote {OUT}/all.faa and {OUT}/genome_groups.tsv")
