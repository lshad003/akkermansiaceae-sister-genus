#!/bin/bash -l
#SBATCH -p epyc -c 16 --mem 48G -t 08:00:00
#SBATCH -J mucinpfam -o /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/mucinpfam_%j.out

CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/mucin_pfam
DB=$CH3/db/pfam37/panel.hmm
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
mkdir -p $W

module unload hmmer 2>/dev/null
module load hmmer/3.3.2
hmmscan -h | head -2

$PY - << 'PYEOF'
import os, csv
CH3="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
REPO="/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus"
W=CH3+"/results/mucin_pfam"
CEN=CH3+"/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"

DIRS=[CH3+"/results/pangenome/"+d for d in ("gff199","akkfam","outgroups","gtdb_akk")]
DIRS+=[CH3+"/results/"+d for d in ("dbcan_scaffold","dbcan_ehi_amphibian","dbcan_ehi_nonamph")]
DIRS+=["/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/bin_proteins"]
idx={}
for d in DIRS:
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d,f))

cand=set(r["genome"] for r in csv.DictReader(open(REPO+"/tables/TableS1_genome_quality.tsv"),delimiter="\t"))
FREE=("Luteolibacter","Haloferula","Rubritalea","Roseibacillus","Roseibacillus_B",
      "SW10","Oceaniferula","UBA956","Persicirhabdus","WTJZ01")
rows=list(csv.DictReader(open(CEN),delimiter="\t"))
sel=[]
for r in rows:
    a=r["accession"]
    if a in cand: grp="candidate"
    elif r["genus"]=="Akkermansia": grp="akkermansia"
    elif r["genus"] in FREE: grp="free"
    else: continue
    p=idx.get(a) or idx.get(a.replace("GB_","").replace("RS_",""))
    if p: sel.append((grp,a,p))
from collections import Counter
print("selected:", dict(Counter(g for g,_,_ in sel)))
with open(W+"/members.tsv","w") as o:
    o.write("group\tgenome\tpath\n")
    for g,a,p in sel: o.write("%s\t%s\t%s\n"%(g,a,p))
PYEOF
[ -s $W/members.tsv ] || { echo "STAGING FAILED"; exit 1; }

: > $W/all.faa
while IFS=$'\t' read -r G A P; do
  [ "$G" = "group" ] && continue
  sed "s/^>/>${A}|/" "$P" >> $W/all.faa
done < $W/members.tsv
echo "proteins: $(grep -c '^>' $W/all.faa)"

hmmscan --domtblout $W/panel.domtbl -E 1e-5 --cpu 16 -o /dev/null $DB $W/all.faa
echo "domtbl lines: $(grep -vc '^#' $W/panel.domtbl)"
rm -f $W/all.faa
