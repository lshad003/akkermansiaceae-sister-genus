#!/bin/bash
#SBATCH --job-name=sgl
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/sgl_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/sgl_%j.err
module load diamond
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/sgl_prevalence
mkdir -p $W
# query: the SGL protein from the type genome
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import sys
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
src=B+"/data/UHM_Akkermansia.filtered.faa"
tgt="UHM1073.23039_R.bin.125|k141_311317_53"
keep=False; out=[]
for line in open(src):
    if line.startswith(">"):
        keep = line[1:].split()[0] == tgt
        if keep: out.append(">SGL_query\n")
    elif keep: out.append(line)
if not out:
    print("REFUSED: query not found"); sys.exit(1)
open(B+"/results/sgl_prevalence/query.faa","w").writelines(out)
print("query length: %d aa" % sum(len(l.strip()) for l in out[1:]))
PY
# search the same unified 7,572-genome database used for zwf, gnd and opcA
U=$CH3/results/pgl_unified
M=$CH3/results/ppp_unified/manifest.tsv
: > $W/all.faa
while IFS=$'\t' read -r B2 DIR; do
  F="$DIR/$B2.faa"
  [ -s "$F" ] && sed "s/^>/>${B2}|/" "$F" >> $W/all.faa
done < "$M"
echo "proteins: $(grep -c '^>' $W/all.faa)"
diamond makedb --in $W/all.faa -d $W/db --quiet
diamond blastp -q $W/query.faa -d $W/db -o $W/sgl_hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 \
  --max-target-seqs 100000 --more-sensitive --threads 16 --quiet
echo "hits: $(wc -l < $W/sgl_hits.tsv)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import csv
from collections import Counter
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/sgl_prevalence"
cen=list(csv.DictReader(open(B+"/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"),delimiter="\t"))
grp,gen={},{}
for r in cen:
    if r["family"]!="Akkermansiaceae": continue
    a=r["accession"]; gen[a]=r["genus"]
    if r["genus"]=="Akkermansia": grp[a]="Akkermansia"
    elif not r["genus"].strip() or r["genus"] in ("unknown","NO_GENUS"):
        if r["host_class"]=="amphibian": grp[a]="candidate"
    else: grp[a]="free-living"
man=[l.split("\t")[0].strip() for l in open(B+"/results/ppp_unified/manifest.tsv") if l.strip()]
tot=Counter(grp[g] for g in man if g in grp)
for floor in (0.0,30.0):
    hit={}
    for line in open(W+"/sgl_hits.tsv"):
        c=line.split("\t")
        if float(c[2])<floor: continue
        g=c[1].split("|",1)[0].split("::")[-1]
        if g in grp: hit.setdefault(grp[g],set()).add(g)
    print("floor %2.0f%%:"%floor, {k:"%d/%d"%(len(hit.get(k,())),tot[k]) for k in ("candidate","Akkermansia","free-living")},
          " free-living genera:", len(set(gen[g] for g in hit.get("free-living",()))))
PY
rm -f $W/db.dmnd $W/all.faa
