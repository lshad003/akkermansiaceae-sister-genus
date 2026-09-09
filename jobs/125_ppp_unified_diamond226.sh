#!/bin/bash
#SBATCH --job-name=pppD2
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pppD2_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pppD2_%j.err
module load diamond/2.2.6
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/ppp_unified_d2
mkdir -p $W
D=$(which diamond); echo "diamond: $D"; $D --version
M=$CH3/results/ppp_unified/manifest.tsv
N=$(wc -l < $M); [ "$N" -eq 7572 ] || { echo "REFUSED: manifest has $N rows"; exit 1; }
: > $W/all.faa
while IFS=$'\t' read -r B DIR; do
  F="$DIR/$B.faa"
  [ -s "$F" ] && sed "s/^>/>${B}|/" "$F" >> $W/all.faa
done < "$M"
echo "proteins: $(grep -c '^>' $W/all.faa)"
diamond makedb --in $W/all.faa -d $W/db --threads 16 --quiet
PG=$CH3/results/pangenome
for Q in ppp_trio control_gh20; do
  diamond blastp -q $PG/$Q.faa -d $W/db -o $W/$Q.tsv \
    --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 \
    --max-target-seqs 100000 --more-sensitive --threads 16 --quiet
  echo "$Q hits: $(wc -l < $W/$Q.tsv)"
done
for Q in pgl sgl; do
  S=$CH3/results/pgl_search/${Q}_query.faa
  [ -s "$S" ] || S=$CH3/results/sgl_prevalence/query.faa
  [ -s "$S" ] || continue
  diamond blastp -q $S -d $W/db -o $W/${Q}.tsv \
    --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 \
    --max-target-seqs 100000 --more-sensitive --threads 16 --quiet
  echo "$Q hits: $(wc -l < $W/${Q}.tsv)"
done
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import csv, os
from collections import Counter, defaultdict
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/ppp_unified_d2"
QN={"UHM979.41089_R.bin.103_CDS_0654":"zwf",
    "UHM1210.23070_R.bin.101_CDS_0293":"gnd",
    "EHM058980_CDS_1432":"opcA"}
cen=list(csv.DictReader(open(B+"/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"),delimiter="\t"))
grp={}
for r in cen:
    if r["family"]!="Akkermansiaceae": continue
    a=r["accession"]
    if r["genus"]=="Akkermansia": grp[a]="Akkermansia"
    elif not r["genus"].strip() or r["genus"] in ("unknown","NO_GENUS"):
        if r["host_class"]=="amphibian": grp[a]="candidate"
    else: grp[a]="free-living"
man=[l.split("\t")[0].strip() for l in open(B+"/results/ppp_unified/manifest.tsv") if l.strip()]
tot=Counter(grp[g] for g in man if g in grp)
print("\ngroups searched:", dict(tot))
def report(path, label, floor=30.0, qmap=None):
    if not os.path.exists(path): return
    hit=defaultdict(set)
    for line in open(path):
        c=line.rstrip("\n").split("\t")
        if float(c[2])<floor: continue
        q=qmap.get(c[0], c[0]) if qmap else label
        g=c[1].split("|",1)[0].split("::")[-1]
        if g in grp: hit[q].add((grp[g],g))
    for q in sorted(hit):
        d=Counter(k for k,_ in hit[q])
        print("  %-14s %s" % (q, {k:"%d/%d"%(len(set(x for kk,x in hit[q] if kk==k)),tot[k]) for k in ("candidate","Akkermansia","free-living")}))
print("\nat a 30 percent identity floor:")
report(W+"/ppp_trio.tsv", None, 30.0, QN)
report(W+"/control_gh20.tsv", "GH20", 30.0)
report(W+"/pgl.tsv", "pgl", 30.0)
report(W+"/sgl.tsv", "SGL", 30.0)
PY
rm -f $W/db.dmnd $W/all.faa
