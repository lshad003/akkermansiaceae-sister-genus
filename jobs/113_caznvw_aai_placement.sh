#!/bin/bash
#SBATCH --job-name=caznvw
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/caznvw_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/caznvw_%j.err
module load diamond
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/r232_new_genera
Q=$W/faa/GCA_964576595.1.faa
IN=$CH3/results/pangenome/orthofinder_in
echo "AAI of CAZNVW01 to each of the 17 candidate species-cluster representatives"
echo "and to Akkermansia, for reference"
for F in $IN/NOVEL_sp*.faa $IN/AKK_clade01*.faa; do
  B=$(basename $F .faa)
  diamond makedb --in $F -d $W/s --quiet
  diamond makedb --in $Q -d $W/q --quiet
  diamond blastp -q $Q -d $W/s -o $W/f.tsv --outfmt 6 qseqid sseqid pident --max-target-seqs 1 --more-sensitive --evalue 1e-5 --quiet --threads 8
  diamond blastp -q $F -d $W/q -o $W/r.tsv --outfmt 6 qseqid sseqid pident --max-target-seqs 1 --more-sensitive --evalue 1e-5 --quiet --threads 8
  /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - "$B" << 'PY'
import sys, statistics as st
W="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/r232_new_genera"
def best(p):
    b={}
    for l in open(p):
        f=l.split("\t"); v=float(f[2])
        if f[0] not in b or v>b[f[0]][1]: b[f[0]]=(f[1],v)
    return b
fwd,rev=best(W+"/f.tsv"),best(W+"/r.tsv")
ids=[v for q,(s,v) in fwd.items() if rev.get(s,(None,))[0]==q]
print("  %-34s AAI %6.2f  n_RBH %d" % (sys.argv[1], st.mean(ids) if ids else float("nan"), len(ids)))
PY
  rm -f $W/s.dmnd $W/q.dmnd $W/f.tsv $W/r.tsv
done
