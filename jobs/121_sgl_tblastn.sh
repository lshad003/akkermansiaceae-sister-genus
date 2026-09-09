#!/bin/bash
#SBATCH --job-name=sglTB
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/sglTB_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/sglTB_%j.err
module load ncbi-blast
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/sgl_prevalence
G=$CH3/results/ppp_tblastn/genomes.tsv
[ -s "$W/query.faa" ] || { echo "REFUSED: missing $W/query.faa"; exit 1; }
[ -s "$G" ] || { echo "REFUSED: missing $G"; exit 1; }
cat $W/query.faa > $W/tb_queries.faa
# add GH20 as the positive control, same as job 116
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/sgl_prevalence"
out=[]; keep=False; taken=False
for line in open(B+"/results/pangenome/control_gh20.faa"):
    if line.startswith(">"):
        keep = not taken; taken = True
        if keep: out.append(">GH20_control\n")
    elif keep: out.append(line)
open(W+"/tb_queries.faa","a").writelines(out)
print("queries:", [l.strip()[1:] for l in open(W+"/tb_queries.faa") if l.startswith(">")])
PY
: > $W/akk154.fna
while IFS=$'\t' read -r NAME PATHF; do sed "s/^>/>${NAME}|/" "$PATHF" >> $W/akk154.fna; done < "$G"
echo "contigs: $(grep -c '^>' $W/akk154.fna)"
makeblastdb -in $W/akk154.fna -dbtype nucl -out $W/db154 > /dev/null
tblastn -query $W/tb_queries.faa -db $W/db154 -out $W/sgl_tblastn_hits.tsv \
  -outfmt "6 qseqid sseqid pident length evalue" -evalue 1e-3 -num_threads 16 -max_target_seqs 100000
echo
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
from collections import defaultdict
W="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/sgl_prevalence"
tot=set(l.split("\t")[0] for l in open("/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/ppp_tblastn/genomes.tsv") if l.strip())
h=defaultdict(set); best=defaultdict(float)
for line in open(W+"/sgl_tblastn_hits.tsv"):
    c=line.split("\t"); g=c[1].split("|",1)[0]
    h[c[0]].add(g); best[c[0]]=max(best[c[0]], float(c[2]))
for q in ("SGL_query","GH20_control"):
    print("  %-14s %3d of %d genomes   best identity %.1f" % (q, len(h.get(q,())), len(tot), best.get(q,0.0)))
PY
rm -f $W/db154.n* $W/akk154.fna
