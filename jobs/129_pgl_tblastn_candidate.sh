# TBLASTN of the pgl orthogroup OG0001762 against all 105 candidate assemblies, GH20 control
# Source: ch3-chitin-evolution/scripts/pgl_tblastn_candidate.sh
# Output: results/pgl_tblastn_candidate/pgl_tblastn_hits.tsv
#!/bin/bash
#SBATCH --job-name=pglTB
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pglTB_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pglTB_%j.err
module load ncbi-blast
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/pgl_tblastn_candidate
mkdir -p $W
cp $CH3/results/pgl_search/pgl_query.faa $W/queries.faa
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/pgl_tblastn_candidate"
out=[]; keep=False; taken=False
for line in open(B+"/results/pangenome/control_gh20.faa"):
    if line.startswith(">"):
        keep = not taken; taken = True
        if keep: out.append(">GH20_control\n")
    elif keep: out.append(line)
open(W+"/queries.faa","a").writelines(out)
print("queries:", [l.strip()[1:] for l in open(W+"/queries.faa") if l.startswith(">")])
PY
# the 105 candidate assemblies
N=0
: > $W/cand105.fna
while read -r P; do
  [ -s "$P" ] || { echo "REFUSED: missing $P"; exit 1; }
  G=$(basename "$P" .fa)
  sed "s/^>/>${G}|/" "$P" >> $W/cand105.fna
  N=$((N+1))
done < $CH3/results/akkfam_derep/NOVEL.list
echo "candidate genomes: $N   contigs: $(grep -c '^>' $W/cand105.fna)"
[ "$N" -eq 105 ] || { echo "REFUSED: expected 105 genomes, got $N"; exit 1; }
makeblastdb -in $W/cand105.fna -dbtype nucl -out $W/db105 > /dev/null
tblastn -query $W/queries.faa -db $W/db105 -out $W/pgl_tblastn_hits.tsv \
  -outfmt "6 qseqid sseqid pident length evalue" -evalue 1e-3 -num_threads 16 -max_target_seqs 100000
echo
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
from collections import defaultdict
W="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pgl_tblastn_candidate"
h=defaultdict(set); best=defaultdict(float)
for line in open(W+"/pgl_tblastn_hits.tsv"):
    c=line.split("\t"); h[c[0]].add(c[1].split("|",1)[0]); best[c[0]]=max(best[c[0]], float(c[2]))
for q in ("OG0001762","GH20_control"):
    print("  %-14s %3d of 105 candidate genomes   best identity %.1f" % (q, len(h.get(q,())), best.get(q,0.0)))
PY
rm -f $W/db105.n* $W/cand105.fna
