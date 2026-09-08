#!/bin/bash
#SBATCH --job-name=pglU
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pglU_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pglU_%j.err
module load diamond
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/ppp_unified
U=$CH3/results/pgl_unified
mkdir -p $U
M=$W/manifest.tsv
N=$(wc -l < $M)
if [ "$N" -ne 7572 ]; then echo "REFUSED: manifest has $N rows, expected 7572"; exit 1; fi
echo "rebuilding the unified database from $N manifest rows"
: > $U/all.faa
while IFS=$'\t' read -r B DIR; do
  F="$DIR/$B.faa"
  if [ -s "$F" ]; then
    sed "s/^>/>${B}|/" "$F" >> $U/all.faa
  else
    echo "MISSING $F" >> $U/missing.txt
  fi
done < "$M"
echo "proteins: $(grep -c '^>' $U/all.faa)"
echo "missing proteomes: $(wc -l < $U/missing.txt 2>/dev/null || echo 0)"
diamond makedb --in $U/all.faa -d $U/db --threads 16 --quiet
for Q in pgl zwf; do
  diamond blastp -q $CH3/results/pgl_search/${Q}_query.faa -d $U/db -o $U/${Q}_hits.tsv \
    --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 \
    --max-target-seqs 100000 --more-sensitive --threads 16 --quiet
  echo "$Q hits: $(wc -l < $U/${Q}_hits.tsv)"
done
rm -f $U/db.dmnd $U/all.faa
