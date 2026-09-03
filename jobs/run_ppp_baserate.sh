#!/bin/bash -l
#SBATCH --job-name=ppp_base
#SBATCH -p batch -N 1 -n 1 -c 16 --mem=64G --time=8:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ppp_base_%j.out
#SBATCH -e /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ppp_base_%j.err

R=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$R/results/ppp_baserate
S=$R/results/dbcan_scaffold
PG=$R/results/pangenome
mkdir -p "$W"

D=$(ls /bigdata/stajichlab/lshad003/condaenvs/*/bin/diamond 2>/dev/null | head -1)
echo "diamond: $D"
[ -x "$D" ] || { echo "no diamond"; exit 1; }

echo "proteomes: $(ls $S/*.faa 2>/dev/null | wc -l)"

# one concatenated db, headers prefixed with the genome id
if [ ! -s "$W/all.faa" ]; then
  for f in "$S"/*.faa; do
    B=$(basename "$f" .faa)
    awk -v g="$B" '/^>/{print ">"g"|"substr($0,2);next}{print}' "$f"
  done > "$W/all.faa"
fi
echo "concatenated proteins: $(grep -c '^>' $W/all.faa)"

"$D" makedb --in "$W/all.faa" -d "$W/alldb" --threads 16 --quiet

for Q in ppp_trio control_gh20; do
  [ -s "$PG/$Q.faa" ] || { echo "MISSING $PG/$Q.faa"; continue; }
  echo "searching $Q ..."
  "$D" blastp -q "$PG/$Q.faa" -d "$W/alldb" -o "$W/$Q.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 \
    --threads 16 --max-target-seqs 1000000 --quiet
  echo "  hits: $(wc -l < $W/$Q.tsv)"
done

/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $R/scripts/ppp_baserate_report.py
