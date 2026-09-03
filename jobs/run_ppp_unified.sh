#!/bin/bash -l
#SBATCH --job-name=ppp_unified
#SBATCH -p batch -N 1 -n 1 -c 16 --mem=64G --time=8:00:00
#SBATCH -o /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ppp_unified_%j.out
#SBATCH -e /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ppp_unified_%j.err

# One search, one database: the 187 Akkermansia, the 105 candidates, and the
# free-living Akkermansiaceae. Removes every cross-run comparison from the result.
R=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$R/results/ppp_unified
PG=$R/results/pangenome
mkdir -p "$W"

D=$(ls /bigdata/stajichlab/lshad003/condaenvs/*/bin/diamond 2>/dev/null | head -1)
echo "diamond: $D"; [ -x "$D" ] || { echo "no diamond"; exit 1; }

# gather proteomes from every dir that holds them, dedup by basename
if [ ! -s "$W/all.faa" ]; then
  : > "$W/manifest.tsv"
  for DIR in $PG/gff199 $PG/akkfam $PG/outgroups $PG/gtdb_akk \
             $R/results/dbcan_scaffold $R/results/dbcan_ehi_amphibian \
             $R/results/dbcan_ehi_nonamph; do
    [ -d "$DIR" ] || continue
    for f in "$DIR"/*.faa; do
      [ -e "$f" ] || continue
      B=$(basename "$f" .faa)
      grep -q -m1 -P "^\Q$B\E\t" "$W/manifest.tsv" 2>/dev/null && continue
      printf '%s\t%s\n' "$B" "$DIR" >> "$W/manifest.tsv"
    done
  done
  echo "distinct genomes: $(wc -l < $W/manifest.tsv)"
  while IFS=$'\t' read -r B DIR; do
    awk -v g="$B" '/^>/{print ">"g"|"substr($0,2);next}{print}' "$DIR/$B.faa"
  done < "$W/manifest.tsv" > "$W/all.faa"
fi
echo "proteins: $(grep -c '^>' $W/all.faa)"

"$D" makedb --in "$W/all.faa" -d "$W/db" --threads 16 --quiet
for Q in ppp_trio control_gh20; do
  [ -s "$PG/$Q.faa" ] || { echo "MISSING $PG/$Q.faa"; continue; }
  "$D" blastp -q "$PG/$Q.faa" -d "$W/db" -o "$W/$Q.tsv" \
    --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 \
    --threads 16 --max-target-seqs 1000000 --quiet
  echo "$Q hits: $(wc -l < $W/$Q.tsv)"
done

/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 $R/scripts/ppp_unified_report.py
