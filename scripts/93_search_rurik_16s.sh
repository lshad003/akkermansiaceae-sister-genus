#!/bin/bash
# Candidate 16S matched against an independent amplicon survey
# Source: ch3-chitin-evolution/scripts/search_rurik_16s.sh
# Output: results/rurik_16s/hits_genus.tsv
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
OUT=$CH3/results/rurik_16s
DB=$OUT/rurik_otus.fasta

module load ncbi-blast 2>/dev/null || module load blast 2>/dev/null || echo "no blast module, trying PATH"
which makeblastdb blastn

echo ""
echo "=== building db ==="
makeblastdb -in $DB -dbtype nucl -out $OUT/rurikdb -logfile $OUT/makedb.log
echo "done"
echo ""

echo "=== our verified genus 16S against it ==="
blastn -query $CH3/results/rrna16s/verified_genus_16S.fna -db $OUT/rurikdb \
  -outfmt "6 qseqid sseqid pident length qstart qend evalue bitscore" \
  -perc_identity 90 -max_target_seqs 5000 -num_threads 8 -out $OUT/hits_genus.tsv
echo "hits at 90% or above: $(grep -c '' $OUT/hits_genus.tsv)"
echo ""
echo "identity distribution of those hits:"
cut -f3 $OUT/hits_genus.tsv | sort -n | awk '{a[NR]=$1} END{if(NR)print "  min",a[1]," median",a[int(NR/2)]," max",a[NR]}'
echo ""
echo "=== hits at 97% or above, i.e. likely our genus ==="
awk -F'\t' '$3>=97 && $4>=200' $OUT/hits_genus.tsv | sort -k3,3nr | head -20
echo "--- count:"
awk -F'\t' '$3>=97 && $4>=200' $OUT/hits_genus.tsv | wc -l
echo "--- distinct sequences matched at 97%:"
awk -F'\t' '$3>=97 && $4>=200{print $2}' $OUT/hits_genus.tsv | sort -u | wc -l
echo ""
echo "=== control: how many match at 99% or above ==="
awk -F'\t' '$3>=99 && $4>=200{print $2}' $OUT/hits_genus.tsv | sort -u | wc -l
