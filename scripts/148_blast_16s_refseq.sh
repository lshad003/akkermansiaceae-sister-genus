#!/bin/bash
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
Q=$CH3/results/rrna16s/verified_genus_16S.fna
DB=$CH3/db/16S_refseq/16S_ribosomal_RNA
OUT=$CH3/results/rrna16s/vs_16S_refseq.tsv

[ -s "$Q" ] || { echo "REFUSED: missing $Q"; exit 1; }
echo "queries: $(grep -c '^>' $Q)"

module load ncbi-blast 2>/dev/null || module load blast 2>/dev/null
blastn -version

export BLASTDB=$CH3/db/16S_refseq
blastn -query $Q -db $DB \
  -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore staxid ssciname stitle" \
  -max_target_seqs 20 -num_threads 4 -out $OUT

echo "hit lines: $(wc -l < $OUT)"
echo
echo "=== best hit per query ==="
sort -k1,1 -k3,3gr $OUT | awk -F'\t' '!seen[$1]++ {printf "%-32s %6.2f%%  %5d bp  %s\n", $1, $3, $4, $10}'
echo
echo "=== top 5 overall by identity ==="
sort -k3,3gr $OUT | head -5 | cut -f1,3,4,10,11 | cut -c1-150
