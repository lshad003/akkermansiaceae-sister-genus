#!/bin/bash
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
OUT=$CH3/results/rrna16s
module load ncbi-blast 2>/dev/null
echo "=== 16S sequences on file ==="
grep -c ">" $OUT/verified_genus_16S.fna $OUT/candidate_16S_full.fna
grep ">" $OUT/verified_genus_16S.fna
echo
echo "=== full-length 16S identity vs A. muciniphila, genus threshold is about 94.5 ==="
blastn -query $OUT/verified_genus_16S.fna \
  -subject $CH3/results/rurik_16s/akkref/akk_16S.fna \
  -outfmt "6 qseqid pident length qlen slen evalue" -out $OUT/vs_amuc.tsv
column -t $OUT/vs_amuc.tsv
