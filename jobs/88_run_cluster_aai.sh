#!/bin/bash
# Cluster-to-cluster identity searches submitted
# Source: ch3-chitin-evolution/scripts/run_cluster_aai.sh
# Output: results/cluster_aai/cluster_aai_matrix.tsv

#SBATCH --job-name=claai
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/claai_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/claai_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -p short

CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
FAA=$CH3/results/pangenome/gff199
OUT=$CH3/results/cluster_aai
mkdir -p $OUT/db $OUT/hits

module load diamond

REPS="C286:EHM059405 C287:EHM047778 C288:UHM1073.23039_R.bin.125 C289:EHM034774 C290:UHM979.23045_R.bin.163 C291:UHM1171.23069_R.bin.59 C292:UHM967.23060_R.bin.25 C293:UHM989.23064_R.bin.100 C294:UHM1210.41100_R.bin.142 C295:UHM1327.41103_R.bin.200 C296:UHM1327.41103_R.bin.48 C297:EHM059333 C298:STP248.12601_R.bin.118 C299:UHM1088.23067_R.bin.54 C300:UHM298.41072_R.bin.90 C301:UHM905.23056_R.bin.2 C302:UHM973.23044_R.bin.53"

echo "=== rebuilding any missing databases ==="
for R in $REPS; do
  C=${R%%:*}; G=${R##*:}
  if [ ! -s "$OUT/db/$C.dmnd" ]; then
    echo "building $C from $G"
    diamond makedb --in $FAA/$G.faa -d $OUT/db/$C --threads 8
  fi
done
echo ""
echo "databases on disk:"
ls $OUT/db/*.dmnd | wc -l
echo "expected 17"
echo ""

echo "=== filling in missing hit files ==="
N=0
for A in $REPS; do
  CA=${A%%:*}; GA=${A##*:}
  for B in $REPS; do
    CB=${B%%:*}
    if [ "$CA" = "$CB" ]; then continue; fi
    F=$OUT/hits/${CA}_vs_${CB}.tsv
    if [ ! -s "$F" ]; then
      diamond blastp -q $FAA/$GA.faa -d $OUT/db/$CB -o $F \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --max-target-seqs 1 --evalue 1e-5 --quiet --threads 8
      N=$((N+1))
    fi
  done
done
echo "searches run this pass: $N"
echo ""
ls $OUT/hits | wc -l
echo "expected 272"
