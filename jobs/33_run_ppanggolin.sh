#!/bin/bash
# Pangenome partitioning
# Source: ch3-chitin-evolution/scripts/run_ppanggolin.sh
# Output: results/pangenome/ppang_out

#SBATCH --job-name=ppang
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ppang_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ppang_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH -w x01
module load workspace/scratch
module load mmseqs2/15-6f452
PP=/opt/linux/rocky/8.x/x86_64/pkgs/ppanggolin/2.2.3/env/bin/ppanggolin
P=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome
echo "[START] $(date) host $(hostname)"
rm -rf $P/ppang_out_id40
$PP workflow --anno $P/ppanggolin_anno.tsv --fasta $P/ppanggolin_input.tsv \
  --output $P/ppang_out_id40 --cpu 16 --tmpdir $SCRATCH --identity 0.4 --coverage 0.5
echo "[EXIT] $?"
echo "[DONE] $(date)"
ls -la $P/ppang_out 2>/dev/null
