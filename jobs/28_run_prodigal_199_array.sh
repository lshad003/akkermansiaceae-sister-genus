#!/bin/bash
# Uniform gene calls, candidate and sister set
# Source: ch3-chitin-evolution/scripts/run_prodigal_199_array.sh
# Output: results/pangenome/gff199

#SBATCH --job-name=prod199a
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/prod199a_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/prod199a_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --array=1-199

module load prodigal
FA=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/amphibia_gtdbtk_input
P=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome
OUT=$P/gff199
mkdir -p $OUT

G=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $P/ppanggolin_input.tsv | cut -f1)
if [ -z "$G" ]; then echo "no genome at line ${SLURM_ARRAY_TASK_ID}"; exit 0; fi
if [ -s "$OUT/${G}.gff" ]; then echo "already done: $G"; exit 0; fi

echo "[$(date)] $G"
prodigal -i $FA/${G}.fa -o $OUT/${G}.gff -a $OUT/${G}.faa -p meta -f gff -q
echo "[EXIT] $?"
