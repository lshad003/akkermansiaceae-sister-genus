#!/bin/bash
# Uniform gene calls, family reference set
# Source: ch3-chitin-evolution/scripts/run_prodigal_akkfam2.sh
# Output: results/pangenome/akkfam

#SBATCH --job-name=prodfam2
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/prodfam2_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/prodfam2_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=short

module load workspace/scratch
module load prodigal
O=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/akkfam
A=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/list2.tsv | cut -f1)
P=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $O/list2.tsv | cut -f3)
[ -z "$A" ] && exit 0
[ -s "$O/${A}.faa" ] && { echo "done: $A"; exit 0; }
echo "[$(date)] $A"
zcat "$P" > $SCRATCH/${A}.fna
prodigal -i $SCRATCH/${A}.fna -o $O/${A}.gff -a $O/${A}.faa -p meta -f gff -q
echo "[EXIT] $? proteins=$(grep -c '^>' $O/${A}.faa)"
