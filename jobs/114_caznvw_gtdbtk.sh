#!/bin/bash
#SBATCH --job-name=caznvwT
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/caznvwT_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/caznvwT_%j.err
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/caznvw_tree
mkdir -p $W/genome
cp $CH3/results/r232_new_genera/fna/GCA_964576595.1.fna $W/genome/CAZNVW01.fna
module load gtdbtk
gtdbtk classify_wf --genome_dir $W/genome --out_dir $W/gtdbtk \
  --extension fna --cpus 16 --skip_ani_screen
echo "=== classification ==="
cat $W/gtdbtk/classify/gtdbtk.bac120.summary.tsv 2>/dev/null | cut -f1,2,14,19,20
