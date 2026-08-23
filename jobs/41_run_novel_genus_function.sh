#!/bin/bash
# Functional panel submitted
# Source: ch3-chitin-evolution/scripts/run_novel_genus_function.sh
# Output: stdout

#SBATCH --job-name=novel_genus_fn
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/novel_genus_fn_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/novel_genus_fn_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=short

mkdir -p /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs
echo "[START] $(date)  host $(hostname)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/novel_genus_function.py
echo "[DONE] $(date)"
