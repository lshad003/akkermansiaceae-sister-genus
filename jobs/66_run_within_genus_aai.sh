#!/bin/bash
# Within-genus identity calculation submitted
# Source: ch3-chitin-evolution/scripts/run_within_genus_aai.sh
# Output: results/within_genus_aai/within_genus_aai_to_type.tsv

#SBATCH --job-name=aai105
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/aai105_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/aai105_%j.err

module load diamond

/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/within_genus_aai.py
