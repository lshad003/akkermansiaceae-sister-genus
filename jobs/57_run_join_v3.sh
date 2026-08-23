#!/bin/bash
# Family join submitted
# Source: ch3-chitin-evolution/scripts/run_join_v3.sh
# Output: stdout

#SBATCH --job-name=join_v3
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/join_v3_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/join_v3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --partition=short
echo "[START] $(date) host $(hostname)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/join_herptile_verru_family_v3.py
echo "[DONE] $(date)"
