#!/bin/bash
# Identity and POCP calculation submitted
# Source: ch3-chitin-evolution/scripts/run_aai_pocp.sh
# Output: stdout

#SBATCH --job-name=aaipocp
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/aaipocp_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/aaipocp_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=short
echo "[START] $(date) host $(hostname)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/aai_pocp.py
echo "[EXIT] $?"
echo "[DONE] $(date)"
