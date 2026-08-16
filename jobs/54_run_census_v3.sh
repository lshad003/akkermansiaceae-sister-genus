#!/bin/bash
# Census submitted
# Source: ch3-chitin-evolution/scripts/run_census_v3.sh
# Output: stdout

#SBATCH --job-name=census_v3
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/census_v3_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/census_v3_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --partition=short

echo "[START] $(date) host $(hostname)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
  /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/gh75_verru_census_v3.py
echo "[DONE] $(date)"
