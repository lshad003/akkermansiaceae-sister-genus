#!/bin/bash
#SBATCH --job-name=pocp2
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pocp2_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pocp2_%j.err
module load diamond
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
S=$CH3/scripts
D=$(which diamond)
echo "diamond: $D"; $D --version
for F in pocp_fixed.py pocp_calibrate.py pocp_stability.py; do
  [ -s "$S/$F" ] || { echo "note: $S/$F not present, skipped"; continue; }
  N="${F%.py}_d2.py"
  sed 's|D="/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond"|D=os.environ.get("DIAMOND_BIN","diamond")|' "$S/$F" > "$S/$N"
  grep -q "^import os" "$S/$N" || sed -i '1i import os' "$S/$N"
  sed -i 's|results/aai_pocp|results/aai_pocp_d2|g' "$S/$N"
  echo "=== running $N ==="
  DIAMOND_BIN="$D" /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 "$S/$N"
done
