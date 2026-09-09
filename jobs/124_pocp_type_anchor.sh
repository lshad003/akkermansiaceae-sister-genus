#!/bin/bash
#SBATCH --job-name=pocpT
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pocpT_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pocpT_%j.err
module load diamond/2.2.6
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
S=$CH3/scripts
echo "diamond: $(which diamond)"; diamond --version
mkdir -p $CH3/results/aai_pocp_type
sed -e 's|D="/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond"|D=os.environ.get("DIAMOND_BIN","diamond")|' \
    -e 's|TYPE="EHM058340"|TYPE="UHM1073.23039_R.bin.125"|' \
    -e 's|results/aai_pocp|results/aai_pocp_type|g' \
    "$S/pocp_fixed.py" > "$S/pocp_fixed_type.py"
grep -q "^import os" "$S/pocp_fixed_type.py" || sed -i '1i import os' "$S/pocp_fixed_type.py"
grep -n 'TYPE=' "$S/pocp_fixed_type.py" | head -2
DIAMOND_BIN=$(which diamond) /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 "$S/pocp_fixed_type.py"
