#!/bin/bash
#SBATCH --job-name=pinver
#SBATCH --partition=short
#SBATCH --time=0:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pinver_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pinver_%j.err
echo "=== default module ==="
module load diamond
which diamond; diamond --version
echo
echo "=== versions available ==="
module avail diamond 2>&1 | tail -3
echo
echo "=== hardcoded path used by the original POCP scripts ==="
/bigdata/stajichlab/lshad003/condaenvs/diamond/bin/diamond --version
echo
echo "=== blast for tblastn ==="
module load ncbi-blast
which tblastn; tblastn -version | head -1
