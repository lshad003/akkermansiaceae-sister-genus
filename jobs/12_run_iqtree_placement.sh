#!/bin/bash
# Placement re-inferred by maximum likelihood with branch support
# Source: ch3-chitin-evolution/scripts/run_iqtree_placement.sh
# Output: results/novel_akk_tree/iqtree_placement

#SBATCH --job-name=iqplace
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/iqtree_placement_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/iqtree_placement_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --partition=epyc

module load iqtree/2.2.2.6
BIN=/opt/linux/rocky/8.x/x86_64/pkgs/iqtree/2.2.2.6/bin
echo "--- contents of module bin ---"
ls $BIN
IQ=""
for c in $BIN/iqtree2 $BIN/iqtree $BIN/iqtree3; do
  if [ -x "$c" ]; then IQ=$c; break; fi
done
if [ -z "$IQ" ]; then echo "no iqtree binary found in $BIN"; exit 0; fi
echo "using: $IQ"
$IQ --version | head -2

B=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
IN=$B/results/novel_akk_tree/akk_placement.faa
O=$B/results/novel_akk_tree/iqtree_placement
mkdir -p $O
echo "sequences: $(grep -c '^>' $IN)"
CPU=${SLURM_CPUS_ON_NODE:-16}

$IQ -s $IN -m LG+G4 -B 1000 -alrt 1000 -T $CPU -seed 20260723 \
    -pre $O/akk_placement_iqtree
echo "EXIT=$?"
ls -la $O/
