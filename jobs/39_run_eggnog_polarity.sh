#!/bin/bash
# Functional annotation of the polarized orthogroups
# Source: ch3-chitin-evolution/scripts/run_eggnog_polarity.sh
# Output: results/pangenome/eggnog_polarity

#SBATCH --job-name=eggpol
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/eggnog_polarity_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/eggnog_polarity_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH -w x01
module load workspace/scratch
EM=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin/emapper.py
export EGGNOG_DATA_DIR=/srv/projects/db/eggNOG/LATEST
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin:$PATH
P=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome
R=$P/polarity_reps
O=$P/eggnog_polarity
mkdir -p $O
echo "[START] $(date) host $(hostname)"
$EM --version
for C in A_lost_in_Akkermansia C_Akkermansia_enriched D_novel_specific
do
  echo "=== $C ==="
  $EM -i $R/${C}.faa -o $C \
    --output_dir $O --cpu 8 -m diamond \
    --data_dir $EGGNOG_DATA_DIR \
    --tax_scope prokaryota_broad --temp_dir $SCRATCH
  echo "[EXIT $C] $?"
done
echo "[DONE] $(date)"
ls -la $O/
