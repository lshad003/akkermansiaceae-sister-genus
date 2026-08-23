#!/bin/bash
# Retained orthogroups annotated as the comparison background
# Source: ch3-chitin-evolution/scripts/run_eggnog_polarity_shared.sh
# Output: results/pangenome/eggnog_polarity/shared_gut.emapper.annotations

# Functional annotation of the orthogroups retained across both gut genera.
# This is the BACKGROUND set. Without it, category proportions describe the losses but
# cannot establish enrichment.
# Output: results/pangenome/eggnog_polarity
#SBATCH --job-name=eggpolS
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/eggnog_polarityS_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/eggnog_polarityS_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#
# Same emapper invocation, database and tax scope as jobs/39 and the category B run,
# so all five categories stay comparable. Category B (1,047) finished well inside
# this window at 16 CPU.

module load workspace/scratch
EM=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin/emapper.py
export EGGNOG_DATA_DIR=/srv/projects/db/eggNOG/LATEST
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin:$PATH

P=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome
R=$P/polarity_reps
O=$P/eggnog_polarity
C=shared_gut

if [ ! -s "$R/${C}.faa" ]; then
    echo "[FATAL] missing $R/${C}.faa" >&2
    exit 1
fi
echo "[INFO] sequences: $(grep -c '>' $R/${C}.faa)"
mkdir -p $O
echo "[START] $(date) host $(hostname)"
$EM --version
$EM -i $R/${C}.faa -o $C \
  --output_dir $O --cpu 16 -m diamond \
  --data_dir $EGGNOG_DATA_DIR \
  --tax_scope prokaryota_broad --temp_dir $SCRATCH
echo "[EXIT] $?"
echo "[DONE] $(date)"
ls -la $O/
