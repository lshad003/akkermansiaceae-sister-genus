#!/bin/bash
# Ancestral losses annotated functionally
# Source: ch3-chitin-evolution/scripts/run_eggnog_polarity_B.sh
# Output: results/pangenome/eggnog_polarity/B_lost_at_gut_ancestor.emapper.annotations

# Functional annotation of the orthogroups lost at the shared gut ancestor
# Output: results/pangenome/eggnog_polarity
#SBATCH --job-name=eggpolB
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/eggnog_polarityB_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/eggnog_polarityB_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
#
# Category B holds about 1,047 orthogroups against 57 for category A, so this is
# roughly eighteen times the work of jobs/39_run_eggnog_polarity.sh. Same emapper
# invocation, same database, same tax scope, so the annotations stay comparable.

module load workspace/scratch
EM=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin/emapper.py
export EGGNOG_DATA_DIR=/srv/projects/db/eggNOG/LATEST
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/eggnog-mapper/2.1.9/env/bin:$PATH

P=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome
R=$P/polarity_reps
O=$P/eggnog_polarity
C=B_lost_at_gut_ancestor

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
