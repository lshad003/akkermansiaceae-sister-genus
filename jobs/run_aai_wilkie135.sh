#!/bin/bash -l
#SBATCH --job-name=wilkie_aai
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/wilkie_aai_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/wilkie_aai_%j.err

R=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3

echo "host    : $(hostname)"
echo "started : $(date)"
echo "job     : $SLURM_JOB_ID"
echo "cpus    : $SLURM_CPUS_PER_TASK"
echo

$PY $R/scripts/aai_wilkie135_run.py
RC=$?

echo
echo "exit    : $RC"
echo "finished: $(date)"
echo
if [ -s "$R/results/wilkie_aai/aai_vs_type.tsv" ]; then
  echo "=== top 20 by AAI ==="
  head -21 "$R/results/wilkie_aai/aai_vs_type.tsv"
  echo
  echo "rows: $(($(wc -l < "$R/results/wilkie_aai/aai_vs_type.tsv") - 1))"
else
  echo "NO OUTPUT TABLE. see the .err file."
fi
