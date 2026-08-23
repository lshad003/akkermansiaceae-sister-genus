#!/bin/bash -l
# Chimerism screening submitted for the candidate and sister genomes
# Source: ch3-chitin-evolution/scripts/run_gunc_199.sh
# Output: results/gunc_199/GUNC.progenomes_2.1.maxCSS_level.tsv

# GUNC chimerism screening of the 199 Prodigal proteomes (105 candidate genus + 94 amphibian
# Akkermansia) against proGenomes 2.1.
# Output: results/gunc_199/
#SBATCH --job-name=gunc199
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/gunc199_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/gunc199_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --partition=batch
#
# RATIONALE: 199 proteomes, roughly a sixth of the 1,171 in the Bacillota run, so 6h is generous.
# 64G and 24 threads carried over unchanged from the run that worked.
#
# Gene calling is SKIPPED. jobs/28_run_prodigal_199_array.sh line 27 used "prodigal -p meta",
# which is GUNC's own default, so --gene_calls is valid here.
#
# DO NOT "module load diamond". GUNC 1.0.6 pins DIAMOND to exactly 2.0.4 and ships that binary.
# Loading the diamond module puts 2.1.24 ahead of it on PATH and GUNC aborts.
#
# Input is read straight from results/pangenome/gff199, which also holds 199 .gff files;
# --file_suffix .faa selects only the proteomes. Nothing is staged, because /bigdata is at 98%.
#
# Log paths are ABSOLUTE. NO module purge. NO set -ue.

WORKDIR=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
IN=${WORKDIR}/results/pangenome/gff199
OUT=${WORKDIR}/results/gunc_199
DB=/srv/projects/db/GUNC/gunc_db_progenomes2.1.dmnd

source /etc/profile.d/modules.sh || true
module load gunc

if ! command -v gunc > /dev/null; then
    echo "[FATAL] gunc not on PATH" >&2
    exit 1
fi
if ! command -v diamond > /dev/null; then
    echo "[FATAL] diamond not on PATH" >&2
    exit 1
fi

echo "[INFO] gunc:    $(command -v gunc)"
echo "[INFO] diamond: $(command -v diamond)"
gunc --version
DV=$(diamond --version 2>&1 | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
echo "[INFO] diamond version: ${DV}"
if [ "${DV}" != "2.0.4" ]; then
    echo "[FATAL] GUNC 1.0.6 requires diamond 2.0.4, found ${DV}." >&2
    echo "        Something put another diamond ahead of GUNC's own on PATH." >&2
    exit 1
fi

if [ ! -s "${DB}" ]; then
    echo "[FATAL] GUNC database missing: ${DB}" >&2
    exit 1
fi
echo "[INFO] db: ${DB}"

N=$(ls -1 ${IN}/*.faa 2>/dev/null | wc -l)
if [ "${N}" -ne 199 ]; then
    echo "[FATAL] expected 199 .faa in ${IN}, found ${N}" >&2
    exit 1
fi
echo "[INFO] proteomes: ${N}"

mkdir -p ${OUT}
export TMPDIR=${WORKDIR}/results/gunc_199/tmp
mkdir -p ${TMPDIR}

echo "[START] $(date) host=$(hostname)"
gunc run \
    --input_dir ${IN} \
    --gene_calls \
    --file_suffix .faa \
    --db_file ${DB} \
    --out_dir ${OUT} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --temp_dir ${TMPDIR} \
    --detailed_output
RC=$?
echo "[DONE] $(date) rc=${RC}"
if [ ${RC} -ne 0 ]; then
    echo "[FATAL] gunc rc=${RC}" >&2
    exit 1
fi

SUM=$(ls -1 ${OUT}/GUNC.*.maxCSS_level.tsv 2>/dev/null | head -1)
if [ -s "${SUM}" ]; then
    echo "[INFO] summary: ${SUM}"
    echo "[INFO] rows: $(( $(wc -l < ${SUM}) - 1 )) of ${N}"
    echo "[INFO] pass.GUNC at default CSS 0.45:"
    awk -F'\t' 'NR>1{print $NF}' ${SUM} | sort | uniq -c
    echo "[INFO] taxonomic level of max CSS:"
    awk -F'\t' 'NR>1{print $2}' ${SUM} | sort | uniq -c | sort -rn
else
    echo "[WARN] no maxCSS_level.tsv in ${OUT}"
    ls -1 ${OUT} | head -10
fi
echo "GUNC_FINISHED"
