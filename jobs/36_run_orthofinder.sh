#!/usr/bin/bash -l
# Orthogroup inference
# Source: ch3-chitin-evolution/scripts/run_orthofinder.sh
# Output: results/pangenome/orthofinder_out/run1

#SBATCH -p intel -N 1 -c 32 --mem 200gb --time 48:00:00 -J orthof --out /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/orthofinder.log

IN=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/orthofinder_in
OUT=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/orthofinder_out
mkdir -p $OUT
mkdir -p /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs

module load orthofinder/2.5.5
echo "orthofinder: $(which orthofinder)"
orthofinder -h 2>&1 | head -3
echo "diamond: $(which diamond)"
echo "mcl: $(which mcl)"

echo "input proteomes: $(ls $IN/*.faa | wc -l)"
CPU=${SLURM_CPUS_ON_NODE:-16}
echo "threads: $CPU"

orthofinder -f $IN -o $OUT/run1 -t $CPU -a 8 -S diamond -n novelakk
echo "EXIT=$?"
