#!/bin/bash
#SBATCH --job-name=ftplace
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ftplace_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/ftplace_%j.err
module load fasttree/2.1.11
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/novel_akk_tree
ALN=$W/akk_placement.faa
[ -s "$ALN" ] || { echo "REFUSED: missing $ALN"; exit 1; }
echo "sequences: $(grep -c '^>' $ALN)"
FastTreeMP -lg -gamma $ALN > $W/akk_placement_rerun.nwk 2> $W/fasttree_rerun.log
echo "EXIT=$?"
tail -3 $W/fasttree_rerun.log
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 -c "
import re
for p in ('$W/akk_placement.nwk','$W/akk_placement_rerun.nwk'):
    t=open(p).read()
    print(p.split('/')[-1], 'tips:', len(re.findall(r'[(,]([^(),:]+):', t)))"
