#!/bin/bash
# Placement tree inferred
# Source: ch3-chitin-evolution/scripts/rebuild_placement_tree.sh
# Output: results/novel_akk_tree/akk_placement.nwk
BASE=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
FAA=$BASE/results/novel_akk_tree/akk_placement.faa
OUT=$BASE/results/novel_akk_tree/akk_placement_rebuilt.nwk
FT=$(command -v FastTreeMP || command -v fasttreeMP || command -v FastTree || command -v fasttree)
echo "[START] $(date)"
echo "alignment: $FAA  seqs=$(grep -c '^>' $FAA)"
echo "FastTree binary: $FT"
if [ -z "$FT" ]; then echo "NO FastTree on PATH -- tell Claude, will load module"; exit 1; fi
$FT -lg -gamma $FAA > $OUT
echo "EXIT=$?  tree tips=$(grep -o ',' $OUT | wc -l)"
echo "wrote $OUT"
echo "[DONE] $(date)"
