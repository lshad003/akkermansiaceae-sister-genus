#!/bin/bash
module load diamond
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/novinvenio_check
mkdir -p $W
FAA=$CH3/data/UHM_Akkermansia.filtered.faa
[ -s "$FAA" ] || { echo "REFUSED: missing $FAA"; exit 1; }
echo "novelty proteins in file: $(grep -c '>' $FAA)"
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import sys
P="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/ppp_trio.faa"
O="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/novinvenio_check/ppp_queries.faa"
WANT={"UHM979.41089_R.bin.103_CDS_0654":"zwf",
      "UHM1210.23070_R.bin.101_CDS_0293":"gnd",
      "EHM058980_CDS_1432":"opcA"}
keep=None; out=[]
for line in open(P):
    if line.startswith(">"):
        h=line[1:].split()[0]
        keep=WANT.get(h)
        if keep: out.append(">"+keep+"\n")
    elif keep: out.append(line)
open(O,"w").writelines(out)
print("queries written:", sum(1 for l in out if l.startswith(">")))
PY
diamond makedb --in $FAA -d $W/nov --quiet
diamond blastp -q $W/ppp_queries.faa -d $W/nov -o $W/ppp_vs_novelty.tsv \
  --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 --more-sensitive --quiet --threads 4
echo
echo "=== PPP queries vs the 651 novelty candidates ==="
if [ -s $W/ppp_vs_novelty.tsv ]; then column -t $W/ppp_vs_novelty.tsv; else echo "NO HITS"; fi
rm -f $W/nov.dmnd
