#!/bin/bash
#SBATCH --job-name=r232chk
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/r232chk_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/r232chk_%j.err
module load diamond
module load prodigal 2>/dev/null
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/r232_new_genera
mkdir -p $W/fna $W/faa
for A in GCA_044373535.1 GCA_964576595.1 GCA_044373465.1; do
  P1=$(echo $A | cut -c5-7); P2=$(echo $A | cut -c8-10); P3=$(echo $A | cut -c11-13)
  URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/${P1}/${P2}/${P3}"
  DIR=$(curl -s ${URL}/ | grep -o "${A}[^\"/]*" | head -1)
  if [ -z "$DIR" ]; then echo "REFUSED: could not resolve $A at $URL"; continue; fi
  echo "fetching $A as $DIR"
  curl -sL "${URL}/${DIR}/${DIR}_genomic.fna.gz" -o $W/fna/${A}.fna.gz
  gunzip -f $W/fna/${A}.fna.gz
  echo "  contigs: $(grep -c '>' $W/fna/${A}.fna)"
  prodigal -i $W/fna/${A}.fna -a $W/faa/${A}.faa -p meta -q > /dev/null 2>&1
  echo "  proteins: $(grep -c '>' $W/faa/${A}.faa)"
done
TYPE=$CH3/results/pangenome/gff199/UHM1073.23039_R.bin.125.faa
echo
echo "=== AAI to the proposed type genome, reciprocal best hits ==="
for A in GCA_044373535.1 GCA_964576595.1 GCA_044373465.1; do
  Q=$W/faa/${A}.faa
  [ -s "$Q" ] || { echo "$A  NO PROTEOME"; continue; }
  diamond makedb --in $TYPE -d $W/t --quiet
  diamond makedb --in $Q    -d $W/q --quiet
  diamond blastp -q $Q -d $W/t -o $W/f.tsv --outfmt 6 qseqid sseqid pident --max-target-seqs 1 --more-sensitive --evalue 1e-5 --quiet --threads 8
  diamond blastp -q $TYPE -d $W/q -o $W/r.tsv --outfmt 6 qseqid sseqid pident --max-target-seqs 1 --more-sensitive --evalue 1e-5 --quiet --threads 8
  /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - "$A" << 'PY'
import sys, statistics as st
W="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/r232_new_genera"
def best(p):
    b={}
    for l in open(p):
        f=l.split("\t"); v=float(f[2])
        if f[0] not in b or v>b[f[0]][1]: b[f[0]]=(f[1],v)
    return b
fwd,rev=best(W+"/f.tsv"),best(W+"/r.tsv")
ids=[v for q,(s,v) in fwd.items() if rev.get(s,(None,))[0]==q]
print("%-20s AAI %6.2f   n_RBH %d" % (sys.argv[1], st.mean(ids) if ids else float("nan"), len(ids)))
PY
  rm -f $W/t.dmnd $W/q.dmnd $W/f.tsv $W/r.tsv
done
echo
echo "reference: within-genus AAI to this type genome is 60.38 to 99.97; other Akkermansiaceae genera are 48.42 to 50.87"
