#!/bin/bash
# Pathway search rerun against a genus-defined Akkermansia set with a positive control
# Source: ch3-chitin-evolution/scripts/rebuild_ppp_akkermansia.sh
# Output: results/pangenome/ppp_redo/akkermansia_list.txt

#SBATCH --job-name=pppredo
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pppredo_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pppredo_%j.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -p short

CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
PG=$CH3/results/pangenome
OUT=$PG/ppp_redo
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
mkdir -p $OUT

module load diamond

echo "[START] $(date)"
$PY << 'PYEOF'
import os, csv
CH3="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
PG=CH3+"results/pangenome/"
OUT=PG+"ppp_redo/"
cen=CH3+"results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
gen={}
for r in csv.DictReader(open(cen),delimiter="\t"):
    gen[r["accession"]]=(r["genus"] or "").strip()
have={}
for d in ("akkfam","gtdb_akk","gff199"):
    p=PG+d
    if os.path.isdir(p):
        for x in os.listdir(p):
            if x.endswith(".faa"):
                have.setdefault(x[:-4],p+"/"+x)
akk=[a for a in have if gen.get(a)=="Akkermansia"]
print("proteomes available:",len(have))
print("GTDB genus Akkermansia among them:",len(akk))
import collections
print("genus of the rest, top 8:",dict(collections.Counter(
    gen.get(a,"NOT_IN_CENSUS") for a in have if gen.get(a)!="Akkermansia").most_common(8)))
with open(OUT+"akkermansia_only.faa","w") as o:
    for a in sorted(akk):
        for line in open(have[a]):
            o.write((">"+a+"|"+line[1:]) if line.startswith(">") else line)
open(OUT+"akkermansia_list.txt","w").write("\n".join(sorted(akk))+"\n")
print("wrote akkermansia_only.faa")
PYEOF

echo ""
echo "proteins pooled: $(grep -c '>' $OUT/akkermansia_only.faa)"
diamond makedb --in $OUT/akkermansia_only.faa -d $OUT/akkdb --threads 16

echo ""
echo "=== POSITIVE CONTROL, GH20 ==="
diamond blastp -q $PG/control_gh20.faa -d $OUT/akkdb -o $OUT/control.tsv \
  --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 --threads 16 --max-target-seqs 5000
echo "control hits: $(grep -c '' $OUT/control.tsv)"
echo "genomes hit: $(cut -f2 $OUT/control.tsv | cut -d'|' -f1 | sort -u | wc -l) of $(grep -c '' $OUT/akkermansia_list.txt)"

echo ""
echo "=== PPP, strict ==="
diamond blastp -q $PG/ppp_trio.faa -d $OUT/akkdb -o $OUT/ppp.tsv \
  --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 --threads 16 --max-target-seqs 5000
echo "PPP hits: $(grep -c '' $OUT/ppp.tsv)"
cut -f2 $OUT/ppp.tsv 2>/dev/null | cut -d'|' -f1 | sort -u | wc -l
head -20 $OUT/ppp.tsv

echo ""
echo "=== PPP, very sensitive, evalue 1 ==="
diamond blastp -q $PG/ppp_trio.faa -d $OUT/akkdb -o $OUT/ppp_relaxed.tsv \
  --outfmt 6 qseqid sseqid pident length evalue --evalue 1 --threads 16 --max-target-seqs 5000 --very-sensitive
echo "relaxed hits: $(grep -c '' $OUT/ppp_relaxed.tsv)"
head -20 $OUT/ppp_relaxed.tsv
echo ""
echo "[DONE] $(date)"
