#!/bin/bash
#SBATCH --job-name=val651
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/val651_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/val651_%j.err
module load diamond
module load ncbi-blast
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/validate_651
mkdir -p $W
Q=$CH3/data/UHM_Akkermansia.filtered.faa
[ -s "$Q" ] || Q=$CH3/data/UHM_Akkermansia_filtered.faa
[ -s "$Q" ] || { echo "REFUSED: novelty FASTA not found in $CH3/data"; ls $CH3/data | grep -i filtered; exit 1; }
echo "query proteins: $(grep -c '>' $Q)"

# STEP 1: protein search against all 331 Akkermansia proteomes
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import csv, os, sys
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/validate_651"
man={}
for l in open(B+"/results/ppp_unified/manifest.tsv"):
    if l.strip():
        g,d=l.rstrip("\n").split("\t"); man[g]=d
cen=list(csv.DictReader(open(B+"/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"),delimiter="\t"))
akk=[r["accession"] for r in cen if r["genus"]=="Akkermansia" and r["accession"] in man]
if len(akk)!=331:
    print("REFUSED: %d Akkermansia in manifest, expected 331" % len(akk)); sys.exit(1)
n=0
with open(W+"/akk331.faa","w") as o:
    for g in akk:
        p=os.path.join(man[g], g+".faa")
        if not os.path.exists(p):
            print("REFUSED: missing proteome %s" % p); sys.exit(1)
        for line in open(p):
            o.write(">"+g+"|"+line[1:].split()[0]+"\n" if line.startswith(">") else line)
        n+=1
print("proteomes concatenated: %d" % n)
PY
[ -s $W/akk331.faa ] || { echo "REFUSED: no Akkermansia proteome file"; exit 1; }
echo "Akkermansia proteins: $(grep -c '>' $W/akk331.faa)"
diamond makedb --in $W/akk331.faa -d $W/akk331 --quiet
diamond blastp -q $Q -d $W/akk331 -o $W/step1_protein_hits.tsv \
  --outfmt 6 qseqid sseqid pident length evalue --evalue 1e-5 --more-sensitive \
  --max-target-seqs 5000 --threads 16 --quiet
echo "step 1 protein hits: $(wc -l < $W/step1_protein_hits.tsv)"

# STEP 2: TBLASTN of the survivors against the 154 assemblies
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - "$Q" << 'PY'
import sys
W="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/validate_651"
hit=set(l.split("\t")[0] for l in open(W+"/step1_protein_hits.tsv"))
keep=None; n=0; k=0
with open(W+"/survivors.faa","w") as o:
    for line in open(sys.argv[1]):
        if line.startswith(">"):
            h=line[1:].split()[0]; n+=1
            keep = h not in hit
            if keep: k+=1
        if keep: o.write(line)
print("step 1: %d of %d had NO protein hit in 331 Akkermansia" % (k,n))
PY
G=$CH3/results/ppp_tblastn/genomes.tsv
[ -s "$G" ] || { echo "REFUSED: missing $G, run job 116 first"; exit 1; }
: > $W/akk154.fna
while IFS=$'\t' read -r NAME PATHF; do sed "s/^>/>${NAME}|/" "$PATHF" >> $W/akk154.fna; done < $G
makeblastdb -in $W/akk154.fna -dbtype nucl -out $W/akkdb154 > /dev/null
tblastn -query $W/survivors.faa -db $W/akkdb154 -out $W/step2_tblastn_hits.tsv \
  -outfmt "6 qseqid sseqid pident length evalue" -evalue 1e-3 -num_threads 16 -max_target_seqs 100000
echo "step 2 tblastn hits: $(wc -l < $W/step2_tblastn_hits.tsv)"

# STEP 3: report
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import csv
from collections import Counter
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/validate_651"
surv=[l[1:].split()[0] for l in open(W+"/survivors.faa") if l.startswith(">")]
tb=set(l.split("\t")[0] for l in open(W+"/step2_tblastn_hits.tsv"))
final=[s for s in surv if s not in tb]
print()
print("SURVIVORS after protein search vs 331: %d" % len(surv))
print("of those, with a TBLASTN hit in 154 assemblies: %d" % len(set(surv)&tb))
print("FINAL candidate-specific set: %d" % len(final))
ann={}
for r in csv.DictReader(open(B+"/data/UHM_Akkermansia_novelties.tsv"),delimiter="\t"):
    ann[r["Protein ID"].split(":",1)[1]] = (r.get("Pfam domains") or "", r.get("Product") or "")
with open(W+"/final_candidate_specific.tsv","w") as o:
    o.write("protein\tpfam\tproduct\n")
    for f in final:
        p,pr = ann.get(f,("",""))
        o.write("%s\t%s\t%s\n" % (f,p,pr))
pf=Counter()
for f in final:
    for d in ann.get(f,("",""))[0].split(","):
        if d.strip(): pf[d.strip()]+=1
print("annotated: %d of %d" % (sum(1 for f in final if ann.get(f,("",""))[0]), len(final)))
print("top Pfam:", pf.most_common(15))
print("wrote", W+"/final_candidate_specific.tsv")
PY
rm -f $W/akk331.dmnd $W/akkdb154.n* $W/akk154.fna $W/akk331.faa
