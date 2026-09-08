#!/bin/bash
#SBATCH --job-name=pppTB
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pppTB_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/pppTB_%j.err
module load ncbi-blast
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/ppp_tblastn
mkdir -p $W
# 1. queries: zwf, gnd, opcA, plus GH20 as a positive control
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
import csv, os, sys
B="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=B+"/results/ppp_tblastn"
WANT={"UHM979.41089_R.bin.103_CDS_0654":"zwf",
      "UHM1210.23070_R.bin.101_CDS_0293":"gnd",
      "EHM058980_CDS_1432":"opcA"}
out=[]; keep=None
for line in open(B+"/results/pangenome/ppp_trio.faa"):
    if line.startswith(">"):
        keep=WANT.get(line[1:].split()[0])
        if keep: out.append(">"+keep+"\n")
    elif keep: out.append(line)
# positive control from the GH20 query file
keep=None; n0=len(out)
for line in open(B+"/results/pangenome/control_gh20.faa"):
    if line.startswith(">"):
        keep = (len(out)==n0)          # take only the first sequence
        if keep: out.append(">GH20_control\n")
    elif keep: out.append(line)
open(W+"/queries.faa","w").writelines(out)
print("queries:", [l.strip()[1:] for l in out if l.startswith(">")])

# 2. genome list: Akkermansia in the manifest that have an assembly
man=set(l.split("\t")[0].strip() for l in open(B+"/results/ppp_unified/manifest.tsv") if l.strip())
cen=list(csv.DictReader(open(B+"/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"),delimiter="\t"))
akk=set(r["accession"] for r in cen if r["genus"]=="Akkermansia" and r["accession"] in man)
paths={}
for d in ["/data/amphibia_gtdbtk_input","/results/novel_akk_ani/ref_fna"]:
    p=B+d
    if os.path.isdir(p):
        for f in os.listdir(p):
            for e in (".fa",".fna",".fasta"):
                if f.endswith(e) and f[:-len(e)] in akk:
                    paths.setdefault(f[:-len(e)], os.path.join(p,f))
if len(paths)!=154:
    print("REFUSED: found %d assemblies, expected 154" % len(paths)); sys.exit(1)
with open(W+"/genomes.tsv","w") as fh:
    for g,p in sorted(paths.items()): fh.write(g+"\t"+p+"\n")
print("genomes to search:", len(paths))
PY
[ -s $W/genomes.tsv ] || { echo "REFUSED: no genome list"; exit 1; }
# 3. one combined nucleotide database with genome name prefixed on every contig
: > $W/all_akk.fna
while IFS=$'\t' read -r G P; do
  sed "s/^>/>${G}|/" "$P" >> $W/all_akk.fna
done < $W/genomes.tsv
echo "contigs in database: $(grep -c '^>' $W/all_akk.fna)"
makeblastdb -in $W/all_akk.fna -dbtype nucl -out $W/akkdb > /dev/null
# 4. tblastn, permissive on purpose: we want to find anything that is there
tblastn -query $W/queries.faa -db $W/akkdb -out $W/ppp_tblastn_hits.tsv \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -evalue 1e-3 -num_threads 16 -max_target_seqs 100000
echo
echo "=== hits per query, and how many distinct genomes ==="
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PY'
from collections import defaultdict
W="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/ppp_tblastn"
tot=set(l.split("\t")[0] for l in open(W+"/genomes.tsv") if l.strip())
h=defaultdict(set); best=defaultdict(float)
for line in open(W+"/ppp_tblastn_hits.tsv"):
    c=line.split("\t"); g=c[1].split("|",1)[0]
    h[c[0]].add(g); best[c[0]]=max(best[c[0]], float(c[2]))
for q in ("zwf","gnd","opcA","GH20_control"):
    print("  %-14s %3d of %d genomes   best identity %.1f" % (q, len(h.get(q,())), len(tot), best.get(q,0.0)))
PY
rm -f $W/akkdb.n* $W/all_akk.fna
