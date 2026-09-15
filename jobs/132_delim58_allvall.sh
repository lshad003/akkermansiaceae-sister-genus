#!/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64G -t 12:00:00
#SBATCH -J delim58 -o /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/delim58_%j.out

CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/delim58
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
mkdir -p $W/faa $W/hits

module unload diamond 2>/dev/null
module load diamond/2.2.6
diamond --version

$PY - << 'PYEOF'
import os, csv
CH3="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W=CH3+"/results/delim58"
DIRS=[CH3+"/results/pangenome/"+d for d in ("gff199","akkfam","outgroups","gtdb_akk")]
DIRS+=[CH3+"/results/"+d for d in ("dbcan_scaffold","dbcan_ehi_amphibian","dbcan_ehi_nonamph")]
idx={}
for d in DIRS:
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d,f))

sel=[]
for r in csv.DictReader(open(CH3+"/results/cluster_aai/cluster_representatives.tsv"),delimiter="\t"):
    sel.append(("CAND_"+r["cluster_id"], r["representative"]))
for r in csv.DictReader(open(W+"/orthofinder_map.tsv"),delimiter="\t"):
    if r["group"]=="AKK":
        sel.append((r["label"].split("_n")[0], r["genome"]))
for r in csv.DictReader(open(CH3+"/results/aai_pocp_type/pocp_fixed.tsv"),delimiter="\t"):
    if r["genus"]!="Akkermansia":
        sel.append((r["genus"], r["rep"]))

labs=[x[0] for x in sel]
assert len(labs)==len(set(labs)), "REFUSED: duplicate labels"
print("selected: %d  (CAND %d, AKK %d, FREE %d)" % (len(sel),
      sum(1 for l,_ in sel if l.startswith("CAND_")),
      sum(1 for l,_ in sel if l.startswith("AKK_clade")),
      sum(1 for l,_ in sel if not l.startswith(("CAND_","AKK_clade")))))
if len(sel)!=58:
    print("REFUSED: expected 58, got", len(sel)); raise SystemExit(1)

miss=[]
with open(W+"/members.tsv","w") as o:
    o.write("label\tgenome\tpath\n")
    for lab,g in sel:
        p=idx.get(g) or idx.get(g.replace("GB_","").replace("RS_",""))
        if not p: miss.append((lab,g)); continue
        o.write("%s\t%s\t%s\n"%(lab,g,p))
        os.system("cp -f '%s' '%s/faa/%s.faa'"%(p,W,lab))
if miss:
    print("REFUSED, missing:",miss); raise SystemExit(1)
print("staged", len(sel))
PYEOF
[ $? -eq 0 ] || { echo "STAGING FAILED"; exit 1; }

cd $W/faa
for A in *.faa; do diamond makedb --in $A -d ${A%.faa} --threads 16 --quiet; done
for A in *.faa; do
  for B in *.faa; do
    [ "$A" = "$B" ] && continue
    O=$W/hits/${A%.faa}__${B%.faa}.tsv
    [ -s "$O" ] && continue
    diamond blastp -q $A -d ${B%.faa} -o $O \
      --outfmt 6 qseqid sseqid pident length qlen evalue bitscore \
      --evalue 1e-5 --max-target-seqs 100000 --more-sensitive --threads 16 --quiet
  done
done
rm -f $W/faa/*.dmnd
echo "PAIRS: $(ls $W/hits | wc -l)   expected 3306"
