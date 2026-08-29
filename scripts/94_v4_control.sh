#!/bin/bash
# Marker resolution between the candidate genus and Akkermansia confirmed
# Source: ch3-chitin-evolution/scripts/v4_control.sh
# Output: results/rurik_16s/hits_akk.tsv
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
OUT=$CH3/results/rurik_16s
mkdir -p $OUT/akkref
cd $OUT/akkref

module load barrnap 2>/dev/null
module load ncbi-blast 2>/dev/null

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/020/225/GCF_000020225.1_ASM2022v1/GCF_000020225.1_ASM2022v1_genomic.fna.gz"
echo "=== downloading A. muciniphila ATCC BAA-835 ==="
curl -sL -o amuc.fna.gz "$URL"
ls -l amuc.fna.gz
gunzip -f amuc.fna.gz
echo "size: $(grep -v '>' amuc.fna | tr -d '\n' | wc -c) bp"
echo ""

echo "=== 16S from it ==="
barrnap --quiet amuc.fna 2>/dev/null | grep 16S | head -3
/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 - << 'PYEOF'
import subprocess, os
d="/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/rurik_16s/akkref/"
def readfa(p):
    r={};n=None;b=[]
    for l in open(p):
        if l.startswith(">"):
            if n: r[n]="".join(b)
            n=l[1:].split()[0];b=[]
        else: b.append(l.strip())
    if n: r[n]="".join(b)
    return r
def rc(s): return s[::-1].translate(str.maketrans("ACGTacgt","TGCAtgca"))
seqs=readfa(d+"amuc.fna")
out=subprocess.run(["barrnap","--quiet",d+"amuc.fna"],capture_output=True,text=True).stdout
best=""
for line in out.splitlines():
    if "\t" not in line or "16S" not in line: continue
    f=line.split("\t")
    if f[0] not in seqs: continue
    s=seqs[f[0]][int(f[3])-1:int(f[4])]
    if f[6]=="-": s=rc(s)
    if len(s)>len(best): best=s
open(d+"akk_16S.fna","w").write(">Akkermansia_muciniphila_ATCC_BAA835 len=%d\n%s\n"%(len(best),best))
print("Akkermansia 16S extracted: %d bp"%len(best))
PYEOF
echo ""

echo "=== Akkermansia 16S against the same 150,028 sequences ==="
blastn -query $OUT/akkref/akk_16S.fna -db $OUT/rurikdb \
  -outfmt "6 qseqid sseqid pident length" -perc_identity 90 \
  -max_target_seqs 5000 -num_threads 8 -out $OUT/hits_akk.tsv
awk -F'\t' '$4>=200{if($3>=99)a++; else if($3>=97)b++; else if($3>=95)c++; else d++} END{print "  99-100%:",a+0; print "  97-99% :",b+0; print "  95-97% :",c+0; print "  90-95% :",d+0}' $OUT/hits_akk.tsv
echo ""

echo "=== THE DECIDING TEST: do our 29 also match Akkermansia at 97%+ ==="
awk -F'\t' '$3>=97 && $4>=200{print $2}' $OUT/hits_akk.tsv | sort -u > $OUT/a97.txt
sort -u $OUT/matched_ids.txt > $OUT/g97.txt
echo "  ours only:        $(comm -23 $OUT/g97.txt $OUT/a97.txt | wc -l)"
echo "  Akkermansia only: $(comm -13 $OUT/g97.txt $OUT/a97.txt | wc -l)"
echo "  BOTH, ambiguous:  $(comm -12 $OUT/g97.txt $OUT/a97.txt | wc -l)"
echo ""
echo "=== are the three abundant sequences ambiguous ==="
for S in 1f17ffedbaa6b0d4a32ab3777204c1c7 d723db35121fdf9c889e9d3b2916e09a 58ab061b3a1471d68263fc7f7b7ce590; do
  G=$(awk -F'\t' -v s=$S '$2==s{print $3}' $OUT/hits_genus.tsv | sort -nr | head -1)
  A=$(awk -F'\t' -v s=$S '$2==s{print $3}' $OUT/hits_akk.tsv | sort -nr | head -1)
  printf "  %-34s ours %-8s Akkermansia %s\n" $S "${G:-none}" "${A:-none}"
done
