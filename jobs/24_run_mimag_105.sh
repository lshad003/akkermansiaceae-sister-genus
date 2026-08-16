#!/bin/bash
# rRNA and tRNA features called for the candidate set
# Source: ch3-chitin-evolution/scripts/run_mimag_105.sh
# Output: results/mimag

#SBATCH --job-name=mimag
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/mimag3_%A_%a.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/mimag3_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --array=1-105

export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/hmmer/3.3.2/bin:$PATH
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin:$PATH
which nhmmer   || { echo "FATAL no nhmmer"; exit 1; }
which bedtools || { echo "FATAL no bedtools"; exit 1; }
BAR=/opt/linux/rocky/8.x/x86_64/pkgs/barrnap/0.9/bin/barrnap
TRN=/opt/linux/rocky/8.x/x86_64/pkgs/trnascan-se/2.0.12/bin/tRNAscan-SE
FA=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data/amphibia_gtdbtk_input
O=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/mimag
mkdir -p $O
L=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/akk_composition/novel_akk_genome_list.txt
G=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $L | xargs -n1 basename | sed 's/\.fa$//')
[ -z "$G" ] && exit 0
echo "[$(date)] $G"
$BAR --kingdom bac "$FA/${G}.fa" > $O/${G}.rrna.gff
echo "barrnap EXIT=$?"
rm -f $O/${G}.trna.txt
$TRN -B -Q -o $O/${G}.trna.txt "$FA/${G}.fa" > /dev/null
echo "trnascan EXIT=$?"
echo "RESULT $G 5S=$(grep -c '5S_rRNA' $O/${G}.rrna.gff) 16S=$(grep -c '16S_rRNA' $O/${G}.rrna.gff) 23S=$(grep -c '23S_rRNA' $O/${G}.rrna.gff)"
