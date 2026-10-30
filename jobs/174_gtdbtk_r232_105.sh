#!/bin/bash
#SBATCH --job-name=r232cls
#SBATCH --partition=batch
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=180G
#SBATCH --output=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/r232cls_%j.out
#SBATCH --error=/bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/r232cls_%j.err
# GTDB-Tk r232 classify_wf on the 105 candidate genomes plus CAZNVW01 as a positive control.
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
LIST=$CH3/results/r232_new_genera/novel105.txt
CTRL=$CH3/results/r232_new_genera/fna/GCA_964576595.1.fna
OUT=$CH3/results/gtdbtk_r232_105
cd "$CH3" || exit 1

if [ -e "$OUT" ]; then echo "REFUSED: $OUT already exists"; exit 1; fi

module load gtdbtk
VER=$(gtdbtk --version 2>&1 | grep -oE "[0-9]+\.[0-9]+\.[0-9]+" | head -1)
MAJ=${VER%%.*}; REST=${VER#*.}; MIN=${REST%%.*}
if [ -z "$VER" ] || [ "$MAJ" -lt 2 ] || { [ "$MAJ" -eq 2 ] && [ "$MIN" -lt 7 ]; }; then
  echo "REFUSED: GTDB-Tk $VER cannot read R232"; exit 1
fi
if ! grep -q "^VERSION_DATA=r232" "$GTDBTK_DATA_PATH/metadata/metadata.txt" 2>/dev/null; then
  echo "REFUSED: $GTDBTK_DATA_PATH is not r232"; exit 1
fi
echo "GTDB-Tk $VER, data $GTDBTK_DATA_PATH"

N=$(grep -c . "$LIST")
if [ "$N" -ne 105 ]; then echo "REFUSED: $LIST has $N genomes, expected 105"; exit 1; fi
MISS=0
while read -r p; do [ -z "$p" ] && continue; [ -s "$p" ] || { echo "MISSING $p"; MISS=$((MISS+1)); }; done < "$LIST"
[ -s "$CTRL" ] || { echo "MISSING $CTRL"; MISS=$((MISS+1)); }
if [ "$MISS" -ne 0 ]; then echo "REFUSED: $MISS input FASTAs missing"; exit 1; fi

mkdir -p "$OUT"
BATCH=$OUT/batchfile.tsv
while read -r p; do
  [ -z "$p" ] && continue
  b=$(basename "$p"); id=${b%.fa}; id=${id%.fna}; id=${id%.fasta}
  printf "%s\t%s\n" "$p" "$id"
done < "$LIST" > "$BATCH"
printf "%s\t%s\n" "$CTRL" "control_CAZNVW01" >> "$BATCH"
echo "batchfile rows: $(wc -l < "$BATCH") (expect 106)"

gtdbtk classify_wf --batchfile "$BATCH" --out_dir "$OUT/classify" --cpus 16 --pplacer_cpus 1

if [ ! -s "$OUT/classify/gtdbtk.bac120.summary.tsv" ]; then echo "FAILED: no bac120 summary"; exit 1; fi

$PY - "$OUT" << 'PYEOF'
import csv, sys, collections
out = sys.argv[1]
TYPE = "UHM1073.23039_R.bin.125"
CTRL = "control_CAZNVW01"
rows = list(csv.DictReader(open(out + "/classify/gtdbtk.bac120.summary.tsv"), delimiter="\t"))
print("rows:", len(rows), "(expect 106: 105 + control)")
def genus(c):
    for f in c.split(";"):
        if f.startswith("g__"):
            return f
    return "g__(none)"
def fl(x):
    try:
        return float(x)
    except Exception:
        return None
calls = []
for r in rows:
    cl = r.get("classification", "NA") or "NA"
    calls.append((r["user_genome"], genus(cl), r.get("classification_method", "NA") or "NA",
                  r.get("red_value", "NA") or "NA", cl))
by = {c[0]: c for c in calls}
for name in (CTRL, TYPE):
    c = by.get(name)
    if c:
        print(name, "->", c[4], "|", c[2], "| RED", c[3])
    else:
        print(name, "-> MISSING")
q = [c for c in calls if c[0] != CTRL]
print("query genomes:", len(q))
for g, n in collections.Counter(c[1] for c in q).most_common():
    reds = [fl(c[3]) for c in q if c[1] == g and fl(c[3]) is not None]
    rr = "RED %.4f-%.4f" % (min(reds), max(reds)) if reds else "RED NA"
    print("  %-22s %3d  %s" % (g if g != "g__" else "g__ (empty)", n, rr))
print("classification_method:", dict(collections.Counter(c[2] for c in q)))
with open(out + "/r232_genus_calls.tsv", "w") as fh:
    fh.write("genome\tgenus\tclassification_method\tred_value\tclassification\n")
    for c in sorted(q):
        fh.write("\t".join(c) + "\n")
print("wrote", out + "/r232_genus_calls.tsv")
PYEOF
echo "DONE"
