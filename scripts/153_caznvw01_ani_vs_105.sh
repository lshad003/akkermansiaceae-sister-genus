#!/bin/bash
CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
SKANI=/bigdata/stajichlab/lshad003/condaenvs/drep/bin/skani
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
CEN=$CH3/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv
GDIR=$CH3/data/amphibia_gtdbtk_input
Q=$CH3/results/caznvw_tree/genome/CAZNVW01.fna
OUTD=$CH3/results/caznvw_ani
mkdir -p $OUTD

[ -s "$Q" ] || { echo "REFUSED: missing $Q"; exit 1; }

$PY - << 'PYEOF' > $OUTD/candidate_105.txt
import csv, os
CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CEN = CH3 + "/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"
GDIR = CH3 + "/data/amphibia_gtdbtk_input"
rows = [r for r in csv.DictReader(open(CEN), delimiter="\t")
        if r["family"] == "Akkermansiaceae" and r["annotated"] == "1"]
def gen(r): return (r["genus"] or "").strip()
cand = [r for r in rows if r["host_class"] == "amphibian" and gen(r) in ("", "unknown", "NO_GENUS")]
import sys
miss = []
for r in cand:
    acc = r.get("accession") or r.get("genome")
    p = None
    for ext in (".fa", ".fna", ".fasta"):
        for a in (acc, acc[3:] if acc[:3] in ("GB_", "RS_") else acc):
            q = os.path.join(GDIR, a + ext)
            if os.path.exists(q): p = q; break
        if p: break
    if p: print(p)
    else: miss.append(acc)
sys.stderr.write("candidates: %d   files found: %d   missing: %d %s\n" % (len(cand), len(cand)-len(miss), len(miss), miss[:5]))
PYEOF

N=$(wc -l < $OUTD/candidate_105.txt)
echo "assemblies listed: $N"
if [ "$N" -ne 105 ]; then echo "REFUSED: expected 105"; exit 1; fi

$SKANI dist -q $Q --rl $OUTD/candidate_105.txt -o $OUTD/caznvw01_ani_vs_105.tsv -t 8 --min-af 0
echo "rows: $(($(wc -l < $OUTD/caznvw01_ani_vs_105.tsv) - 1))"
echo
echo "=== top 5 by ANI ==="
sort -t$'\t' -k3,3gr $OUTD/caznvw01_ani_vs_105.tsv | grep -v "^Ref_file" | head -5 | awk -F'\t' '{printf "%-44s ANI %6.2f  AF_ref %5.1f  AF_query %5.1f\n", $7, $3, $4, $5}'
echo
echo "hits at or above 95 percent ANI: $(awk -F'\t' 'NR>1 && $3>=95' $OUTD/caznvw01_ani_vs_105.tsv | wc -l)"
