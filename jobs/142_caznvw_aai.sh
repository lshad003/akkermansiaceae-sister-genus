#!/bin/bash -l
#SBATCH -p short -c 8 --mem 24G -t 01:00:00
#SBATCH -J caznvwaai -o /bigdata/stajichlab/lshad003/ch3-chitin-evolution/logs/caznvwaai_%j.out

CH3=/bigdata/stajichlab/lshad003/ch3-chitin-evolution
W=$CH3/results/r232_new_genera
PY=/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
mkdir -p $W/aai_tmp

module unload diamond 2>/dev/null
module load diamond/2.2.6
diamond --version

Q=$W/faa/GCA_964576595.1.faa
[ -s "$Q" ] || { echo "REFUSED: missing $Q"; exit 1; }
echo "query proteins: $(grep -c '^>' $Q)"

$PY - << 'PYEOF'
import os, csv, subprocess, shutil, itertools

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/r232_new_genera"
T = W + "/aai_tmp"
Q = W + "/faa/GCA_964576595.1.faa"
D = shutil.which("diamond")

DIRS = [CH3 + "/results/pangenome/" + d for d in ("gff199","akkfam","outgroups","gtdb_akk")]
DIRS += [CH3 + "/results/" + d for d in ("dbcan_scaffold","dbcan_ehi_amphibian","dbcan_ehi_nonamph")]
idx = {}
for d in DIRS:
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.endswith(".faa"): idx.setdefault(f[:-4], os.path.join(d, f))

targets = []
for r in csv.DictReader(open(CH3 + "/results/cluster_aai/cluster_representatives.tsv"), delimiter="\t"):
    targets.append(("CAND_" + r["cluster_id"], r["representative"]))
targets.append(("Akkermansia_ref", "RS_GCF_040616545.1"))

def find(a):
    return idx.get(a) or idx.get(a.replace("GB_","").replace("RS_",""))

miss = [g for _, g in targets if not find(g)]
if miss:
    print("REFUSED, no proteome for:", miss); raise SystemExit(1)
print("targets:", len(targets))

def blast(q, s, tag):
    db = "%s/db_%s" % (T, tag); o = "%s/h_%s.tsv" % (T, tag)
    subprocess.run([D,"makedb","--in",s,"-d",db,"--quiet"], check=True)
    subprocess.run([D,"blastp","-q",q,"-d",db,"-o",o,"--outfmt","6","qseqid","sseqid",
        "pident","length","evalue","bitscore","--evalue","1e-5","--max-target-seqs","1",
        "--more-sensitive","--quiet","--threads","8"], check=True)
    best = {}
    for line in open(o):
        f = line.rstrip("\n").split("\t")
        if len(f) < 6: continue
        q_, s_, pid, bs = f[0], f[1], float(f[2]), float(f[5])
        if q_ not in best or bs > best[q_][2]: best[q_] = (s_, pid, bs)
    os.remove(o); os.remove(db + ".dmnd")
    return best

rows = []
for lab, g in targets:
    p = find(g)
    fwd = blast(Q, p, "f"); rev = blast(p, Q, "r")
    ids = [v[1] for q_, v in fwd.items() if rev.get(v[0]) and rev[v[0]][0] == q_]
    aai = sum(ids)/len(ids) if ids else float("nan")
    rows.append((lab, g, aai, len(ids)))
    print("  %-18s AAI %6.2f  RBH %5d" % (lab, aai, len(ids)))

with open(W + "/caznvw_aai.tsv", "w") as o:
    o.write("target_label\ttarget_genome\tAAI_pct\tn_RBH\n")
    for r in rows:
        o.write("%s\t%s\t%.2f\t%d\n" % r)

cand = [r[2] for r in rows if r[0].startswith("CAND_")]
akk = [r[2] for r in rows if r[0] == "Akkermansia_ref"][0]
print()
print("vs 17 candidate reps: %.2f to %.2f   (chat-only values were 60.5 to 84.8)" % (min(cand), max(cand)))
print("vs Akkermansia ref  : %.2f          (chat-only value was 55.7)" % akk)
print("WROTE", W + "/caznvw_aai.tsv")
PYEOF
rm -rf $W/aai_tmp
