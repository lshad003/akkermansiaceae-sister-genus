#!/usr/bin/env python3
# Census patched to a uniform evidence filter
# Source: ch3-chitin-evolution/scripts/patch_census_v3.py
# Output: scripts/gh75_verru_census_v3.py
import re, os, shutil
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SRC=f"{BASE}/scripts/gh75_verru_census_v2.py"
DST=f"{BASE}/scripts/gh75_verru_census_v3.py"
RUM="/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"

s=open(SRC).read()
orig=s
def sub(old,new,label):
    global s
    if old not in s:
        print(f"  FAIL: could not find pattern for [{label}]"); print("  ---\n"+old+"\n  ---"); raise SystemExit(1)
    if s.count(old)!=1:
        print(f"  FAIL: pattern for [{label}] appears {s.count(old)} times, expected 1"); raise SystemExit(1)
    s=s.replace(old,new); print(f"  ok: {label}")

print("PATCHING census v2 -> v3")

# 1. add the herptile raw dbCAN dir (absolute path)
old1='''DBCAN_DIRS = ["dbcan_verru", "dbcan_bact_refs", "dbcan_ehi_amphibian",
              "dbcan_ehi_nonamph",  # ADDED 2026-07-13: the 232 reptile+mammal EHI Verru (job 26223734)
              "dbcan_endo_allphyla", "dbcan_flavo_refs", "dbcan_scaffold"]'''
new1='''DBCAN_DIRS = ["dbcan_verru", "dbcan_bact_refs", "dbcan_ehi_amphibian",
              "dbcan_ehi_nonamph",  # ADDED 2026-07-13: the 232 reptile+mammal EHI Verru (job 26223734)
              "dbcan_endo_allphyla", "dbcan_flavo_refs", "dbcan_scaffold",
              # ADDED 2026-07-15 (v3): raw dbCAN for the 147 herptile_MAG Verru. Chain verified
              # byte-for-byte vs cazy_ehi_amphibian_array.sh: prodigal -p meta, HMMER 3.3.2,
              # CAZyDB v11.0, hmmscan-parser.sh. Replaces the precomputed table branch.
              "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"]

# v3: the filter, applied UNIFORMLY to every genome. v2 applied NO filter to tsv-read genomes
# but inherited E<=1e-10 AND cov>=0.30 for herptile rows via the precomputed table, so the two
# halves of the amphibian set were filtered differently. That is the bug this fixes.
EVAL_MAX = 1e-10
COV_MIN  = 0.30'''
sub(old1,new1,"DBCAN_DIRS + filter constants")

# 2. allow absolute dirs in the index loop
old2='''for d in DBCAN_DIRS:
    p = f"{BASE}/results/{d}"'''
new2='''for d in DBCAN_DIRS:
    p = d if d.startswith("/") else f"{BASE}/results/{d}"'''
sub(old2,new2,"absolute dir support")

# 3. THE FIX: filter GH75 domains by evalue + coverage
old3='''        with open(tsv) as fh:
            for line in fh:
                if not line.startswith("GH75"):   # cheap prefilter before regex
                    continue
                f = line.rstrip("\\n").split("\\t")
                if len(f) < 3 or not GH75_RE.match(f[0]):
                    continue
                n_dom += 1
                prots.add(f[2])'''
new3='''        with open(tsv) as fh:
            for line in fh:
                if not line.startswith("GH75"):   # cheap prefilter before regex
                    continue
                f = line.rstrip("\\n").split("\\t")
                if len(f) < 10 or not GH75_RE.match(f[0]):
                    continue
                # v3: uniform quality filter. hmmscan-parser.sh columns are
                # 0 FAM.hmm  1 hmm_len  2 prot_id  3 prot_len  4 evalue  ...  9 coverage
                try:
                    ev = float(f[4]); cv = float(f[9])
                except ValueError:
                    continue
                if ev > EVAL_MAX or cv < COV_MIN:
                    continue
                n_dom += 1
                prots.add(f[2])'''
sub(old3,new3,"evalue + coverage filter")

# 4. redirect ALL outputs to a new dir
m=re.search(r'^OUTDIR\s*=\s*.*$', s, re.M)
if not m:
    print("  FAIL: no OUTDIR assignment found"); raise SystemExit(1)
print(f"  found: {m.group(0)}")
s=s[:m.start()]+'OUTDIR = f"{BASE}/results/gh75_census_v3"\nos.makedirs(OUTDIR, exist_ok=True)'+s[m.end():]
print("  ok: OUTDIR -> results/gh75_census_v3")

open(DST,"w").write(s)
print(f"\nwrote {DST}")
print("\n=== changed lines ===")
import difflib
for l in difflib.unified_diff(orig.splitlines(), s.splitlines(), lineterm="", n=1):
    if l.startswith(("+","-")) and not l.startswith(("+++","---")):
        print("  "+l)
