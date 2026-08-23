#!/usr/bin/env python3
# Family join patched to the corrected census
# Source: ch3-chitin-evolution/scripts/patch_join_v3.py
# Output: scripts/join_herptile_verru_family_v3.py
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SRC=f"{BASE}/scripts/join_herptile_verru_family_v2.py"
DST=f"{BASE}/scripts/join_herptile_verru_family_v3.py"
s=open(SRC).read()
def sub(old,new,label):
    global s
    if s.count(old)!=1:
        print(f"FAIL [{label}]: found {s.count(old)}, expected 1"); raise SystemExit(1)
    s=s.replace(old,new); print(f"  ok: {label}")
sub('CEN  = f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_v2.tsv"',
    'CEN  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_v3.tsv"',"CEN -> v3")
sub('OUT  = f"{BASE}/results/gh75_census/gh75_verru_census_per_genome_familyfilled_v2.tsv"',
    'OUT  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"',"OUT -> v3")
open(DST,"w").write(s)
print(f"wrote {DST}")
