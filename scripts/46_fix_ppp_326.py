#!/usr/bin/env python3
# Presence table rebuilt with corrected collection denominators
# Source: ch3-chitin-evolution/scripts/fix_ppp_326.py
# Output: results/pangenome/ppp_presence_verified.tsv
import os
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
P=f"{BASE}/results/pangenome"

# rewrite ppp_presence_verified.tsv with correct denominators from the DB directories
D_AMPH=94
D_AKKFAM=len([f for f in os.listdir(f"{P}/akkfam") if f.endswith(".faa")])   # 172
D_GTDB=len([f for f in os.listdir(f"{P}/gtdb_akk") if f.endswith(".faa")])   # 60
TOTAL=D_AMPH+D_AKKFAM+D_GTDB
print(f"denominators: amph {D_AMPH}, akkfam {D_AKKFAM}, gtdb {D_GTDB}, total {TOTAL}")

# candidate carriers (verified earlier): zwf 100, gnd 101, g6pd 97
cand={"zwf":100,"gnd":101,"g6pd_sub":97}
out=f"{P}/ppp_presence_verified.tsv"
with open(out,"w") as f:
    f.write("gene\tn_cand\td_cand\tn_akkamph\td_akkamph\tn_akkfam\td_akkfam\tn_gtdb\td_gtdb\tn_akk_total\td_akk_total\n")
    for g in ["zwf","gnd","g6pd_sub"]:
        f.write(f"{g}\t{cand[g]}\t105\t0\t{D_AMPH}\t0\t{D_AKKFAM}\t0\t{D_GTDB}\t0\t{TOTAL}\n")
print(f"rewrote {out} with Akkermansia total 0/{TOTAL}")
print(open(out).read())
