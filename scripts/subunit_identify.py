#!/usr/bin/env python3
# What is the 367 aa "G6PD accessory subunit"? Print it for external lookup and
# report what the existing annotations in this project already say about it.
import os, glob, subprocess
from collections import Counter, defaultdict

R = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = R + "/results/ppp_unified"
OUT = W + "/subunit_identity.txt"
QID = "EHM058980_CDS_1432"

fh = open(OUT, "w")
def say(m=""):
    print(m); fh.write(str(m) + "\n"); fh.flush()

# 1. the sequence
seq, keep = [], False
for line in open(R + "/results/pangenome/ppp_trio.faa"):
    if line.startswith(">"):
        keep = QID in line
    elif keep: seq.append(line.strip())
s = "".join(seq)
say("=" * 70); say("THE QUERY"); say("=" * 70)
say(">%s  %d aa" % (QID, len(s)))
for i in range(0, len(s), 60): say(s[i:i+60])

# 2. does any project annotation name it
say(""); say("=" * 70); say("WHAT OUR OWN ANNOTATIONS CALL IT"); say("=" * 70)
hits = 0
for pat in ("/results/**/*.tsv", "/results/**/*.annotations", "/data/**/*.tsv"):
    for p in glob.glob(R + pat, recursive=True):
        if os.path.getsize(p) > 3e9: continue
        try:
            for line in open(p, errors="replace"):
                if QID in line:
                    say("  %s" % p.replace(R + "/", ""))
                    say("     " + line.rstrip()[:260]); hits += 1
                    break
        except OSError: pass
        if hits > 12: break
    if hits > 12: break
if not hits: say("  no annotation row mentions this protein id")

# 3. eggNOG / KEGG assignment if present
say(""); say("=" * 70); say("SEARCHING eggNOG OR KOFAM OUTPUT"); say("=" * 70)
found = False
for p in glob.glob(R + "/results/**/*emapper*", recursive=True) + \
         glob.glob(R + "/results/**/*kofam*", recursive=True) + \
         glob.glob(R + "/results/kofam_nonbac/*", recursive=True)[:3]:
    say("  candidate file: " + p.replace(R + "/", "")); found = True
if not found: say("  none found under results/")

# 4. self-consistency: is the subunit family distinct from zwf?
say(""); say("=" * 70); say("IS THE SUBUNIT FAMILY DISTINCT FROM zwf?"); say("=" * 70)
D = "/bigdata/stajichlab/lshad003/condaenvs/MYNEWENV/bin/diamond"
q = W + "/_pair.faa"
with open(q, "w") as f:
    keep = None
    for line in open(R + "/results/pangenome/ppp_trio.faa"):
        if line.startswith(">"):
            keep = ("zwf" if "UHM979.41089_R.bin.103_CDS_0654" in line else
                    "sub" if QID in line else None)
            if keep: f.write(">%s\n" % keep)
        elif keep: f.write(line)
r = subprocess.run([D, "blastp", "-q", q, "-d", q.replace(".faa", ""), "--quiet"],
                   capture_output=True, text=True)
subprocess.run([D, "makedb", "--in", q, "-d", W + "/_pair", "--quiet"], check=True)
r = subprocess.run([D, "blastp", "-q", q, "-d", W + "/_pair", "--quiet", "--very-sensitive",
                    "-e", "10", "-f", "6", "qseqid", "sseqid", "pident", "length", "evalue"],
                   capture_output=True, text=True)
for line in r.stdout.splitlines():
    say("  " + line)
say("")
say("  If zwf and sub do not align to each other even at evalue 10, they are")
say("  unrelated proteins and the OG01 fusion last night was a clustering artifact.")

say(""); say("=" * 70); say("WHAT YOU NEED TO DO"); say("=" * 70)
say("Paste the sequence above into:")
say("  1. NCBI BLASTp, nr, https://blast.ncbi.nlm.nih.gov")
say("  2. InterProScan, https://www.ebi.ac.uk/interpro/search/sequence/")
say("Look for: OpcA, G6PD assembly protein, DUF1537, or a named domain.")
say("If it returns OpcA, cite Hagen and Meeks 2001 and call it opcA.")
say("If it returns nothing characterised, call it a co-located hypothetical")
say("protein and drop the words accessory subunit, which assert a function.")
for t in (q, W + "/_pair.dmnd"):
    if os.path.exists(t): os.remove(t)
say(""); say("written to " + OUT)
fh.close()
