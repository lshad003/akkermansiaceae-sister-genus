#!/usr/bin/env python3
# Pathway query proteins identified by sequence motif and length
# Source: ch3-chitin-evolution/scripts/id_ppp_trio.py
# Output: results/pangenome/ppp_trio_gene_identity.txt
P = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/results/pangenome/ppp_trio.faa"

name = None
seq = []
recs = []
for line in open(P):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if name:
            recs.append((name, "".join(seq)))
        name = line[1:].strip()
        seq = []
    else:
        seq.append(line.strip())
if name:
    recs.append((name, "".join(seq)))

print("EXPECTED LENGTHS: zwf ~490-510, gnd ~470-490, g6pd accessory subunit ~200-250")
print("")
for n, s in recs:
    if len(s) < 300:
        call = "g6pd_sub (short)"
    elif len(s) >= 495:
        call = "zwf (long)"
    else:
        call = "gnd (mid)"
    print("%-36s %4d aa   -> %s" % (n, len(s), call))
    print("     first 60 aa: %s" % s[:60])
