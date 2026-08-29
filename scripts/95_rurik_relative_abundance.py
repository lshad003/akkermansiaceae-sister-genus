#!/usr/bin/env python3
# Relative abundance of the genus by host order
# Source: ch3-chitin-evolution/scripts/rurik_relative_abundance.py
# Output: results/rurik_16s/relative_abundance.txt
import csv, sys, collections, openpyxl
CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/"
OUT = CH3 + "results/rurik_16s/relative_abundance.txt"
csv.field_size_limit(sys.maxsize)

wb = openpyxl.load_workbook(CH3+"data/251027_16S_metadata_final.xlsx", read_only=True, data_only=True)
ws = wb.worksheets[0]
rows = ws.iter_rows(values_only=True)
hdr = ['' if c is None else str(c) for c in next(rows)]
I = {h:i for i,h in enumerate(hdr)}
meta = {}
for row in rows:
    c = ['' if x is None else str(x).strip() for x in row]
    if I["sample_name"] < len(c):
        meta[c[I["sample_name"]]] = (
            c[I["Clade_Order"]] if I["Clade_Order"] < len(c) else "",
            c[I["env_broad_scale"]] if I["env_broad_scale"] < len(c) else "",
            (c[I["host_genus"]]+" "+c[I["host_species"]]).strip() if I["host_species"] < len(c) else "")
wb.close()

want = set(l.strip() for l in open(CH3+"results/rurik_16s/matched_ids.txt") if l.strip())
with open(CH3+"data/251001_16S_count_table.csv", newline="") as fh:
    r = csv.reader(fh)
    samples = next(r)[1:]
    depth = [0.0]*len(samples)
    ours  = [0.0]*len(samples)
    for row in r:
        k = row[0].strip().strip('"')
        mine = k in want
        for i, v in enumerate(row[1:]):
            try: x = float(v)
            except ValueError: continue
            depth[i] += x
            if mine: ours[i] += x

L=[]
def say(s):
    print(s); L.append(s)

ORD = ("Caudata","Anura","Squamata","Testudines")
by = collections.defaultdict(list)
for i, s in enumerate(samples):
    m = meta.get(s)
    if not m or m[0] not in ORD or depth[i] <= 0:
        continue
    by[m[0]].append((100.0*ours[i]/depth[i], ours[i], depth[i], s, m[1], m[2]))

say("RELATIVE ABUNDANCE OF THE GENUS, percent of reads per sample")
say("")
say("%-12s %6s %8s %9s %9s %9s %9s" % ("order","n","positive","med all","med pos","mean pos","max"))
for o in ORD:
    v = by[o]
    if not v: continue
    pos = sorted(x[0] for x in v if x[0] > 0)
    alls = sorted(x[0] for x in v)
    say("%-12s %6d %8d %8.3f%% %8.3f%% %8.3f%% %8.2f%%" % (
        o, len(v), len(pos),
        alls[len(alls)//2],
        pos[len(pos)//2] if pos else 0,
        (sum(pos)/len(pos)) if pos else 0,
        max(alls)))
say("")

say("HOW ABUNDANT WHERE PRESENT: samples above 1 percent of the community")
for o in ORD:
    v = [x for x in by[o] if x[0] >= 1.0]
    say("   %-12s %3d samples" % (o, len(v)))
say("")

say("TOP 15 SAMPLES BY RELATIVE ABUNDANCE")
allv = sorted((x for o in ORD for x in by[o]), reverse=True)
say("%-20s %-11s %-26s %-8s %8s %10s" % ("sample","order","host","env","pct","reads"))
for pct, o_reads, d, s, env, host in allv[:15]:
    o = meta[s][0]
    say("%-20s %-11s %-26s %-8s %7.2f%% %10.0f" % (s, o, host[:26], env[:8], pct, o_reads))
say("")

say("TOP 10 REPTILE SAMPLES BY RELATIVE ABUNDANCE")
rep = sorted((x for o in ("Squamata","Testudines") for x in by[o]), reverse=True)
for pct, o_reads, d, s, env, host in rep[:10]:
    say("   %-20s %-26s %-8s %7.3f%% %8.0f reads of %.0f" % (s, host[:26], env[:8], pct, o_reads, d))

open(OUT,"w").write("\n".join(L)+"\n")
print(""); print("WROTE:", OUT)
