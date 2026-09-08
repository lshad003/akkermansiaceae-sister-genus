#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
import csv, os
from collections import Counter
D = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution/data"
for f in ("UHM_Akkermansia.losses.tsv", "UHM_Akkermansia.core.tsv"):
    p = os.path.join(D, f)
    if not os.path.exists(p):
        print("MISSING", p); continue
    r = list(csv.DictReader(open(p), delimiter="\t"))
    print("\n===== %s : %d rows" % (f, len(r)))
    print("columns:", list(r[0].keys()) if r else "empty")
    if "losses" in f:
        print("ingroup TBLASTN 'none' (real absence):", sum(1 for x in r if x["Ingroup TBLASTN"] == "none"))
        print("ingroup TBLASTN has hits (prediction gap):", sum(1 for x in r if x["Ingroup TBLASTN"] != "none"))
        pf = Counter()
        for x in r:
            for d in (x.get("Pfam domains") or "").split(","):
                if d.strip(): pf[d.strip()] += 1
        print("top Pfam among losses:", pf.most_common(12))
        print("\nlosses with NO ingroup TBLASTN hit, annotated:")
        for x in r:
            if x["Ingroup TBLASTN"] == "none" and (x.get("Pfam domains") or x.get("Product")):
                print("   %-52s %-12s %s" % (x["Protein ID"][:52], x["Outgroup breadth"], (x.get("Pfam domains") or x.get("Product"))[:60]))
    else:
        print("first 3 rows:")
        for x in r[:3]: print("  ", {k: v for k, v in list(x.items())[:6]})
faa = D + "/UHM_Akkermansia.filtered.faa"
if os.path.exists(faa):
    n = sum(1 for l in open(faa) if l.startswith(">"))
    src = Counter(l[1:].split("|")[0] for l in open(faa) if l.startswith(">"))
    print("\n===== filtered.faa : %d proteins" % n)
    print("by source genome:", dict(src))
