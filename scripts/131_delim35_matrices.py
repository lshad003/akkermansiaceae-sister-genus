import os, csv, itertools

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/delim35"
H = W + "/hits"

mem = list(csv.DictReader(open(W + "/members.tsv"), delimiter="\t"))
labs = [r["label"] for r in mem]
print("genomes:", len(labs))

nprot = {}
for r in mem:
    n = 0
    for line in open(r["path"], errors="replace"):
        if line[:1] == ">": n += 1
    nprot[r["label"]] = n

def load(a, b):
    p = "%s/%s__%s.tsv" % (H, a, b)
    qual = set(); best = {}
    for line in open(p, errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7: continue
        q, s = f[0], f[1]
        pid = float(f[2]); alen = int(f[3]); qlen = int(f[4]); bs = float(f[6])
        if pid > 40.0 and qlen > 0 and (alen / qlen) > 0.5:
            qual.add(q)
        cur = best.get(q)
        if cur is None or bs > cur[2]:
            best[q] = (s, pid, bs)
    return qual, best

POCP = {}; AAI = {}; NRBH = {}
for a, b in itertools.combinations(labs, 2):
    qa, ba = load(a, b)
    qb, bb = load(b, a)
    p = 100.0 * (len(qa) + len(qb)) / (nprot[a] + nprot[b])
    POCP[(a, b)] = POCP[(b, a)] = p
    ids = []
    for q, (s, pid, bs) in ba.items():
        r = bb.get(s)
        if r and r[0] == q:
            ids.append(pid)
    AAI[(a, b)] = AAI[(b, a)] = (sum(ids) / len(ids)) if ids else float("nan")
    NRBH[(a, b)] = NRBH[(b, a)] = len(ids)

def w(path, D):
    with open(path, "w") as o:
        o.write("label\t" + "\t".join(labs) + "\n")
        for a in labs:
            row = ["%.2f" % D[(a, b)] if a != b else "100.00" for b in labs]
            o.write(a + "\t" + "\t".join(row) + "\n")

w(W + "/aai_matrix.tsv", AAI)
w(W + "/pocp_matrix.tsv", POCP)

TYPE = "CAND_C288"
print()
print("VALIDATION, type genome %s vs Akkermansia" % TYPE)
print("  POCP  %.2f   (pocp_fixed.tsv type anchor = 52.52)" % POCP[(TYPE, "Akkermansia")])
print("  AAI   %.2f over %d RBH   (TableS4 = 56.80 over 1442, DIAMOND 2.1.24)"
      % (AAI[(TYPE, "Akkermansia")], NRBH[(TYPE, "Akkermansia")]))
d = abs(POCP[(TYPE, "Akkermansia")] - 52.52)
print("  delta %.2f  %s" % (d, "OK" if d < 0.6 else "REFUSED, does not reproduce pocp_fixed"))

print()
print("TYPE GENOME vs EVERY GENUS, sorted by AAI")
rows = [(b, AAI[(TYPE, b)], NRBH[(TYPE, b)], POCP[(TYPE, b)])
        for b in labs if not b.startswith("CAND_")]
for b, a, n, p in sorted(rows, key=lambda z: -z[1]):
    print("  %-18s AAI %6.2f  RBH %5d  POCP %6.2f" % (b, a, n, p))

cand = [x for x in labs if x.startswith("CAND_")]
wi = [AAI[(a, b)] for a, b in itertools.combinations(cand, 2)]
print()
print("WITHIN-CANDIDATE AAI: %d pairs, %.2f to %.2f" % (len(wi), min(wi), max(wi)))
print("WROTE aai_matrix.tsv and pocp_matrix.tsv")
