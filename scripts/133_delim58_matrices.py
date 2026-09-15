import os, csv, itertools

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/delim58"
H = W + "/hits"

mem = list(csv.DictReader(open(W + "/members.tsv"), delimiter="\t"))
labs = [r["label"] for r in mem]
print("genomes:", len(labs))

nprot = {}
for r in mem:
    nprot[r["label"]] = sum(1 for line in open(r["path"], errors="replace") if line[:1] == ">")

def load(a, b):
    qual = set(); best = {}
    for line in open("%s/%s__%s.tsv" % (H, a, b), errors="replace"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7: continue
        q, s = f[0], f[1]
        pid = float(f[2]); alen = int(f[3]); qlen = int(f[4]); bs = float(f[6])
        ok = pid > 40.0 and qlen > 0 and (alen / qlen) > 0.5
        cur = best.get(q)
        if cur is None or bs > cur[2]:
            best[q] = (s, pid, bs, ok)
    qual = set(q for q in best if best[q][3])
    return qual, best

POCP = {}; AAI = {}; NRBH = {}
for a, b in itertools.combinations(labs, 2):
    qa, ba = load(a, b)
    qb, bb = load(b, a)
    POCP[(a,b)] = POCP[(b,a)] = 100.0*(len(qa)+len(qb))/(nprot[a]+nprot[b])
    ids = [v[1] for q, v in ba.items() if bb.get(v[0]) and bb[v[0]][0] == q]
    AAI[(a,b)] = AAI[(b,a)] = (sum(ids)/len(ids)) if ids else float("nan")
    NRBH[(a,b)] = NRBH[(b,a)] = len(ids)

def w(path, D):
    with open(path, "w") as o:
        o.write("label\t" + "\t".join(labs) + "\n")
        for a in labs:
            o.write(a + "\t" + "\t".join("100.00" if a==b else "%.2f"%D[(a,b)] for b in labs) + "\n")
w(W + "/aai_matrix.tsv", AAI)
w(W + "/pocp_matrix.tsv", POCP)

TYPE = "CAND_C288"
cand = [x for x in labs if x.startswith("CAND_")]
akk  = [x for x in labs if x.startswith("AKK_clade")]
free = [x for x in labs if not x.startswith(("CAND_","AKK_clade"))]

def rng(pairs, D):
    v = [D[p] for p in pairs]
    v = [x for x in v if x == x]
    v.sort()
    return len(v), v[0], v[len(v)//2], v[-1]

print()
print("VALIDATION")
xa = [AAI[(TYPE,b)] for b in akk]; xp = [POCP[(TYPE,b)] for b in akk]
print("  type vs 24 Akkermansia reps: AAI %.2f to %.2f, POCP %.2f to %.2f"
      % (min(xa), max(xa), min(xp), max(xp)))
print("  delim35 single Akk rep gave: AAI 56.77, POCP 52.64")

print()
for name, grp in (("within CANDIDATE", cand), ("within AKKERMANSIA", akk)):
    n,lo,md,hi = rng(list(itertools.combinations(grp,2)), AAI)
    n2,lo2,md2,hi2 = rng(list(itertools.combinations(grp,2)), POCP)
    print("%-20s %3d pairs  AAI %.2f-%.2f (med %.2f)   POCP %.2f-%.2f (med %.2f)"
          % (name, n, lo, hi, md, lo2, hi2, md2))

cross = [(a,b) for a in cand for b in akk]
n,lo,md,hi = rng(cross, AAI); n2,lo2,md2,hi2 = rng(cross, POCP)
print("%-20s %3d pairs  AAI %.2f-%.2f (med %.2f)   POCP %.2f-%.2f (med %.2f)"
      % ("CAND vs AKK", n, lo, hi, md, lo2, hi2, md2))

cf = [(a,b) for a in cand for b in free]
n,lo,md,hi = rng(cf, AAI)
print("%-20s %3d pairs  AAI %.2f-%.2f (med %.2f)" % ("CAND vs free-living", n, lo, hi, md))
print()
print("WROTE aai_matrix.tsv and pocp_matrix.tsv")
