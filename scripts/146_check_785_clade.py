import os, glob

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"

cand = []
for pat in ("results/caznvw_tree/**/*.treefile", "results/caznvw_tree/**/*.nwk",
            "results/**/*785*", "results/novel_akk_tree/**/*785*"):
    cand += glob.glob(os.path.join(CH3, pat), recursive=True)
cand = sorted(set(p for p in cand if os.path.isfile(p)))

print("CANDIDATE TREE FILES")
for p in cand:
    print("   %-62s %9d bytes" % (p.replace(CH3 + "/", ""), os.path.getsize(p)))
if not cand:
    print("   NONE FOUND. The 785-tip tree was not run or is elsewhere.")
    raise SystemExit

tf = None
for p in cand:
    t = open(p, errors="replace").read()
    if p.endswith(".treefile"):
        tf = p
        txt = t
        break
if tf is None:
    print("\nNo file looks like a newick tree.")
    raise SystemExit

print("\nUSING", tf.replace(CH3 + "/", ""))

import re
tips = re.findall(r'[(,]\s*([^(),:]+)\s*:', txt)
tips = [t.strip().strip("'\"") for t in tips if t.strip()]
print("tips:", len(tips))

target = [t for t in tips if "964576595" in t or "CAZNVW" in t.upper()]
print("CAZNVW01 tip:", target if target else "NOT FOUND IN LABELS")
if not target:
    raise SystemExit

def parse(s):
    s = s.strip().rstrip(";")
    pos = [0]
    def node():
        if s[pos[0]] == "(":
            pos[0] += 1
            ch = [node()]
            while s[pos[0]] == ",":
                pos[0] += 1
                ch.append(node())
            pos[0] += 1
            j = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",)":
                pos[0] += 1
            return ("i", s[j:pos[0]], ch)
        j = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",)":
            pos[0] += 1
        lab = s[j:pos[0]].split(":")[0].strip().strip("'\"")
        return ("l", lab, [])
    return node()

root = parse(txt)
tgt = target[0]

def leaves(n):
    if n[0] == "l":
        return [n[1]]
    out = []
    for c in n[2]:
        out += leaves(c)
    return out

path = []
def find(n):
    if n[0] == "l":
        return n[1] == tgt
    for c in n[2]:
        if find(c):
            path.append(n)
            return True
    return False
find(root)

print("\nNESTED CLADES CONTAINING CAZNVW01, smallest first")
for n in path[:6]:
    L = leaves(n)
    novel = [x for x in L if "NOVEL" in x.upper() or x.startswith(("UHM", "EHM", "STP"))]
    akk = [x for x in L if "AKK" in x.upper() and "NOVEL" not in x.upper()]
    print("   size %4d   candidate-like %4d   akkermansia-like %4d   support %s"
          % (len(L), len(novel), len(akk), n[1][:22]))
    if len(L) <= 12:
        print("      ", "; ".join(L))
