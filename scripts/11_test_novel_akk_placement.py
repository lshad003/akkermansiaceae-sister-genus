#!/usr/bin/env python3
# Monophyly, exclusion, and sister relationship tested
# Source: ch3-chitin-evolution/scripts/test_novel_akk_placement.py
# Output: stdout
import re, collections
BASE="/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
NWK=f"{BASE}/results/novel_akk_tree/akk_placement.nwk"
s=open(NWK).read().strip()

# parse newick into nested lists of tip labels
def parse(s):
    i=0
    def node():
        nonlocal i
        kids=[]
        if s[i]=="(":
            i+=1
            while True:
                kids.append(node())
                if s[i]==",": i+=1
                else: break
            i+=1  # ')'
        lab=""
        while i<len(s) and s[i] not in ",():;":
            lab+=s[i]; i+=1
        if s[i:i+1]==":":
            i+=1
            while i<len(s) and s[i] not in ",():;": i+=1
        return (kids,lab.strip())
    return node()
root=parse(s)

def tips(n):
    k,l=n
    if not k: return [l]
    out=[]
    for c in k: out+=tips(c)
    return out
def support(n):
    k,l=n
    if not k: return None
    try: return float(l)
    except: return None

allt=tips(root)
def grp(t): return t.split("|")[0]
print("tips:",len(allt))
print(collections.Counter(grp(t) for t in allt).most_common(6))

novel={t for t in allt if grp(t)=="NOVEL"}
akk={t for t in allt if grp(t)=="Akkermansia"}
print(f"\nNOVEL {len(novel)}  Akkermansia {len(akk)}")

# every clade
clades=[]
def walk(n):
    k,l=n
    if not k: return
    clades.append((set(tips(n)), support(n)))
    for c in k: walk(c)
walk(root)

# 1. smallest clade containing all NOVEL
cand=[(c,sp) for c,sp in clades if novel <= c]
cand.sort(key=lambda x: len(x[0]))
c,sp=cand[0]
intr=c-novel
print(f"\n=== Q1: are the 105 NOVEL monophyletic? ===")
print(f"  smallest clade containing all NOVEL: {len(c)} tips, support {sp}")
print(f"  intruders: {len(intr)}")
if intr:
    print("  intruder groups:",collections.Counter(grp(t) for t in intr).most_common())
else:
    print("  NOVEL IS MONOPHYLETIC")

# 2. is that clade inside Akkermansia, or sister to it?
cand2=[(cc,ss) for cc,ss in clades if akk <= cc]
cand2.sort(key=lambda x: len(x[0]))
ca,sa=cand2[0]
print(f"\n=== Q2: NOVEL vs Akkermansia ===")
print(f"  smallest clade containing all Akkermansia: {len(ca)} tips, support {sa}")
print(f"  NOVEL tips inside it: {len(novel & ca)}/{len(novel)}")
print(f"  non-Akkermansia in it:",collections.Counter(grp(t) for t in ca-akk).most_common(5))

# 3. sister of the NOVEL clade
def find_parent(target):
    best=None
    for cc,ss in clades:
        if target < cc and (best is None or len(cc)<len(best[0])): best=(cc,ss)
    return best
p=find_parent(c)
if p:
    sis=p[0]-c
    print(f"\n=== Q3: sister group of the NOVEL clade ===")
    print(f"  parent {len(p[0])} tips support {p[1]}; sister {len(sis)} tips")
    print("  sister composition:",collections.Counter(grp(t) for t in sis).most_common(5))
