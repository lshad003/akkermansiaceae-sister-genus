#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Threshold sensitivity for the five polarity categories over ALL orthogroups.
# Recomputes occupancy from Orthogroups.GeneCount.tsv with the same grouping as repo script 37.
import sys, csv
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
GC = B + "/results/pangenome/orthofinder_out/run1/Results_novelakk/Orthogroups/Orthogroups.GeneCount.tsv"
OUT = B + "/results/pangenome/polarity_threshold_sweep.tsv"

r = csv.reader(open(GC), delimiter="\t"); hdr = next(r)
sp = hdr[1:-1] if hdr[-1].lower().startswith("total") else hdr[1:]
def grp(s):
    if s.startswith("NOVEL"): return "NOVEL"
    if s.startswith("AKK"): return "AKK"
    return "FREE"
def genus(s): return s.split("_rep")[0]
cols = {g: [i for i, s in enumerate(sp) if grp(s) == g] for g in ("NOVEL", "AKK", "FREE")}
fgen = sorted(set(genus(s) for s in sp if grp(s) == "FREE"))
gidx = {g: [i for i, s in enumerate(sp) if grp(s) == "FREE" and genus(s) == g] for g in fgen}
print("proteomes: NOVEL %d | AKK %d | FREE %d (%d genera)" % (len(cols["NOVEL"]), len(cols["AKK"]), len(cols["FREE"]), len(fgen)))

rows = []
for row in r:
    vals = [int(x) for x in row[1:1 + len(sp)]]
    pres = [1 if v > 0 else 0 for v in vals]
    pn = 100.0 * sum(pres[i] for i in cols["NOVEL"]) / len(cols["NOVEL"])
    pa = 100.0 * sum(pres[i] for i in cols["AKK"]) / len(cols["AKK"])
    ng = sum(1 for g in fgen if any(pres[i] for i in gidx[g]))
    rows.append((pn, pa, ng))
print("orthogroups read:", len(rows))
if len(rows) != 8812: sys.exit("REFUSED: expected 8812 orthogroups")

def count(HI, LO, ANC, FA, FB):
    A  = sum(1 for n, a, g in rows if n >= HI and a <= LO and g >= FA)
    Bc = sum(1 for n, a, g in rows if n <= ANC and a <= ANC and g >= FB)
    C  = sum(1 for n, a, g in rows if a >= HI and n <= LO)
    D  = sum(1 for n, a, g in rows if n >= HI and a <= LO and g == 0)
    SH = sum(1 for n, a, g in rows if n >= HI and a >= HI)
    gap = sum(1 for n, a, g in rows if n >= HI and a <= LO and 0 < g < FA)
    return A, Bc, C, D, SH, gap

base = (70, 10, 5, 3, 5)
bA, bB, bC, bD, bSH, bgap = count(*base)
if (bA, bB, bC, bD, bSH) != (57, 1047, 34, 55, 1026):
    sys.exit("REFUSED: baseline gives A=%d B=%d C=%d D=%d SH=%d, expected 57/1047/34/55/1026" % (bA, bB, bC, bD, bSH))
print("baseline reproduced. A/D box orthogroups with 1 or 2 free-living genera (in neither A nor D): %d" % bgap)

grid = [("baseline", base)]
for lab, HI in (("HI=60", 60), ("HI=80", 80)):   grid.append((lab, (HI, 10, 5, 3, 5)))
for lab, LO in (("LO=5", 5), ("LO=15", 15)):     grid.append((lab, (70, LO, 5, 3, 5)))
grid.append(("ANC=10", (70, 10, 10, 3, 5)))
for lab, FA in (("FA=2", 2), ("FA=4", 4)):       grid.append((lab, (70, 10, 5, FA, 5)))
for lab, FB in (("FB=4", 4), ("FB=6", 6)):       grid.append((lab, (70, 10, 5, 3, FB)))
grid.append(("all_loose",  (60, 15, 10, 2, 4)))
grid.append(("all_strict", (80, 5, 5, 4, 6)))

hdr = ["setting", "HI_present", "LO_absent", "ANC_absent", "free_min_A", "free_min_B",
       "A_akk_branch", "B_shared_ancestor", "C_akk_enriched", "D_cand_specific", "retained", "gap_1to2_free", "B_over_A"]
with open(OUT, "w") as o:
    o.write("\t".join(hdr) + "\n")
    print("\n%-11s %3s %3s %4s %3s %3s | %5s %6s %4s %4s %6s %5s %6s" % ("setting", "HI", "LO", "ANC", "FA", "FB", "A", "B", "C", "D", "ret", "gap", "B/A"))
    for lab, p in grid:
        A, Bc, C, D, SH, gap = count(*p)
        ratio = Bc / A if A else float("nan")
        o.write("\t".join(map(str, [lab, *p, A, Bc, C, D, SH, gap, "%.1f" % ratio])) + "\n")
        print("%-11s %3d %3d %4d %3d %3d | %5d %6d %4d %4d %6d %5d %6.1f" % (lab, *p, A, Bc, C, D, SH, gap, ratio))
print("\nwrote", OUT)
