import csv, itertools, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/delim58"
OUT = W + "/FigureS_aai_pocp_scatter"
TYPE = "CAND_C288"

labs = [r["label"] for r in csv.DictReader(open(W + "/members.tsv"), delimiter="\t")]
def mat(p):
    rows = list(csv.reader(open(p), delimiter="\t")); h = rows[0][1:]
    return {(r[0], h[j]): float(v) for r in rows[1:] for j, v in enumerate(r[1:])}
A = mat(W + "/aai_matrix.tsv"); P = mat(W + "/pocp_matrix.tsv")

def cls(x):
    if x.startswith("CAND_"): return "C"
    if x.startswith("AKK_clade"): return "A"
    return "O"
GRP = [("CC", "candidate / candidate", "#1b7837", 16),
       ("AA", "Akkermansia s.l. / Akkermansia s.l.", "#8c510a", 16),
       ("CA", "candidate / Akkermansia s.l.", "#b2182b", 16),
       ("CO", "candidate / other genera", "#4575b4", 13),
       ("AO", "Akkermansia s.l. / other genera", "#762a83", 11),
       ("OO", "other genera / other genera", "#999999", 11)]

pts = {k: ([], []) for k, _, _, _ in GRP}
for a, b in itertools.combinations(labs, 2):
    t = "".join(sorted(cls(a) + cls(b)))
    k = {"CC": "CC", "AA": "AA", "AC": "CA", "CO": "CO", "AO": "AO", "OO": "OO"}[t]
    pts[k][0].append(A[(a, b)]); pts[k][1].append(P[(a, b)])

fig, ax = plt.subplots(figsize=(5.4, 4.8))
for k, lab, col, s in GRP:
    x, y = pts[k]
    ax.scatter(x, y, s=s, alpha=0.55, linewidth=0, color=col,
               label="%s (n=%d)" % (lab, len(x)))
ax.scatter([A[(TYPE, "AKK_clade01")]], [P[(TYPE, "AKK_clade01")]], s=0, alpha=0)
ax.axvline(65, ls="--", lw=0.9, color="#666666")
ax.axhline(50, ls="--", lw=0.9, color="#666666")
ax.text(65.8, 16, "65% AAI", fontsize=6.5, color="#666666")
ax.text(46.5, 51, "50% POCP", fontsize=6.5, color="#666666")
ax.set_xlabel("AAI (%)", fontsize=9); ax.set_ylabel("POCP (%)", fontsize=9)
ax.tick_params(labelsize=8)
ax.legend(fontsize=6.2, loc="upper left", frameon=False)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
plt.tight_layout()
plt.savefig(OUT + ".pdf", bbox_inches="tight")
plt.savefig(OUT + ".png", dpi=300, bbox_inches="tight")

print("pairs by class:", {k: len(pts[k][0]) for k, _, _, _ in GRP},
      "total", sum(len(pts[k][0]) for k, _, _, _ in GRP))
x, y = pts["CA"]
print("candidate/Akkermansia: AAI %.2f-%.2f  POCP %.2f-%.2f" % (min(x), max(x), min(y), max(y)))
print("CA pairs above 50%% POCP: %d of %d" % (sum(1 for v in y if v >= 50), len(y)))
x0, y0 = pts["CC"]
print("candidate/candidate AAI %.2f-%.2f; below 65%%: %d of %d"
      % (min(x0), max(x0), sum(1 for v in x0 if v < 65), len(x0)))
x2, y2 = pts["AA"]
print("Akk/Akk pairs below 50%% POCP: %d of %d" % (sum(1 for v in y2 if v < 50), len(y2)))
print("WROTE", OUT + ".pdf and .png")
