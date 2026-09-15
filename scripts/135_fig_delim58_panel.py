import csv, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ete3 import Tree

CH3 = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
W = CH3 + "/results/delim58"
OUT = W + "/Figure_genus_delimitation.pdf"
OUTPNG = W + "/Figure_genus_delimitation.png"
TYPE = "CAND_C288"

order = [r["label"] for r in csv.DictReader(open(W + "/delim58_tree_order.tsv"), delimiter="\t")]
n = len(order)

def mat(path):
    rows = list(csv.reader(open(path), delimiter="\t"))
    hdr = rows[0][1:]
    D = {(r[0], hdr[j]): float(v) for r in rows[1:] for j, v in enumerate(r[1:])}
    return np.array([[D[(a, b)] for b in order] for a in order])

AAI = mat(W + "/aai_matrix.tsv")
POCP = mat(W + "/pocp_matrix.tsv")

t = Tree(W + "/delim58_tree.nwk", format=1)
ypos = {lab: i for i, lab in enumerate(order)}
for nd in t.traverse("preorder"):
    nd.x = (nd.up.x + nd.dist) if nd.up else 0.0
for nd in t.traverse("postorder"):
    nd.y = ypos[nd.name] if nd.is_leaf() else sum(c.y for c in nd.children)/len(nd.children)

def blk(pref):
    ix = [ypos[l] for l in order if l.startswith(pref)]
    return min(ix), max(ix)
c0, c1 = blk("CAND_")
a0, a1 = blk("AKK_clade")

fig = plt.figure(figsize=(13.5, 6.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0], wspace=0.30)

axT = fig.add_subplot(gs[0, 0])
for nd in t.traverse():
    if nd.up:
        axT.plot([nd.up.x, nd.x], [nd.y, nd.y], lw=0.7, color="#333333")
    if not nd.is_leaf():
        ys = [c.y for c in nd.children]
        axT.plot([nd.x, nd.x], [min(ys), max(ys)], lw=0.7, color="#333333")
xm = max(nd.x for nd in t.traverse())
for nd in t.iter_leaves():
    lab = nd.name
    col = "#1b7837" if lab.startswith("CAND_") else ("#8c510a" if lab.startswith("AKK_clade") else "#000000")
    txt = lab.replace("CAND_", "").replace("AKK_clade", "Akk ")
    axT.plot([nd.x, xm * 1.02], [nd.y, nd.y], lw=0.4, ls=":", color="#999999")
    axT.text(xm * 1.04, nd.y, txt + (" (type)" if lab == TYPE else ""),
             fontsize=5.6, va="center", color=col,
             fontweight="bold" if lab == TYPE else "normal")
axT.add_patch(plt.Rectangle((0, c0-0.5), xm*1.03, c1-c0+1, facecolor="#1b7837", alpha=0.07, zorder=0))
axT.add_patch(plt.Rectangle((0, a0-0.5), xm*1.03, a1-a0+1, facecolor="#8c510a", alpha=0.07, zorder=0))
try:
    sup = {r["node"]: r["support"] for r in csv.DictReader(open(W + "/node_support.tsv"), delimiter="\t")}
except Exception:
    sup = {}
cand_n = t.get_common_ancestor([l for l in t.iter_leaves() if l.name.startswith("CAND_")])
akk_n  = t.get_common_ancestor([l for l in t.iter_leaves() if l.name.startswith("AKK_clade")])
shar_n = t.get_common_ancestor([l for l in t.iter_leaves()
                                if l.name.startswith(("CAND_", "AKK_clade"))])
for nd, tag, col in ((cand_n, "CAND", "#1b7837"), (akk_n, "AKK", "#8c510a"),
                     (shar_n, "SHARED", "#000000")):
    v = sup.get(tag, "NA")
    if v == "NA":
        continue
    axT.scatter([nd.x], [nd.y], s=16, color=col, zorder=5)
    axT.text(nd.x - xm*0.015, nd.y, v, fontsize=5.4, ha="right", va="center", color=col)

sb = round(xm/5, 2)
axT.plot([0, sb], [n+0.8, n+0.8], lw=1.4, color="#000000")
axT.text(sb/2, n+1.9, "%.2f subs/site" % sb, fontsize=6, ha="center")
axT.set_ylim(n+3, -1.2); axT.set_xlim(-xm*0.02, xm*1.45)
axT.axis("off")
axT.set_title("A  bac120 phylogeny", fontsize=9, loc="left")

def heat(ax, M, title, cmap, vmin, vmax, thr, thrlab):
    im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks([]); ax.set_yticks(range(n))
    ax.set_yticklabels([l.replace("CAND_", "").replace("AKK_clade", "Akk ") for l in order], fontsize=4.6)
    for (b0, b1, col) in ((c0, c1, "#1b7837"), (a0, a1, "#8c510a")):
        ax.add_patch(plt.Rectangle((b0-0.5, b0-0.5), b1-b0+1, b1-b0+1,
                     fill=False, edgecolor=col, lw=1.2))
    ax.set_title(title, fontsize=9, loc="left")
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6)
    cb.ax.axhline(thr, color="white", lw=1.6)
    cb.ax.axhline(thr, color="#b2182b", lw=0.9, ls="--")
    cb.ax.text(4.2, thr, thrlab, fontsize=5.6, color="#b2182b", va="center",
               transform=cb.ax.get_yaxis_transform())

heat(fig.add_subplot(gs[0, 1]), AAI, "B  AAI (%)", "viridis", 45, 100, 65, "65%")
heat(fig.add_subplot(gs[0, 2]), POCP, "C  POCP (%)", "magma", 25, 95, 50, "50%")

plt.savefig(OUT, bbox_inches="tight")
plt.savefig(OUTPNG, dpi=300, bbox_inches="tight")
print("WROTE", OUT)
print("WROTE", OUTPNG)
