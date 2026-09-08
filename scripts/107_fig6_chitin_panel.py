#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Figure 6: chitin-related CAZyme prevalence across five lineage-matched groups.
# Input:  results/novel_akk_tree/chitin_panel_5groups.tsv
# Output: results/figures/Figure6_chitin_panel.pdf and .png
import csv, os, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SRC = B + "/results/novel_akk_tree/chitin_panel_5groups.tsv"
OUTD = B + "/results/figures"
os.makedirs(OUTD, exist_ok=True)
rows = [{k: (v.strip() if isinstance(v, str) else v) for k, v in r.items()}
        for r in csv.DictReader(open(SRC), delimiter="\t")]
if len(rows) != 15:
    print("REFUSED: expected 15 families, got %d" % len(rows)); sys.exit(1)

COLS = [("pct_novel", "Candidate genus\n(105)"),
        ("pct_akk_amph", "Akkermansia\namphibian (94)"),
        ("pct_podarcis", "Akkermansia\nreptile (137)"),
        ("pct_mammal", "Akkermansia\nmammal (71)"),
        ("pct_gtdb", "Akkermansia\nGTDB (42)")]
MODNAME = {"A": "attack", "B": "deacetylation", "C": "chitosan",
           "D": "terminal", "control": "control"}

# drop families that are zero in every group, keep CE11 last as the control
keep = [r for r in rows if any(float(r[c]) > 0 for c, _ in COLS)]
ctrl = [r for r in keep if r["module"].startswith("control")]
body = [r for r in keep if not r["module"].startswith("control")]
body.sort(key=lambda r: -float(r["pct_novel"]))
order = body + ctrl
print("families plotted: %d of 15 (all-zero rows dropped: %s)" %
      (len(order), [r["family"] for r in rows if r not in keep]))

M = np.array([[float(r[c]) for c, _ in COLS] for r in order])
fig, ax = plt.subplots(figsize=(7.6, 0.42 * len(order) + 2.2))
im = ax.imshow(M, cmap="RdYlBu_r", vmin=0, vmax=100, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, "%.0f" % M[i, j], ha="center", va="center", fontsize=8,
                color="white" if M[i, j] > 60 or M[i, j] < 12 else "black")
ax.set_xticks(range(len(COLS)))
ax.set_xticklabels([n for _, n in COLS], fontsize=8)
ax.set_yticks(range(len(order)))
ax.set_yticklabels(["%s  (%s)" % (r["family"], MODNAME.get(r["module"], r["module"]))
                    for r in order], fontsize=8)
if ctrl:
    ax.axhline(len(body) - 0.5, color="#333333", lw=1.2)
ax.set_title("Prevalence of chitin-related CAZyme families\nacross the candidate genus and Akkermansia",
             fontsize=10, pad=10)
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
cb.set_label("% genomes carrying the family", fontsize=8)
cb.ax.tick_params(labelsize=8)
for ext in ("pdf", "png"):
    p = "%s/Figure6_chitin_panel.%s" % (OUTD, ext)
    fig.savefig(p, dpi=300, bbox_inches="tight"); print("wrote", p)
