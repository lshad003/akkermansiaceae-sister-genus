import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
SRC  = BASE + "/results/novel_akk_tree/novel_genus_function.tsv"
FREE = BASE + "/results/novel_akk_tree/cazy_polarized.tsv"
OUT  = BASE + "/results/novel_akk_tree/Figure4_mucin_prevalence"

rows = {r["family"]: r for r in csv.DictReader(open(SRC), delimiter="\t")}
fr   = {r["family"]: r for r in csv.DictReader(open(FREE), delimiter="\t")}

mucin = [f for f in rows if rows[f].get("in_mucin_panel") == "1"]
mucin.sort(key=lambda f: float(rows[f]["pct_novel"]), reverse=True)
miss = [f for f in mucin if f not in fr]
if miss:
    print("REFUSED: no pct_FREE for", miss); raise SystemExit(1)
print("families:", len(mucin))

cols = [("pct_novel",         "Ca. Novel genus\namphibian (105)"),
        ("pct_akk_amph",      "Akkermansia\namphibian (94)"),
        ("pct_akk_podarcis",  "Akkermansia\nreptile (137)"),
        ("pct_akk_mammal",    "Akkermansia\nmammal (71)"),
        ("pct_akk_gtdb",      "Akkermansia\nGTDB (42)")]

M = np.array([[float(rows[f][c]) for c, _ in cols] + [float(fr[f]["pct_FREE"])]
              for f in mucin])
labels = [l for _, l in cols] + ["Free-living\nAkkermansiaceae (280)"]

fig, ax = plt.subplots(figsize=(10.5, 8))
im = ax.imshow(M, cmap="RdYlBu_r", vmin=0, vmax=100, aspect="auto")
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=9)
ax.set_yticks(range(len(mucin))); ax.set_yticklabels(mucin, fontsize=9)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        ax.text(j, i, "%.0f" % v, ha="center", va="center", fontsize=8,
                color="white" if (v > 75 or v < 25) else "black")
ax.axvline(4.5, color="black", lw=1.6)
ax.set_title("Prevalence of mucin-associated CAZyme families", fontsize=11)
cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("% genomes", fontsize=9)
plt.tight_layout()
plt.savefig(OUT + ".pdf", bbox_inches="tight")
plt.savefig(OUT + ".png", dpi=300, bbox_inches="tight")

print()
for f in mucin:
    print("  %-8s novel %5.1f   free %5.1f" % (f, float(rows[f]["pct_novel"]), float(fr[f]["pct_FREE"])))
print()
print("WROTE", OUT + ".pdf and .png")
