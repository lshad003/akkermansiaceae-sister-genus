#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
# Per species cluster: n, completeness, raw size, completeness-corrected size, GC. From TableS1.
import statistics as st
S1 = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus/tables/TableS1_genome_quality.tsv"
rows = []
with open(S1) as f:
    h = f.readline().rstrip("\n").split("\t")
    ix = {k: h.index(k) for k in ("genome", "completeness", "species_cluster", "size_Mb", "gc_percent")}
    for line in f:
        x = line.rstrip("\n").split("\t")
        rows.append((x[ix["genome"]], float(x[ix["completeness"]]), x[ix["species_cluster"]],
                     float(x[ix["size_Mb"]]), float(x[ix["gc_percent"]])))
print("genomes:", len(rows), " clusters:", len({r[2] for r in rows}))
print("%-8s %3s %8s %8s %8s %8s %8s" % ("cluster", "n", "min_cmp", "med_cmp", "med_Mb", "corr_Mb", "med_GC"))
for c in sorted({r[2] for r in rows}):
    g = [r for r in rows if r[2] == c]
    print("%-8s %3d %8.1f %8.1f %8.2f %8.2f %8.2f" % (c, len(g), min(r[1] for r in g), st.median(r[1] for r in g),
          st.median(r[3] for r in g), st.median(r[3] / (r[1] / 100) for r in g), st.median(r[4] for r in g)))
print("\nsingletons and their completeness:")
for c in sorted({r[2] for r in rows}):
    g = [r for r in rows if r[2] == c]
    if len(g) == 1: print("  %s %s %.1f" % (c, g[0][0], g[0][1]))
hi = [r for r in rows if r[1] >= 90]
print("\nall 105: size median %.2f (%.2f to %.2f), GC median %.2f (%.2f to %.2f)" % (
    st.median(r[3] for r in rows), min(r[3] for r in rows), max(r[3] for r in rows),
    st.median(r[4] for r in rows), min(r[4] for r in rows), max(r[4] for r in rows)))
print(">=90 complete, n=%d: size median %.2f (%.2f to %.2f)" % (len(hi),
    st.median(r[3] for r in hi), min(r[3] for r in hi), max(r[3] for r in hi)))
lo = min(rows, key=lambda r: r[3]); print("smallest genome: %s %.2f Mb at %.1f complete" % (lo[0], lo[3], lo[1]))
c286 = [r for r in rows if r[2] == "C286"]
print("C286 n=%d: raw median %.2f Mb, corrected %.2f Mb, GC %.2f, min completeness %.1f" % (len(c286),
    st.median(r[3] for r in c286), st.median(r[3] / (r[1] / 100) for r in c286), st.median(r[4] for r in c286), min(r[1] for r in c286)))
