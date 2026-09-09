#!/bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3
import csv, os, statistics as st
B = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
S1 = "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus/tables/TableS1_genome_quality.tsv"
cl = {r["genome"]: r["species_cluster"] for r in csv.DictReader(open(S1), delimiter="\t")}
sz = {r["genome"]: float(r["size_Mb"]) for r in csv.DictReader(open(S1), delimiter="\t")}
cnt = {r["genome"]: (float(r["n_families"]), float(r["n_proteins"])) for r in
       csv.DictReader(open(B + "/results/novel_akk_tree/cazy_per_genome_counts.tsv"), delimiter="\t")
       if r["group"] == "novel"}
miss = [g for g in cnt if g not in cl]
if miss: print("NOT IN TableS1:", miss[:5])
c286 = [g for g in cnt if cl.get(g) == "C286"]
rest = [g for g in cnt if cl.get(g) not in (None, "C286")]
print("C286 n=%d   rest n=%d" % (len(c286), len(rest)))
def rep(lab, gs):
    f = [cnt[g][0] for g in gs]; p = [cnt[g][1] for g in gs]; s = [sz[g] for g in gs]
    fd = [cnt[g][0] / sz[g] for g in gs]; pd = [cnt[g][1] / sz[g] for g in gs]
    print("  %-6s median families %.1f  proteins %.1f  size %.2f Mb  fam/Mb %.1f  prot/Mb %.1f" %
          (lab, st.median(f), st.median(p), st.median(s), st.median(fd), st.median(pd)))
rep("C286", c286); rep("rest", rest)
try:
    from scipy.stats import mannwhitneyu
    for lab, key in (("families per Mb", 0), ("proteins per Mb", 1)):
        a = [cnt[g][key] / sz[g] for g in c286]; b = [cnt[g][key] / sz[g] for g in rest]
        print("  %-18s Mann-Whitney p = %.4g" % (lab, mannwhitneyu(a, b, alternative="two-sided").pvalue))
except ImportError:
    print("  scipy not available, p-values not recomputed")
print("\n=== C286 host species, from TableS3 ===")
os.system("awk -F'\t' 'NR==1 || $9 ~ /EHM/' " + "/bigdata/stajichlab/lshad003/akkermansiaceae-sister-genus/tables/TableS3_host_metadata.tsv | head -3")
print("\n=== files that may hold the 9-overlap null ===")
os.system("ls " + B + "/results/novel_akk_tree/ | grep -i -E 'null|overlap|nest|reduced'")
