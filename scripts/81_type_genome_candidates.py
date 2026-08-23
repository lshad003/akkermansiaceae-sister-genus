#!/usr/bin/env python3
# Type genome candidates ranked against the description criteria
# Source: ch3-chitin-evolution/scripts/type_genome_candidates.py
# Output: results/novel_akk_tree/type_genome_candidates.tsv
import csv, os, statistics, sys

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
MD   = f"{BASE}/results/mimag"
FA   = f"{BASE}/data/amphibia_gtdbtk_input"
CEN  = f"{BASE}/results/gh75_census_v3/gh75_verru_census_per_genome_familyfilled_v3.tsv"
SGC  = f"{BASE}/results/novel_akk_tree/novel_size_gc.tsv"
DER  = f"{BASE}/results/akkfam_derep/clusters.tsv"
OUT  = f"{BASE}/results/novel_akk_tree/type_genome_candidates.tsv"

# ---- parsing identical to scripts/verify_mimag_count2.py
def rrna(acc):
    p = f"{MD}/{acc}.rrna.gff"
    got = set()
    if not os.path.exists(p):
        return got
    for l in open(p, errors="ignore"):
        for k in ("5S", "16S", "23S"):
            if k + "_rRNA" in l or "Name=" + k in l or ("product=" + k) in l:
                got.add(k)
    return got

def trna(acc):
    p = f"{MD}/{acc}.trna.txt"
    got = set()
    if not os.path.exists(p):
        return got
    for l in open(p, errors="ignore"):
        f = l.split()
        if len(f) >= 5 and f[4].isalpha() and len(f[4]) == 3:
            got.add(f[4])
    return got

sz, gc = {}, {}
for r in csv.DictReader(open(SGC), delimiter="\t"):
    sz[r["genome"]] = float(r["genome_size_bp"]) / 1e6
    gc[r["genome"]] = float(r["gc_percent"])

# contamination is NOT in the census. UHM values come from the herptile QC table
# (CheckM v1 lineage_wf), EHI values from the EHI genome list.
HERP = f"{BASE}/data/herptile_cazyme_taxonomy_joined.tsv"
EHIL = f"{BASE}/data/ehi_amphibian_genome_list.tsv"
cont_ext = {}
for path, keycol in ((HERP, "bin_id"), (EHIL, None)):
    if not os.path.exists(path):
        print("MISSING contamination source:", path)
        continue
    rdr = csv.DictReader(open(path), delimiter="\t")
    cols = rdr.fieldnames
    kc = keycol if keycol in (cols or []) else (cols[0] if cols else None)
    cc = next((c for c in (cols or []) if c.lower() == "contamination"), None)
    if kc is None or cc is None:
        print("no usable columns in", path, cols)
        continue
    n = 0
    for r in rdr:
        v = (r.get(cc) or "").strip()
        try:
            cont_ext[r[kc].strip()] = float(v)
            n += 1
        except ValueError:
            pass
    print("contamination read from %s: %d rows" % (os.path.basename(path), n))

comp, cont, dset = {}, {}, {}
for r in csv.DictReader(open(CEN), delimiter="\t"):
    if r["family"] != "Akkermansiaceae":
        continue
    g = (r["genus"] or "").strip()
    if (not g or g in ("unknown", "NO_GENUS")) and r["host_class"] == "amphibian":
        try:
            comp[r["accession"]] = float(r["completeness"])
        except ValueError:
            pass
        cont[r["accession"]] = cont_ext.get(r["accession"], float("nan"))
        dset[r["accession"]] = r["from_dataset"]

cluster = {}
red = set()
for r in csv.DictReader(open(DER), delimiter="\t"):
    if r["genus"] != "NOVEL":
        continue
    for m in r["members"].split(","):
        cluster[m] = r["cluster_id"]
    if int(r["n_genomes"]) == 18:
        red = set(r["members"].split(","))

def contigs(acc):
    for ext in (".fa", ".fna", ".fasta"):
        p = f"{FA}/{acc}{ext}"
        if os.path.exists(p):
            return sum(1 for l in open(p, errors="ignore") if l.startswith(">"))
    return None

genomes = sorted(g for g in comp if g in sz)
print("genomes with completeness and size:", len(genomes))
med_sz = statistics.median([sz[g] for g in genomes])
med_gc = statistics.median([gc[g] for g in genomes])
print("genus median: %.2f Mb, %.2f%% GC" % (med_sz, med_gc))

rows = []
for g in genomes:
    rr = rrna(g); tr = trna(g); nc = contigs(g)
    dist = abs(sz[g] - med_sz) / med_sz + abs(gc[g] - med_gc) / med_gc
    rows.append(dict(genome=g,
                     source="EHI" if g.startswith("EHM") else "UHM",
                     comp=comp[g], cont=cont.get(g, float("nan")),
                     contigs=nc if nc is not None else -1,
                     r5="+" if "5S" in rr else "-",
                     r16="+" if "16S" in rr else "-",
                     r23="+" if "23S" in rr else "-",
                     ntrna=len(tr),
                     cluster=cluster.get(g, "?"),
                     c286="YES" if g in red else "no",
                     size=sz[g], gc=gc[g], dist=dist))

def qualifies(r):
    return (r["comp"] >= 90 and r["cont"] == r["cont"] and r["cont"] < 5
            and r["c286"] == "no" and r["r5"] == "+" and r["r16"] == "+"
            and r["r23"] == "+" and r["ntrna"] >= 18)

hdr = "%-26s %-5s %7s %6s %7s %-6s %5s %-6s %-5s %6s %6s %7s"
print()
print("=== EHM058340, the proposed type ===")
print(hdr % ("genome","src","comp","cont","contigs","5/16/23","tRNA","clust","C286","Mb","GC","medDist"))
for r in rows:
    if r["genome"] == "EHM058340":
        print(hdr % (r["genome"], r["source"], "%.2f" % r["comp"], "%.2f" % r["cont"],
                     r["contigs"], "%s%s%s" % (r["r5"], r["r16"], r["r23"]), r["ntrna"],
                     r["cluster"], r["c286"], "%.2f" % r["size"], "%.2f" % r["gc"],
                     "%.3f" % r["dist"]))

uhm = [r for r in rows if r["source"] == "UHM" and qualifies(r)]
uhm.sort(key=lambda r: (-r["comp"], r["cont"], r["contigs"], r["dist"]))
print()
print("=== UHM candidates meeting every criterion, ranked ===")
print("criteria: completeness >= 90, contamination < 5, outside C286, 5S+16S+23S present, tRNA types >= 18")
print("ranked by completeness, then contamination, then contigs, then distance to the genus median")
print(hdr % ("genome","src","comp","cont","contigs","5/16/23","tRNA","clust","C286","Mb","GC","medDist"))
for r in uhm[:12]:
    print(hdr % (r["genome"], r["source"], "%.2f" % r["comp"], "%.2f" % r["cont"],
                 r["contigs"], "%s%s%s" % (r["r5"], r["r16"], r["r23"]), r["ntrna"],
                 r["cluster"], r["c286"], "%.2f" % r["size"], "%.2f" % r["gc"],
                 "%.3f" % r["dist"]))
print()
print("UHM genomes meeting every criterion: %d of %d UHM genomes"
      % (len(uhm), sum(1 for r in rows if r["source"] == "UHM")))

with open(OUT, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["genome","source","completeness","contamination","contigs",
                "rRNA_5S","rRNA_16S","rRNA_23S","tRNA_types","species_cluster",
                "in_reduced_C286","size_Mb","gc_percent","dist_to_genus_median",
                "meets_all_criteria"])
    for r in sorted(rows, key=lambda x: (-x["comp"], x["cont"])):
        w.writerow([r["genome"], r["source"], "%.2f" % r["comp"], "%.2f" % r["cont"],
                    r["contigs"], r["r5"], r["r16"], r["r23"], r["ntrna"], r["cluster"],
                    r["c286"], "%.2f" % r["size"], "%.2f" % r["gc"], "%.3f" % r["dist"],
                    "yes" if qualifies(r) else "no"])
print("wrote", OUT)
