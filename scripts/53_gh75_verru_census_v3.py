#!/usr/bin/env python3
# Family prevalence census under the uniform filter
# Source: ch3-chitin-evolution/scripts/gh75_verru_census_v3.py
# Output: results/gh75_census_v3/gh75_verru_census_per_genome_v3.tsv
"""
STEP 5: GH75 prevalence census across ALL annotated Verrucomicrobiota.

Run:
  /bigdata/stajichlab/lshad003/condaenvs/rf_py39/bin/python3 \
      /bigdata/stajichlab/lshad003/ch3-chitin-evolution/scripts/gh75_verru_census.py

Runtime ~2-4 min (IO bound: streams ~1.5 GB of dbCAN tsv, only parses GH75 lines).

WHAT THIS ANSWERS
  Is GH75 amphibian/gut-concentrated, or spread across Verrucomicrobiota?

THREE THINGS THIS SCRIPT REFUSES TO DO
  1. Refuses to report 0% GH75 for a dataset that was never annotated. Feng (36 Verru)
     and Teullet (3 Verru) have NO dbCAN output anywhere; 232 of 368 EHI Verru are
     non-amphibian and were never annotated. Absent != negative. Those are reported in
     a COVERAGE block and excluded from every prevalence denominator.
  2. Refuses to call GH75 gut-concentrated without testing whether the gut effect is
     just clade. Akkermansiaceae is gut-restricted by biology. If GH75 is an
     Akkermansiaceae trait, it will LOOK gut-enriched even with zero ecological
     signal. The within-clade / outside-clade contrast is the actual test.
  3. Refuses to report prevalence without the completeness + MAG-vs-isolate strata.
     GH75 absence in a 62%-complete MAG is not evidence of absence. If prevalence
     rises with completeness, every raw prevalence below is a floor, not an estimate.

SOURCES OF GH75 CALLS (method-matched: all dbCAN db-cazy/11.0 + hmmer3.3.2 + hmmscan-parser)
  - results/dbcan_verru/          2,965 (this array)
  - results/dbcan_{bact_refs,ehi_amphibian,endo_allphyla,flavo_refs,scaffold}/  prior, same chain
  - data/herptile_cazyme_taxonomy_joined.tsv   147 herptile Verru (precomputed GH75_count;
    no per-genome tsv on disk, this table is the project standard used for Q1/Q3)

COMPLETENESS
  catalog `completeness` is EMPTY for GTDB_r226_rep (populated only for the 554 MAG-dataset
  rows). Joined from data/gtdb_r226/bac120_metadata_r226.tsv (checkm2_completeness).
"""
import csv, os, re, sys, math
from collections import Counter, defaultdict

BASE = "/bigdata/stajichlab/lshad003/ch3-chitin-evolution"
CAT  = f"{BASE}/data/unified_annotation/unified_genome_annotation.tsv"
HERP = f"{BASE}/data/herptile_cazyme_taxonomy_joined.tsv"
GTDBMETA = f"{BASE}/data/gtdb_r226/bac120_metadata_r226.tsv"
OUTDIR = f"{BASE}/results/gh75_census_v3"
os.makedirs(OUTDIR, exist_ok=True)
DBCAN_DIRS = ["dbcan_verru", "dbcan_bact_refs", "dbcan_ehi_amphibian",
              "dbcan_ehi_nonamph",  # ADDED 2026-07-13: the 232 reptile+mammal EHI Verru (job 26223734)
              "dbcan_endo_allphyla", "dbcan_flavo_refs", "dbcan_scaffold",
              # ADDED 2026-07-15 (v3): raw dbCAN for the 147 herptile_MAG Verru. Chain verified
              # byte-for-byte vs cazy_ehi_amphibian_array.sh: prodigal -p meta, HMMER 3.3.2,
              # CAZyDB v11.0, hmmscan-parser.sh. Replaces the precomputed table branch.
              "/bigdata/stajichlab/lshad003/ruminococcaceae-agent/results/dbcan_allphyla"]

# v3: the filter, applied UNIFORMLY to every genome. v2 applied NO filter to tsv-read genomes
# but inherited E<=1e-10 AND cov>=0.30 for herptile rows via the precomputed table, so the two
# halves of the amphibian set were filtered differently. That is the bug this fixes.
EVAL_MAX = 1e-10
COV_MIN  = 0.30

# exact family token. Must NOT let GH7 or GH750 match, and must allow a GH75_1 subfamily.
GH75_RE = re.compile(r"^GH75(_\d+)?\.hmm$")

os.makedirs(OUTDIR, exist_ok=True)
csv.field_size_limit(10_000_000)

def bare(a):
    return a[3:] if a.startswith(("GB_", "RS_")) else a

def variants(a):
    b = bare(a)
    return [a, b, "GB_" + b, "RS_" + b]

# ----------------------------------------------------------------------------------
# body-site keyword remap. Applied ONLY where catalog body_site == 'unknown'.
# Never overwrites a body_site the catalog already assigned.
# Short/risky tokens use word boundaries: \bgut\b not "gut" (else "gutter"),
# \bcolon\b not "colon" (else "colony"), \bsand\b (else "thousand"), \bice\b (else "juice"),
# \bair\b (else "hair"). "digestive tract" not "digest" (else "anaerobic digestion of waste").
# ----------------------------------------------------------------------------------
GUT_RE = re.compile(
    r"\bgut\b|hindgut|midgut|foregut|intestin|\bfaec|\bfec(es|al)\b|faeces|stool|"
    r"\brumen\b|ruminant|gastrointestinal|gastro-intestinal|c(a)?ecum|c(a)?ecal|"
    r"\bcolon\b|colonic|digestive tract|digestive system|gizzard|\bcrop\b|"
    r"\bdung\b|coprolite|hindgut|\bileum\b|jejunum|duodenum|\bcloaca")
NONGUT_RE = re.compile(
    r"water|seawater|marine|freshwater|\blake\b|\bpond\b|river|ocean|sediment|\bsoil\b|"
    r"\brock\b|groundwater|aquifer|wastewater|sludge|bioreactor|digester|anaerobic digestion|"
    r"hot spring|hydrothermal|cold seep|rhizosphere|\broot\b|phyllosphere|\bleaf\b|\bplant\b|"
    r"coral|algae|macroalgea|macroalgae|seaweed|kelp|wetland|\bpeat\b|\bbog\b|permafrost|"
    r"glacier|\bice\b|\bsand\b|biofilm|microbial mat|biofloc|aquatic|aquaculture|compost|"
    r"landfill|mine drainage|\bcave\b|\bair\b|\bdust\b|sponge|tidal|estuar|saline|brine|"
    r"\bmoss\b|lichen|silage|fermented|\bskin\b|\boral\b|\bmouth\b|saliva|dental|plaque|"
    r"vagina|nasal|\blung\b|respiratory|\bblood\b|urin|\bmilk\b|udder|\bgill\b|\bteat\b")
# deliberately ambiguous -> stay 'unknown' rather than be forced into a bin.
AMBIG_RE = re.compile(r"manure|slurry|\bmucus\b|\bmucosa\b|biosolid")

def remap_body_site(iso):
    s = (iso or "").strip().lower()
    if s in ("", "none", "na", "n/a", "missing", "not applicable", "unknown", "metagenome"):
        return "unknown", "no_isolation_source"
    if AMBIG_RE.search(s):
        return "unknown", "ambiguous_term"
    if GUT_RE.search(s):
        return "gut", "remap_gut"
    if NONGUT_RE.search(s):
        return "non_gut", "remap_non_gut"
    return "unknown", "no_keyword_match"

def fisher_two_sided(a, b, c, d):
    """2x2 exact test. [[a,b],[c,d]]. Dependency-free."""
    n = a + b + c + d
    if n == 0 or (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
        return float("nan")
    def p_of(x):
        return (math.comb(a + b, x) * math.comb(c + d, a + c - x)) / math.comb(n, a + c)
    lo = max(0, a + c - (c + d))
    hi = min(a + b, a + c)
    p_obs = p_of(a)
    tot = 0.0
    for x in range(lo, hi + 1):
        p = p_of(x)
        if p <= p_obs * (1 + 1e-9):
            tot += p
    return min(1.0, tot)

def pct(n, d):
    return "  n/a" if d == 0 else f"{100.0*n/d:5.1f}"

# ----------------------------------------------------------------------------------
# 1. catalog: all Verrucomicrobiota
# ----------------------------------------------------------------------------------
verru = [r for r in csv.DictReader(open(CAT), delimiter="\t")
         if "Verrucomicrobiota" in r["phylum"]]
print("=" * 84)
print("GH75 PREVALENCE CENSUS -- VERRUCOMICROBIOTA")
print("=" * 84)
print(f"Verrucomicrobiota in unified catalog: {len(verru)}")

# ----------------------------------------------------------------------------------
# 2. index every per-genome dbCAN tsv (same annotation chain)
# ----------------------------------------------------------------------------------
idx = {}
for d in DBCAN_DIRS:
    p = d if d.startswith("/") else f"{BASE}/results/{d}"
    if not os.path.isdir(p):
        continue
    for fn in os.listdir(p):
        # guard: .cazyme.tsv is dbcanlight output, a DIFFERENT tool. Never mix chains.
        if fn.endswith(".tsv") and not fn.endswith(".cazyme.tsv"):
            idx.setdefault(fn[:-4], f"{p}/{fn}")

# herptile GH75 calls come from the precomputed join table (no per-genome tsv exists)
herp_gh75 = {}
for r in csv.DictReader(open(HERP), delimiter="\t"):
    v = (r.get("GH75_count") or "").strip()
    try:
        n = int(float(v)) if v not in ("", "NA") else 0
    except ValueError:
        n = 1 if str(r.get("GH75_present", "")).lower() in ("1", "true", "yes") else 0
    herp_gh75[r["bin_id"]] = n

# ----------------------------------------------------------------------------------
# 3. GTDB CheckM completeness (catalog completeness is EMPTY for GTDB reps)
# ----------------------------------------------------------------------------------
gtdb_comp, gtdb_cat = {}, {}
for r in csv.DictReader(open(GTDBMETA), delimiter="\t"):
    a = r["accession"]
    c = (r.get("checkm2_completeness") or r.get("checkm_completeness") or "").strip()
    if c not in ("", "NA", "none"):
        try:
            gtdb_comp[bare(a)] = float(c)
        except ValueError:
            pass
    gtdb_cat[bare(a)] = (r.get("ncbi_genome_category") or "").strip()

# ----------------------------------------------------------------------------------
# 4. per-genome GH75 call
# ----------------------------------------------------------------------------------
rows = []
unannotated = Counter()
scanned = 0
for r in verru:
    acc = r["accession"]
    ds = r["from_dataset"]
    n_dom, n_prot, annotated = 0, 0, False

    tsv = next((idx[v] for v in variants(acc) if v in idx), None)
    if tsv:
        annotated = True
        prots = set()
        with open(tsv) as fh:
            for line in fh:
                if not line.startswith("GH75"):   # cheap prefilter before regex
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10 or not GH75_RE.match(f[0]):
                    continue
                # v3: uniform quality filter. hmmscan-parser.sh columns are
                # 0 FAM.hmm  1 hmm_len  2 prot_id  3 prot_len  4 evalue  ...  9 coverage
                try:
                    ev = float(f[4]); cv = float(f[9])
                except ValueError:
                    continue
                if ev > EVAL_MAX or cv < COV_MIN:
                    continue
                n_dom += 1
                prots.add(f[2])
        n_prot = len(prots)
        scanned += 1
    elif acc in herp_gh75:
        annotated = True
        n_dom = herp_gh75[acc]
        n_prot = n_dom
    else:
        unannotated[ds] += 1

    # completeness: catalog for MAG datasets, GTDB CheckM2 for reps
    comp = ""
    cv = (r.get("completeness") or "").strip()
    if cv not in ("", "NA"):
        comp = cv
    elif bare(acc) in gtdb_comp:
        comp = f"{gtdb_comp[bare(acc)]:.2f}"

    # genome type: MAG vs isolate vs SAG
    gc = (r.get("ncbi_genome_category") or "").strip().lower()
    if not gc or gc == "empty":
        gc = gtdb_cat.get(bare(acc), "").lower()
    if "metagenome" in gc:
        gtype = "MAG"
    elif "single cell" in gc:
        gtype = "SAG"
    elif gc in ("none", ""):
        # 'none' in NCBI = not metagenome-derived = cultured isolate.
        # herptile/EHI/Feng/Teullet MAGs have an empty field but ARE MAGs by construction.
        gtype = "MAG" if ds != "GTDB_r226_rep" else "isolate"
    else:
        gtype = "isolate"

    bs_raw = (r.get("body_site") or "unknown").strip() or "unknown"
    if bs_raw == "unknown":
        bs_final, bs_src = remap_body_site(r.get("ncbi_isolation_source"))
    else:
        bs_final, bs_src = bs_raw, "catalog"

    rows.append(dict(
        accession=acc, from_dataset=ds, annotated=int(annotated),
        gh75_domains=n_dom, gh75_proteins=n_prot, gh75_present=int(n_dom > 0),
        family=r.get("family", ""), genus=r.get("genus", ""),
        body_site_catalog=bs_raw, body_site_final=bs_final, body_site_source=bs_src,
        isolation_source=(r.get("ncbi_isolation_source") or "")[:120],
        genome_type=gtype, completeness=comp,
        host_class=r.get("host_class", ""), host_animal_type=r.get("host_animal_type", ""),
        captivity_status=r.get("captivity_status", ""), sample_id=r.get("sample_id", ""),
    ))

out_tsv = f"{OUTDIR}/gh75_verru_census_per_genome_v2.tsv"
with open(out_tsv, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

ann = [r for r in rows if r["annotated"]]

# ----------------------------------------------------------------------------------
# COVERAGE (read this before any prevalence number)
# ----------------------------------------------------------------------------------
print("\n" + "=" * 84)
print("COVERAGE -- what was actually annotated")
print("=" * 84)
print(f"{'from_dataset':<22}{'catalog':>9}{'annotated':>11}{'NO ANNOT':>10}")
for ds in sorted({r["from_dataset"] for r in rows}):
    tot = sum(1 for r in rows if r["from_dataset"] == ds)
    a = sum(1 for r in rows if r["from_dataset"] == ds and r["annotated"])
    print(f"{ds:<22}{tot:>9}{a:>11}{tot-a:>10}")
print(f"{'TOTAL':<22}{len(rows):>9}{len(ann):>11}{len(rows)-len(ann):>10}")
if unannotated:
    print("\n  WARNING: the following have NO dbCAN output anywhere. They are NOT GH75-negative,")
    print("  they were never looked at. Excluded from every denominator below.")
    for ds, n in unannotated.most_common():
        print(f"    {ds:<22} {n:>5} unannotated")

# ----------------------------------------------------------------------------------
# HEADLINE
# ----------------------------------------------------------------------------------
pos = sum(r["gh75_present"] for r in ann)
print("\n" + "=" * 84)
print("HEADLINE")
print("=" * 84)
print(f"Annotated Verrucomicrobiota:        {len(ann)}")
print(f"Carrying >= 1 GH75 domain:          {pos}  ({pct(pos,len(ann))}%)")
tot_dom = sum(r["gh75_domains"] for r in ann)
print(f"Total GH75 domains:                 {tot_dom}")
if pos:
    print(f"Mean GH75 domains per carrier:      {tot_dom/pos:.2f}")
    print(f"Max GH75 domains in one genome:     {max(r['gh75_domains'] for r in ann)}")

def block(title, keyfn, subset=None, minn=1):
    data = subset if subset is not None else ann
    print("\n" + "-" * 84)
    print(title)
    print("-" * 84)
    print(f"{'group':<34}{'n':>7}{'GH75+':>7}{'%':>8}")
    g = defaultdict(list)
    for r in data:
        g[keyfn(r)].append(r)
    for k in sorted(g, key=lambda k: -len(g[k])):
        v = g[k]
        if len(v) < minn:
            continue
        p = sum(x["gh75_present"] for x in v)
        print(f"{str(k):<34}{len(v):>7}{p:>7}{pct(p,len(v)):>8}")

# ----------------------------------------------------------------------------------
# BY FROM_DATASET
# ----------------------------------------------------------------------------------
block("GH75 PREVALENCE BY from_dataset (annotated only)", lambda r: r["from_dataset"])
print("  note: Feng_2021_chicken and Teullet_2023 are ABSENT above because 0 of their")
print("        Verru were annotated. EHI row = amphibian subset only (136 of 368).")

# ----------------------------------------------------------------------------------
# BODY SITE + REMAP AUDIT
# ----------------------------------------------------------------------------------
block("GH75 PREVALENCE BY body_site (catalog + isolation-source keyword remap)",
      lambda r: r["body_site_final"])
print("\n  remap audit -- where each body_site call came from:")
for src, n in Counter(r["body_site_source"] for r in ann).most_common():
    print(f"    {src:<24} {n:>6}")
print("\n  strings the remap RESCUED from 'unknown' (verify these; a bad keyword here")
print("  moves the headline). Top 12 per class:")
for cls in ("gut", "non_gut"):
    ex = Counter(r["isolation_source"].lower()[:52] for r in ann
                 if r["body_site_source"] == f"remap_{cls}")
    print(f"    --> {cls} ({sum(ex.values())} genomes rescued)")
    for s, n in ex.most_common(12):
        print(f"          {n:>4}  {s}")

# ----------------------------------------------------------------------------------
# MAG vs ISOLATE  (caveat #1)
# ----------------------------------------------------------------------------------
block("GH75 PREVALENCE BY genome_type -- THE MAG-vs-ISOLATE CAVEAT",
      lambda r: r["genome_type"])
mag = [r for r in ann if r["genome_type"] == "MAG"]
iso = [r for r in ann if r["genome_type"] == "isolate"]
if mag and iso:
    pm, pi = sum(r["gh75_present"] for r in mag), sum(r["gh75_present"] for r in iso)
    p = fisher_two_sided(pi, len(iso) - pi, pm, len(mag) - pm)
    print(f"\n  isolates {pct(pi,len(iso))}%  vs  MAGs {pct(pm,len(mag))}%   Fisher p = {p:.3g}")
    print("  If isolates are HIGHER, MAG prevalence is a FLOOR: fragmentary assemblies drop")
    print("  real GH75 genes, and every % above understates the truth. Verrucomicrobiota is")
    print("  overwhelmingly MAG-based, so this caps how hard any absence can be pushed.")

# ----------------------------------------------------------------------------------
# COMPLETENESS STRATIFICATION  (caveat #2)
# ----------------------------------------------------------------------------------
print("\n" + "-" * 84)
print("GH75 PREVALENCE BY COMPLETENESS -- IS ABSENCE REAL, OR JUST MISSING SEQUENCE?")
print("-" * 84)
hasc = [r for r in ann if r["completeness"] not in ("", "NA")]
print(f"  completeness available for {len(hasc)} / {len(ann)} annotated genomes")
if len(hasc) < len(ann):
    miss = Counter(r["from_dataset"] for r in ann if r["completeness"] in ("", "NA"))
    print(f"  missing completeness: {dict(miss)}")
BINS = [(0, 70), (70, 80), (80, 90), (90, 95), (95, 100.01)]
print(f"\n{'completeness bin':<34}{'n':>7}{'GH75+':>7}{'%':>8}")
for lo, hi in BINS:
    v = [r for r in hasc if lo <= float(r["completeness"]) < hi]
    if not v:
        continue
    p = sum(x["gh75_present"] for x in v)
    print(f"{f'{lo}-{hi if hi<=100 else 100}%':<34}{len(v):>7}{p:>7}{pct(p,len(v)):>8}")
hi_only = [r for r in hasc if float(r["completeness"]) >= 90]
if hi_only:
    p = sum(r["gh75_present"] for r in hi_only)
    print(f"\n  RESTRICTED TO >=90% COMPLETE: {p}/{len(hi_only)} = {pct(p,len(hi_only))}% carry GH75")
    print("  Compare to the headline. If it moves, the headline was completeness-limited.")

# ----------------------------------------------------------------------------------
# CLADE -- is GH75 a phylum-wide trait or one family's trait?
# ----------------------------------------------------------------------------------
block("GH75 PREVALENCE BY FAMILY (families with >= 20 annotated genomes)",
      lambda r: r["family"] or "unclassified", minn=20)
akk = [r for r in ann if "Akkermansia" in (r["family"] or "")]
non = [r for r in ann if "Akkermansia" not in (r["family"] or "")]
pa, pn = sum(r["gh75_present"] for r in akk), sum(r["gh75_present"] for r in non)
print(f"\n  Akkermansiaceae:      {pa}/{len(akk)} = {pct(pa,len(akk))}%")
print(f"  all other families:   {pn}/{len(non)} = {pct(pn,len(non))}%")
print(f"  share of ALL GH75+ genomes that are Akkermansiaceae: {pct(pa, pos)}%")
print("  The May analysis claimed GH75 is Akkermansia-specific. If 'all other families'")
print("  is well above zero, that claim is DEAD and GH75 is a phylum-wide trait.")

# ----------------------------------------------------------------------------------
# THE ACTUAL TEST: is the gut signal independent of clade?
# ----------------------------------------------------------------------------------
print("\n" + "=" * 84)
print("IS GUT-ENRICHMENT REAL, OR IS IT JUST AKKERMANSIACEAE PHYLOGENY?")
print("=" * 84)
print("Akkermansiaceae is gut-restricted by biology. If GH75 is an Akkermansiaceae trait,")
print("it will LOOK gut-enriched with zero ecological signal. So: hold clade fixed.")
print("The GTDB-rep rows are species-dereplicated, so a genome-level test is defensible")
print("HERE (unlike herptile MAGs, where replication is host-animal-level, not MAG-level).")
for label, subset in (("ALL annotated", ann),
                      ("OUTSIDE Akkermansiaceae", non),
                      ("WITHIN Akkermansiaceae", akk)):
    g = [r for r in subset if r["body_site_final"] == "gut"]
    ng = [r for r in subset if r["body_site_final"] == "non_gut"]
    if not g or not ng:
        print(f"\n  {label:<26} insufficient gut/non_gut split (gut={len(g)} non_gut={len(ng)})")
        continue
    pg, pnn = sum(r["gh75_present"] for r in g), sum(r["gh75_present"] for r in ng)
    p = fisher_two_sided(pg, len(g) - pg, pnn, len(ng) - pnn)
    print(f"\n  {label}")
    print(f"    gut     {pg:>5}/{len(g):<5} = {pct(pg,len(g))}%")
    print(f"    non_gut {pnn:>5}/{len(ng):<5} = {pct(pnn,len(ng))}%")
    print(f"    Fisher p = {p:.3g}")
print("\n  READ IT THIS WAY: if gut-enrichment VANISHES outside Akkermansiaceae, then GH75")
print("  is clade-concentrated, not gut-concentrated, and 'gut' was a proxy for Akkermansia.")
print("  If it SURVIVES outside Akkermansiaceae, the ecological signal is real.")

# ----------------------------------------------------------------------------------
# AMPHIBIAN SLICE (with the replication + captivity traps stated)
# ----------------------------------------------------------------------------------
print("\n" + "=" * 84)
print("HOST-ASSOCIATED (herptile + EHI) SLICE")
print("=" * 84)
# 'herptile' = amphibians AND reptiles. Do NOT call this slice 'amphibian' without
# splitting host_class, or reptile MAGs get smuggled into an amphibian number.
amph = [r for r in ann if r["from_dataset"] in ("herptile_MAG", "EHI_2025")]
if amph:
    pa2 = sum(r["gh75_present"] for r in amph)
    print(f"  herptile+EHI Verru (annotated): {len(amph)}, GH75+ {pa2} = {pct(pa2,len(amph))}%")
    # corrected animal unit (sample_id is the sequencing run for herptile_MAG; see CHATINDEX
    # "THE UNIT OF REPLICATION"). herptile_MAG -> accession prefix; else -> sample_id.
    def _animal(r):
        return r["accession"].split(".")[0] if r["from_dataset"] == "herptile_MAG" else (r["sample_id"] or "")
    ns = len({_animal(r) for r in amph if _animal(r)})
    print(f"  distinct host animals behind them: {ns}")
    print(f"\n{'host_class (Amphibia vs Reptilia)':<34}{'n':>7}{'GH75+':>7}{'%':>8}")
    for k in sorted({(r["host_class"] or "unknown") for r in amph}):
        v = [r for r in amph if (r["host_class"] or "unknown") == k]
        p = sum(x["gh75_present"] for x in v)
        print(f"{k:<34}{len(v):>7}{p:>7}{pct(p,len(v)):>8}")
    print("  The AMPHIBIA row is the only one that speaks to the amphibian claim.")
    print("  TRAP: these MAGs are NOT independent. Effective replication is host animals")
    print("  (~28 in wild herptile), not MAG count. Do not run a MAG-level test on this")
    print("  slice and call it significant. It is descriptive.")
    print(f"\n{'captivity_status':<34}{'n':>7}{'GH75+':>7}{'%':>8}")
    for k, v in sorted(defaultdict(list, {k: [r for r in amph if (r['captivity_status'] or 'unknown') == k]
                       for k in {(r['captivity_status'] or 'unknown') for r in amph}}).items()):
        p = sum(x["gh75_present"] for x in v)
        print(f"{k:<34}{len(v):>7}{p:>7}{pct(p,len(v)):>8}")
    print("  TRAP: captivity confound. A captive-only signal is a bloom, not a diet trait.")
    gen = Counter(r["genus"] for r in amph if r["gh75_present"])
    print("\n  genus of GH75+ amphibian Verru (is the amphibian signal one clade repeated?):")
    for k, n in gen.most_common(10):
        print(f"    {n:>4}  {k or 'unclassified'}")

print("\n" + "=" * 84)
print(f"per-genome table -> {out_tsv}")
print("=" * 84)
