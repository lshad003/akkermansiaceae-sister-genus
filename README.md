# An unnamed Akkermansiaceae genus sister to Akkermansia: analysis record

This repository documents the analyses used to identify a previously unnamed
Akkermansiaceae genus recovered from wild amphibian gut metagenomes, establish
its position as the sister group of Akkermansia, describe its genomic
properties, and determine what its gene content implies about the origin of
Akkermansia. The gene-content work is split between a repertoire that is
conserved across the gut-associated clade and one that diverges sharply
between the candidate genus and Akkermansia.

Each step below summarizes the analysis, its main result, and the scripts used.

## Steps

### Step 1. Identification of the candidate genome set

Akkermansiaceae genomes carrying no genus assignment are pulled from the
genome catalogue and checked against the classifier output to confirm that the
missing assignment reflects taxonomic novelty rather than a failed call.
Because the catalogue sample identifier is a sequencing run rather than an
animal, host animals are recounted under a corrected unit before any
host-level number is quoted, and captivity status is reconciled across
metadata sources, since a captivity confound would invalidate later
comparisons.

**Main result.** Akkermansiaceae genomes carrying no GTDB genus assignment
number 107 in the combined catalogue. Two of these, both recovered from a
single tortoise, share 53.3% amino acid identity with the type genome (53.26
and 53.28% over 1,193 and 1,037 reciprocal best hits) against 60.3 to 99.4%
for the remaining genomes, a seven-point gap with nothing in between, and they
were excluded on that basis rather than by a decision taken in advance. The
candidate set is therefore 105 MAGs, split 54 from the herptile MAG catalogue
and 51 from the EHI catalogue. All 105 are from amphibians; this reflects the
exclusion together with the amphibian-weighted composition of the source
catalogues, so it describes the sampled range rather than establishing host
specificity. All 105 matched the annotation table with none unmatched, and
median completeness is 96.6% (range 51.2 to 100.0). Under the corrected animal
unit the set represents 66 host animals rather than the larger count implied
by sequencing-run identifiers, an oversplit that inflates the amphibian side
specifically. All 105 genomes are from wild animals. Chimerism screening
against proGenomes 2.1 passes 102 of the 105 at the default threshold, against
86 of 94 amphibian Akkermansia screened in the same run, and returns no genome
meeting a strict chimera call in either group. All 105 candidate genomes fall
below a reference representation score of 0.5, which reflects the absence of a
named representative of this lineage from the reference set rather than
evidence of contamination.

| File | Purpose |
|---|---|
| `scripts/01_check_akk_genus.py` | Genus composition of the Akkermansiaceae reference set |
| `scripts/02_resolve_akk_genus.py` | Genus calls resolved for amphibian Akkermansiaceae |
| `scripts/03_find_51_genus.py` | Genus field checked for the EHI subset |
| `scripts/04_check_51_gtdb_string.py` | GTDB taxonomy strings checked for the EHI subset |
| `scripts/05_recount_host_animals_corrected_unit.py` | Host animals recounted under the corrected animal unit |
| `scripts/06_verify_novel_qc.py` | Completeness and contamination verified for the candidate set |
| `scripts/07_reconcile_captivity.py` | Wild and captive status reconciled across metadata sources |
| `scripts/64_tortoise_exclusion_aai.py` | Non-amphibian genomes excluded by amino acid identity to the type genome |
| `jobs/70_run_gunc_199.sh` | Chimerism screening submitted for the candidate and sister genomes |
| `scripts/71_gunc_audit_199.py` | Chimerism screening results audited by arm |

### Step 2. Phylogenetic placement

Placement is inferred from the concatenated bac120 marker alignment, extracted
from the classifier output in the form that retains reference genomes, since
the user-only alignment cannot place anything. Three questions are asked
directly of the tree: whether the candidate genomes form one clade, whether
any of them fall inside Akkermansia, and what the sister group of the
candidate clade is. The tree is then re-inferred by maximum likelihood with
branch support, because the first inference carries approximate support only.

**Main result.** In a 784-tip alignment the 105 candidate genomes form a
single clade of exactly 105 tips with zero intruders, none of them falls
inside the 187-tip Akkermansia clade, and the parent of the candidate clade
subtends exactly 292 tips, comprising the 105 candidates and the 187
Akkermansia genomes and nothing else. The sister group of the candidate clade
is therefore exactly Akkermansia. Two independent inference methods return
identical topology at all three nodes: approximate likelihood support of 1.0
throughout, and maximum-likelihood support of SH-aLRT 100 with 1000 ultrafast
bootstrap replicates at 100 for the candidate clade, the Akkermansia clade,
and their common parent.

| File | Purpose |
|---|---|
| `scripts/08_check_bac120_msa.py` | Reference-containing bac120 alignment confirmed as the placement input |
| `scripts/09_build_akk_placement_tree.py` | Placement alignment extracted |
| `scripts/10_rebuild_placement_tree.sh` | Placement tree inferred |
| `scripts/11_test_novel_akk_placement.py` | Monophyly, exclusion, and sister relationship tested |
| `jobs/12_run_iqtree_placement.sh` | Placement re-inferred by maximum likelihood with branch support |
| `scripts/13_read_iqtree_support2.py` | Support values decoded at the three key nodes |
| `scripts/14_compare_iqtree_fasttree.py` | Topology compared between the two inference methods |
| `scripts/15_fig1_akk_clades_iqtree.py` | Tree and clade colour annotations for the placement figure |

### Step 3. Genus delimitation and genome description

Delimitation is assessed on three independent axes, since no single threshold
settles a genus boundary. Amino acid identity and percentage of conserved
proteins are computed against a representative of every neighbouring genus,
with the conserved-protein statistic recomputed properly and then tested for
stability across alternative representatives, because a single-reference value
sitting near the decision threshold cannot carry a delimitation. Nucleotide
identity is used for species-level structure within the candidate set only,
after a control establishes what the tool can and cannot resolve at this
distance. Genome size, GC content, and the features required by the MIMAG
standard are then measured for the description.

**Main result.** Amino acid identity between the candidate genus and
Akkermansia is 56.6% over 1,445 reciprocal best hits, against 48.4 to 50.7%
for the other neighbouring genera, placing the candidate closer to
Akkermansia than any other genus while remaining well outside it. The
conserved-protein statistic is 51.8% against a single reference but ranges
from 46.3 to 54.2% across six alternative representatives and declines with
genome size, so it straddles the conventional threshold and cannot carry the
delimitation on its own; the delimitation therefore rests on relative
evolutionary divergence, topology, and amino acid identity. A nucleotide
identity control is decisive about its own limits: 90 of 94 genuinely named
amphibian Akkermansia also returned no hit against the same reference set, so
the absence of hits from the candidate genomes is a tool floor and carries no
information. Within the candidate set, 105 genomes resolve into 17 clusters at
95% identity across 11,025 pairs. Median genome size is 3.04 Mb (1.42 to 4.01)
and median GC is 49.05% (43.46 to 52.72). Five of the 105 genomes meet the
MIMAG high-quality standard. The designated type genome is 3.29 Mb at 49.07%
GC, 100% complete with 0.17% contamination, carries 5S, 16S, and 23S rRNA and
21 tRNA amino acid types, and derives from a wild newt; it ranks 23rd of 105
by distance to the set median, so it is a high-quality representative rather
than the most typical genome. Reduced genomes are distributed according to
how reduction is defined: under a rule flagging any cluster below either
threshold, median GC under 45.5% or median size under 2.6 Mb, seven clusters
totalling 34 genomes qualify, but only one is reduced in both dimensions,
holding 18 genomes at 2.19 Mb and 43.7% GC. The remaining six clusters, 16
genomes between them, each fall below one threshold only. Amino acid identity
within the candidate genus, measured from every genome to the type genome,
ranges from 60.3 to 99.4% with a median of 89.9%.

| File | Purpose |
|---|---|
| `jobs/16_run_aai_pocp.sh` | Identity and POCP calculation submitted |
| `scripts/17_aai_pocp.py` | Amino acid identity and an approximate POCP computed by reciprocal search |
| `scripts/18_pocp_fixed.py` | POCP recomputed against a fixed reference |
| `scripts/19_pocp_calibrate.py` | POCP calibrated on known genus pairs |
| `scripts/20_pocp_stability.py` | POCP stability across alternative representatives |
| `scripts/21_novel_akk_ani.py` | Nucleotide identity triangle and species clustering |
| `scripts/22_novel_akk_control.py` | Screening floor measured on named reference genomes |
| `scripts/23_compute_novel_size_gc.py` | Genome size and GC computed from assemblies |
| `jobs/24_run_mimag_105.sh` | rRNA and tRNA features called for the candidate set |
| `scripts/25_verify_mimag_count2.py` | High-quality genome count verified against the MIMAG standard |
| `scripts/26_check_type_genome_median.py` | Type genome checked against the set median |
| `scripts/27_check_reduced_cluster.py` | Reduced-genome clusters characterized |
| `jobs/66_run_within_genus_aai.sh` | Within-genus identity calculation submitted |
| `scripts/65_within_genus_aai.py` | Amino acid identity of every candidate genome to the type genome |

Note: the POCP reported is the one from `scripts/18_pocp_fixed.py`. The value
produced by `scripts/17_aai_pocp.py` is an approximation, as its own comment
records, and is superseded for that statistic; its amino acid identity output
stands.

### Step 4. Gene calling, pangenome, and orthology polarization

All genomes are re-called with one gene caller under identical settings, so
that no comparison inherits differences between upstream annotation pipelines.
Protein clusters are then partitioned across the candidate and sister sets,
and orthogroups are inferred across a representative set that includes
free-living outgroups. Including the outgroups is what makes the comparison
polarizable: a family absent from Akkermansia but present in both the
candidate genus and the free-living genera indicates loss on the Akkermansia
branch, whereas absence from the whole gut clade indicates loss earlier. The
categories are then annotated functionally on an independent layer.

**Main result.** Of 2,219 orthogroups assigned a polarity, the dominant event
is not on the Akkermansia branch: 1,047 were lost at the ancestor shared by
the candidate genus and Akkermansia, against 57 lost specifically on the
Akkermansia branch, while 1,026 are retained across the gut clade, 55 are
specific to the candidate genus, and 34 are enriched in Akkermansia.
Functional annotation returned assignments for 81 orthogroups across the three
non-background categories (52 Akkermansia-branch losses, 16 candidate-specific,
13 Akkermansia-enriched). Among the annotated Akkermansia-branch losses are both
glucose-6-phosphate dehydrogenase and 6-phosphogluconate dehydrogenase,
recovering by an independent route the pathway loss established in Step 6. Pangenome family counts are not
interpreted as biological quantities, because clustering identity strongly
affects them across a two-genus span.

| File | Purpose |
|---|---|
| `jobs/28_run_prodigal_199_array.sh` | Uniform gene calls, candidate and sister set |
| `jobs/29_run_prodigal_akkfam2.sh` | Uniform gene calls, family reference set |
| `jobs/30_run_prodigal_gtdbakk.sh` | Uniform gene calls, GTDB set |
| `jobs/31_run_prodigal_outgroups.sh` | Uniform gene calls, free-living outgroups |
| `scripts/32_pangenome_prep.py` | Group assignments and pooled proteins assembled |
| `jobs/33_run_ppanggolin.sh` | Pangenome partitioning |
| `scripts/34_ppang_analyze_id40.py` | Partitions summarized at the corrected identity threshold |
| `scripts/35_build_orthofinder_set.py` | Representative proteome set assembled |
| `jobs/36_run_orthofinder.sh` | Orthogroup inference |
| `scripts/37_orthogroup_polarity.py` | Orthogroups polarized against free-living outgroups |
| `scripts/38_extract_polarity_reps.py` | Representative proteins extracted for annotation |
| `jobs/39_run_eggnog_polarity.sh` | Functional annotation of the polarized orthogroups |
| `scripts/40_summarize_polarity_eggnog.py` | Annotated categories summarized |

### Step 5. The conserved mucin-degradation repertoire

Mucin-degrading enzyme families are scored across the candidate genus and four
lineage-matched Akkermansia groups under one uniform evidence filter. The
comparison is directed at a specific published claim, that mucin degradation
was gained at the last common ancestor of Akkermansia and constitutes the
innovation founding the genus. That claim was made using free-living
Verrucomicrobiales as outgroups, all of which fall outside the 292-tip node
identified in Step 2, so the branch it describes skips the entire sister
lineage tested here.

**Main result.** The candidate genus carries the mucin-degradation toolkit in
full. Across 12 families clearing a 10% floor it carries every one at 87 to
100% prevalence, with a median of 94.8% against 97 to 99% in the Akkermansia
groups, and GH2, GH20, GH95, and GH109 at 100%. Four families reported absent
from both the Akkermansia ancestor and Akkermansia genomes are absent here
too, which confirms the two datasets are comparable. All five families
reported as gained at the last Akkermansia common ancestor are present in the
sister genus: GH20 at 100%, GH29 at 92.4%, GH33 at 97.1%, GH35 at 96.2%, and
GH95 at 100%. The parsimonious reading is that this repertoire was inherited
from the common ancestor of the gut-associated Akkermansiaceae rather than
gained on the Akkermansia branch.

| File | Purpose |
|---|---|
| `jobs/41_run_novel_genus_function.sh` | Functional panel submitted |
| `scripts/42_novel_genus_function.py` | Mucin family prevalence across groups |
| `scripts/43_build_mucin_matrix.py` | Mucin presence matrix assembled |
| `scripts/44_fig3_mucin_5group.py` | Mucin conservation heatmap |

Output: `figures/Figure_mucin_conservation.pdf`

### Step 6. Loss of the oxidative pentose phosphate pathway in Akkermansia

The three genes of the oxidative pentose phosphate pathway are scored across
the candidate genus, three separate Akkermansia collections, and the
free-living genera. Denominators are taken from database directory counts
rather than from hit parsing, since a parse wobble in the denominator would
change a prevalence even when the numerator is zero. Physical linkage is then
measured, because two genes that are adjacent and co-oriented across most
genomes carrying both are better evidence of a functional unit than presence
alone, and an assembly-quality control tests whether absence could be an
artifact of genome recovery.

**Main result.** The pathway is present throughout the candidate genus, with
glucose-6-phosphate dehydrogenase in 100 of 105 genomes, 6-phosphogluconate
dehydrogenase in 101 of 105, and the dehydrogenase accessory subunit in 97 of
105, and it is present across the free-living genera. It is absent from every
Akkermansia genome tested: 0 of 94 amphibian, 0 of 172 family reference, and 0
of 60 GTDB genomes, 0 of 326 in total. Because the pathway is present in both
the sister genus and the free-living outgroups, the distribution is read
directly as loss on the Akkermansia branch, with no model required. In the 95
candidate genomes carrying both, the dehydrogenase and its accessory subunit
lie on the same contig in 91, with a median intergenic gap of 20 bp (range 12
to 60, all under 100 bp) and co-orientation in all 91, while
6-phosphogluconate dehydrogenase is unlinked. Assembly quality does not
explain the absence: the gut genomes are 98.0% MAGs at median completeness
96.4%, while the free-living genomes that retain the pathway are 90.7% MAGs at
median completeness 93.1%, so the set lacking the pathway is the
better-assembled one.

| File | Purpose |
|---|---|
| `scripts/45_build_ppp_presence.py` | Pathway genes scored across collections |
| `scripts/46_fix_ppp_326.py` | Presence table rebuilt with corrected collection denominators |
| `scripts/47_check_akk_denominator.py` | Collection denominators verified against database directories |
| `scripts/48_measure_ppp_operon2.py` | Operon adjacency and intergenic gap measured |
| `scripts/49_check_ppp_strand.py` | Operon co-orientation checked |
| `scripts/50_fig_ppp_final.py` | Pathway prevalence heatmap |
| `scripts/51_fig_ppp_operon_final.py` | Operon diagram |
| `scripts/61_check_mag_vs_isolate.py` | Assembly-type and completeness control on absence claims |

Output: `figures/Figure_ppp_loss.pdf`, `figures/Figure_ppp_operon.pdf`

### Step 7. Divergence in the carbohydrate-active enzyme repertoire

Enzyme family prevalence is censused across all annotated Verrucomicrobiota
under one uniform evidence filter, after an earlier census was found to have
applied different filters to two halves of the same dataset. Families are then
compared across five lineage-matched groups, so that any difference is within
Akkermansiaceae rather than between genera of different habitat, and finally
polarized against the free-living genera in the same way as the orthogroups in
Step 4.

**Main result.** The candidate genus and Akkermansia differ sharply in
chitin-related families, in opposite directions for the two principal
families. GH75 is present in 80.0% of candidate genomes against 20.2% of
amphibian Akkermansia, 6.6% of wall lizard Akkermansia, 67.6% of mammal
Akkermansia, and 73.8% of GTDB Akkermansia, while GH18 runs the other way at
11.4% in the candidate genus against 84.0%, 91.2%, 59.2%, and 61.9%. GH46 is
present in 29.5% of candidate genomes and absent from every Akkermansia group.
A housekeeping control family is flat across the same groups, so the contrast
is not an assembly artifact. Against the free-living genera, GH75 is common
throughout (61.1 to 100%) and GH18 is not (0 to 40.7%), and a set of further
families including GH92, GH38, GH139, GH120, PL33, GH154, CBM91, and GH141 is
carried by 53 to 67% of candidate genomes and by the free-living genera while
sitting between 0.3% and 1.4% across Akkermansia. Gain of families in the
candidate genus is not reported here: two criteria applied to the same 294
families return different counts, and neither is adopted. The direction of the
GH75
difference is reported as a prevalence contrast and not as a loss, because a
formal reconstruction of that family did not meet its pre-registered
criterion. The five groups cover 344 of the 346 annotated Akkermansia genomes.

One species cluster within the candidate genus, 18 genomes recovered only from
newts, is reduced in both genome size and GC content relative to the rest of
the genus. Enzyme family richness is proportional to genome size: the cluster
carries a median 44 families per genome against 64 for the remaining 87, but
the two groups are indistinguishable once size is accounted for, at 20.0
against 20.4 families per Mb. Enzyme protein copy-number density is instead
markedly lower: 87 enzyme proteins against 171, which is 51% of the comparison
group against 71% of its genome size, and 39.4 against 54.5 proteins per Mb.
Reduction in this lineage has removed gene copies rather than gene families. This places three reductions
at three depths in the same clade, the third having already produced a narrower
repertoire than Akkermansia itself. Every family absent from Akkermansia but
retained across the candidate genus is also absent from this cluster, but the
nesting is largely expected rather than diagnostic: both sets occupy the same
narrow prevalence band, and a null matched on prevalence recovers a mean of
6.65 of the 9 observed overlaps. The analysis covers enzyme families only, since
the orthogroup work was conducted on cluster-level representative proteomes
rather than individual genomes.

| File | Purpose |
|---|---|
| `scripts/52_patch_census_v3.py` | Census patched to a uniform evidence filter |
| `scripts/53_gh75_verru_census_v3.py` | Family prevalence census under the uniform filter |
| `jobs/54_run_census_v3.sh` | Census submitted |
| `scripts/55_patch_join_v3.py` | Family join patched to the corrected census |
| `scripts/56_join_herptile_verru_family_v3.py` | Family assignments joined onto the census |
| `jobs/57_run_join_v3.sh` | Family join submitted |
| `scripts/58_chitin_panel_5groups.py` | Chitin panel prevalence across lineage-matched groups |
| `scripts/59_cazy_polarize_freeliving.py` | Families polarized against free-living outgroups |
| `scripts/60_cazy_by_freeliving_genus.py` | Family prevalence by free-living genus |
| `scripts/62_fig_cazy_heatmap.py` | Enzyme family heatmap |
| `scripts/63_fig_cazy_panels.py` | Enzyme family panels |
| `scripts/67_build_cazy_per_genome_matrix.py` | Enzyme family presence assembled per genome from the annotation output |
| `scripts/68_reduced_cluster_cazy_loss.py` | Enzyme families lost in the reduced species cluster |
| `scripts/69_reduced_cluster_matched_null.py` | Loss overlap tested against a prevalence-matched null |
| `scripts/72_reduced_cluster_proportionality.py` | Enzyme repertoire of the reduced cluster tested against genome size |

Output: `figures/Figure_cazy_heatmap.pdf`, `figures/Figure_cazy_polarity_panels.pdf`

## Software

IQ-TREE 2.2.2.6, FastTree, GTDB-Tk (GTDB r226), Prodigal, GUNC 1.0.6 with
proGenomes 2.1, DIAMOND 2.0.4 under GUNC, PPanGGOLiN 2.2.3,
OrthoFinder, MMseqs2, DIAMOND, skani 0.2.2, dbCAN, HMMER 3.3.2, eggNOG-mapper
2.1.9, barrnap 0.9, tRNAscan-SE 2.0.12, CheckM, ete3, Python 3.9.

## Repository layout

    jobs/      SLURM submission scripts
    scripts/   analysis code
    figures/   final manuscript figures
    tables/    supplementary tables
    config/    input lists and batch files
    docs/      extended notes
