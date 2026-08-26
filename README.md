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

Genomes with no GTDB genus assignment are pulled from the catalogue and checked against the classifier, to confirm the missing assignment means novelty
rather than a failed call; the evidence is reported in Step 3. Host animals are recounted under a corrected unit, since the catalogue
identifier is a sequencing run and not an animal, and captivity status is reconciled
across metadata sources, since a captivity confound would invalidate later comparisons.

**The set.** Within Akkermansiaceae, 107 genomes carry no genus assignment. Two, both from one tortoise, sit at
53.3% amino acid identity to the type genome (53.26 and 53.28% over 1,193 and 1,037
reciprocal best hits) against 60.3 to 99.4% for the rest. The seven-point gap is empty,
so they were excluded on evidence rather than by a decision taken in advance. That
leaves 105 MAGs: 54 from the herptile catalogue, 51 from EHI, representing 66 host
animals.

**Host range.** All 105 are amphibian, and this was tested rather than assumed. Across the
full census, 107 Akkermansiaceae genomes carry no genus assignment: 105 from amphibians and the
two tortoise genomes excluded above. *Akkermansia* itself is recovered from 350 genomes in the
same census, spanning 137 reptile, 94 amphibian, 71 mammal, five bird and one fish, so the family
is sampled across vertebrate classes while the candidate genus is not.

The Earth Hologenome Initiative collection provides an independent test. Of its 17,666
MAGs, spanning 33 phyla and comprising 8,208 from reptiles, 6,244 from mammals, 3,020 from
amphibians and 194 from birds, 368 are Verrucomicrobiota, of which 140 come from reptiles. 256
are Akkermansiaceae and 51 carry no genus assignment. All 51 are amphibian, and all 51 fall within
the candidate genus by amino acid identity to the type genome. Amphibian sampling in that
collection is narrow, however, covering only *Lissotriton helveticus* and *Calotriton asper*, so
the result bounds the absence of the genus from other host classes rather than establishing its
breadth across amphibians.

**Quality.** Median completeness 96.6% (51.2 to 100.0). These values are of mixed
provenance: CheckM v1 for the 54 herptile MAGs, upstream quality control for the 51 EHI
genomes, and the catalogue records no tool column. All 105 matched the annotation table,
none unmatched. All 105 are from wild animals.

**Chimerism.** 102 of 105 pass GUNC at the default threshold, against 86 of 94 amphibian
*Akkermansia* screened in the same run, and no genome in either group meets a strict
chimera call. All 105 fall below a reference representation score of 0.5, which reflects
having no named representative in the reference set rather than contamination. Within each
arm, genomes that fail carry fewer enzyme families than genomes that pass (candidate 52.0
against 61.0, *Akkermansia* 45.5 against 48.0), the opposite of what chimeric inflation of
the repertoire would produce.

| File | Purpose |
|---|---|
| `scripts/01_check_akk_genus.py` | Genus composition of the Akkermansiaceae reference set |
| `scripts/02_resolve_akk_genus.py` | Genus calls resolved for amphibian Akkermansiaceae |
| `scripts/03_find_51_genus.py` | Genus field checked for the EHI subset |
| `scripts/04_check_51_gtdb_string.py` | GTDB taxonomy strings checked for the EHI subset |
| `scripts/05_recount_host_animals_corrected_unit.py` | Host animals recounted under the corrected animal unit |
| `scripts/06_verify_novel_qc.py` | Completeness and contamination verified for the candidate set |
| `scripts/07_reconcile_captivity.py` | Wild and captive status reconciled across metadata sources |
| `scripts/82_check_novel_107_vs_105.py` | Candidate set counted before and after the non-amphibian exclusion |
| `scripts/90_host_range_test.py` | Host range tested against the full census and the independent EHI collection |
| `scripts/64_tortoise_exclusion_aai.py` | Non-amphibian genomes excluded by amino acid identity to the type genome |
| `jobs/70_run_gunc_199.sh` | Chimerism screening submitted for the candidate and sister genomes |
| `scripts/71_gunc_audit_199.py` | Chimerism screening results audited by arm |

Output: `tables/TableS1_genome_quality.tsv`, `tables/TableS2_chimerism.tsv`

### Step 2. Phylogenetic placement

Placement uses the concatenated bac120 marker alignment from GTDB-Tk v2.4.1 (r226), in
the form that retains reference genomes; the user-only alignment has no references and
cannot place anything. Masking is inherited from GTDB-Tk, not chosen here: the 784-sequence
placement alignment is a subset of the full 144,146-sequence alignment at the same width,
5,036 columns. Three questions are put to the tree directly rather than read off a figure:
whether the candidates form one clade, whether any sits inside *Akkermansia*, and what the
sister group is. The tree is then re-inferred by maximum likelihood with branch support,
since the first inference gives approximate support only.

**What is in the tree.** 784 tips: 105 candidates, 187 *Akkermansia*, 290 from other
Akkermansiaceae genera (*Luteolibacter* 146, SW10 33, *Haloferula* 27, *Oceaniferula* 19,
*Roseibacillus_B* 18, *Rubritalea* 16, and a tail), and 202 from two outgroup families,
Verrucomicrobiaceae and Chthoniobacteraceae.

**Main result.** The 105 candidates form one clade of exactly 105 tips with no intruders,
and none falls inside *Akkermansia*. *Akkermansia* is monophyletic here: the ancestor of all
187 *Akkermansia* tips subtends those 187 and nothing else. The parent of the candidate clade
subtends exactly 292 tips, the 105 candidates plus the 187 *Akkermansia*, with no other
lineage inside it. The sister group is therefore the whole of *Akkermansia*, not part of it.
Both trees were inferred unrooted; the sister statement rests on the root falling outside the
292-tip node, which the two outgroup families and the 290 free-living Akkermansiaceae place
there.

**Support.** Two programs on the same alignment give identical topology at all three nodes.
FastTree (LG, gamma rates) gives approximate likelihood support 1.0 throughout. IQ-TREE
2.2.2.6 (LG+G4, seed 20260723) gives SH-aLRT 100 and ultrafast bootstrap 100, 1000 replicates
each, for the candidate clade, the *Akkermansia* clade, and their parent. They share a model
and an alignment, so this checks tree search and support estimation, not independence.

**Bearing on Step 5.** The published claim addressed there used free-living Verrucomicrobiales
as outgroups. All 290 free-living Akkermansiaceae tips and all 202 outgroup-family tips fall
outside the 292-tip node, so that branch bypasses the candidate genus.

| File | Purpose |
|---|---|
| `scripts/08_check_bac120_msa.py` | Reference-containing bac120 alignment confirmed as the placement input |
| `scripts/09_build_akk_placement_tree.py` | Placement alignment extracted |
| `scripts/10_rebuild_placement_tree.sh` | Placement tree inferred |
| `scripts/11_test_novel_akk_placement.py` | Monophyly, exclusion, and sister relationship tested |
| `jobs/12_run_iqtree_placement.sh` | Placement re-inferred by maximum likelihood with branch support |
| `scripts/13_read_iqtree_support2.py` | Support values decoded at the three key nodes |
| `scripts/14_compare_iqtree_fasttree.py` | Topology compared between the two inference methods |
| `scripts/15_fig1_akk_clades_iqtree.py` | Tree collapsed and coloured for the figure |
| `scripts/73_fig1_render.py` | Placement figure rendered |

Output: `figures/Figure1_placement.pdf`

### Step 3. Genus delimitation and genome description

Delimitation rests on three independent axes, since no single threshold settles a genus
boundary. Nucleotide identity is used only for species structure within the set, after a
control establishes what the tool can resolve at this distance.

**Genus-level novelty.** For the 54 genomes classified in this study, the classifier
assigned family Akkermansiaceae and left the genus empty by relative evolutionary
divergence: RED 0.8262 to 0.8306, median 0.8289, with no closest reference genome for any
of them. All 54 carry the same classification method, taxonomic novelty determined using
RED. The empty genus is therefore a positive determination of novelty at genus level, not a
failed or truncated call. The 51 EHI genomes carry family and genus assignments in the
catalogue but no RED value in that table, so this range covers the 54 classified here.

**Amino acid identity.** 56.6% to *Akkermansia* over 1,445 reciprocal best hits, against
48.4 to 50.7% for the other neighbouring genera. This establishes the ordering rather than the
boundary: *Akkermansia* is the closest named genus by a clear margin, and every other
neighbour is further away. The boundary itself is set by relative evolutionary divergence and
topology, since published amino acid identity cut-offs for the genus rank are broad and vary
between clades.

**Conserved proteins.** 51.8% against a single reference, but 46.3 to 54.2% across six
alternatives, and the value declines with genome size. It straddles the conventional 50% threshold and cannot carry the delimitation
alone. The threshold does discriminate within this family, however: 28 accepted
Akkermansiaceae genus pairs return 20.4 to 48.6%, none of them reaching 50, so the candidate
value sits above every accepted genus pair tested here. The delimitation rests on relative
evolutionary divergence, topology, and amino acid identity, with this one metric dissenting.

**Nucleotide identity control.** No candidate genome returned a nucleotide identity hit
against the named reference set, which would ordinarily suggest extreme divergence. It does
not: 90 of 94 genuinely named amphibian *Akkermansia* returned no hit against the same
references. The absence is a limit of the screening tool at this distance and carries no
information about the candidate genomes, so nucleotide identity is used here only for species
structure within the set.

**Species structure.** 105 genomes resolve into 17 clusters at 95% identity across 5,460
distinct pairs. Within-genus amino acid identity to the type genome runs 60.3 to 99.4%, median
89.9%.

**Genome description.** Median size 3.04 Mb (1.42 to 4.01), median GC 49.05% (43.46 to
52.72). Three genomes meet the MIMAG high-quality standard. Five satisfy it on the
stated criteria, but the ribosomal RNA of two proves not to belong to the genome: both carry only
a short partial 16S fragment, 456 and 526 bp, on contigs of 9.5 kb and 1.9 kb, sharing 0.03 and
0.09 of their 12-mers with verified genus sequences. The count is limited by ribosomal RNA
recovery rather than by completeness, since the standard requires 5S, 16S and 23S to be assembled
and short-read metagenome assemblies recover ribosomal operons poorly.

**Type genome.** 3.29 Mb, 49.07% GC, 100% complete, 0.17% contamination, 5S, 16S and 23S
rRNA, 21 tRNA amino acid types, from a wild newt. It carries two complete ribosomal operons on
separate contigs. The operon on the 65 kb contig is the genome's own, sharing 0.79 of its 12-mers
with verified genus sequences; the operon on a 7.7 kb contig shares 0.17 and is a misbinned
fragment, which is excluded from the deposited assembly. It ranks 23rd of 105 by distance to the set median, so it is a
high-quality representative rather than the most typical genome. It belongs to a
15-genome species cluster, not to the reduced cluster described below.

**Cluster-to-cluster divergence.** A concurrent genome-scale revision of the family
proposes genus thresholds of about 65% amino acid identity and 50% conserved proteins, and
applies them to divide *Akkermansia* on inter-clade identities of 60.8 to 63.7%. Measured
against the same thresholds, the candidate genus sits at 56.6% to *Akkermansia*, below both,
and its conserved-protein value of 51.8% is comparable to the 50.9 to 53.3% on which part of
that division rests.

Applied within the candidate set, the threshold identifies internal structure. Amino acid
identity between representatives of the 17 species clusters spans 59.95 to 95.79% across 136
pairs. Three units fall below 65% to the remainder: one cluster of 18 genomes, a group of
three clusters totalling 7 genomes, and one divergent genome. Values between 60 and 65% form
a continuum rather than a discrete boundary, and the largest gaps in the matrix, 7.05 and 6.3
points, fall at high identity where species boundaries lie. The *Akkermansia* genomes examined
in that revision spanned 60.8 to 99.96% and were treated as one genus, although those
estimates were generated with EzAAI whereas the values here use reciprocal best hits under
DIAMOND, so the two are not strictly calibrated against each other. The clade is retained as a
single genus on the basis of topology and relative evolutionary divergence, with identity
metrics as quantitative support; whether the most divergent clusters warrant separate
recognition is left to denser sampling.

**Reduced genomes.** The count depends on the criterion. Under an either-or rule, median
GC below 45.5% or median size below 2.6 Mb, seven clusters and 34 genomes qualify. Only
one cluster is reduced in both dimensions: 18 genomes at 2.19 Mb and 43.7% GC. The other
six, 16 genomes between them, fall below one threshold only.

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
| `scripts/91_verify_hq_16s.py` | Ribosomal RNA of the high-quality genomes checked for provenance |
| `scripts/26_check_type_genome_median.py` | Type genome checked against the set median |
| `scripts/27_check_reduced_cluster.py` | Reduced-genome clusters characterized |
| `jobs/66_run_within_genus_aai.sh` | Within-genus identity calculation submitted |
| `scripts/65_within_genus_aai.py` | Amino acid identity of every candidate genome to the type genome |
| `scripts/81_type_genome_candidates.py` | Type genome candidates ranked against the description criteria |
| `scripts/87_pick_cluster_reps.py` | Species cluster representatives selected by completeness |
| `jobs/88_run_cluster_aai.sh` | Cluster-to-cluster identity searches submitted |
| `scripts/89_compute_cluster_aai.py` | Cluster-to-cluster identity matrix computed from reciprocal best hits |
| `scripts/83_step3_red_values.py` | Genus-level novelty confirmed from classifier divergence values |

Note: the POCP reported is the one from `scripts/18_pocp_fixed.py`. The value
produced by `scripts/17_aai_pocp.py` is an approximation, as its own comment
records, and is superseded for that statistic; its amino acid identity output
stands.

### Step 4. Gene calling, pangenome, and orthology polarization

All genomes are re-called with one gene caller, so no comparison inherits differences
between upstream pipelines. Orthogroups are then inferred across 60 representative proteomes including free-living
outgroups. Each proteome represents a species cluster rather than a single genome, so every
percentage below is a percentage of clusters. The outgroups make the comparison polarizable: absent
from *Akkermansia* but present in both the candidate genus and the free-living genera
means loss on the *Akkermansia* branch, while absence from the whole gut clade means
loss earlier.

**Where the losses fall.** 2,219 of the 8,812 orthogroups met one of the polarization
criteria. Each criterion is explicit: loss on the *Akkermansia* branch requires presence in at
least 70% of candidate clusters, at most 10% of *Akkermansia*, and at least three free-living
genera; loss at the shared ancestor requires at most 5% in both gut genera and presence in at
least five free-living genera. The dominant pattern is not on the *Akkermansia* branch: 1,047 orthogroups have a
distribution consistent with loss on the branch ancestral to both gut genera, against 57
consistent with loss on the *Akkermansia* branch. These are polarity assignments from
occupancy thresholds, not directly reconstructed events.
1,026 are retained across the gut clade, 55 are candidate-specific, 34 enriched in
*Akkermansia*.

**Annotation.** All five categories were annotated on the same pipeline: 996 of 1,026
retained (97.1%), 909 of 1,047 ancestral losses (86.8%), 52 of 57 branch losses (91.2%),
13 of 34 enriched, 16 of 55 candidate-specific. Categories were compared against the
retained set by Fisher's exact test with Benjamini-Hochberg correction.

**Control.** Translation and ribosomal function is depleted among the ancestral losses
(1.9% against 13.8%, ratio 0.14, q = 1e-22) and absent from the branch losses, as expected
of genes under strong constraint.

**The two events differ.** Loss on the *Akkermansia* branch is concentrated in trafficking
and secretion: 10 of the 52 annotated orthogroups, 19.2% against 3.9% of the background, ratio
4.91, q = 0.001. Motility is elevated but does not survive correction (4 orthogroups, q =
0.06). These ratios rest on small counts and should be read as such. Loss at the shared ancestor is broad: secondary metabolites
(ratio 4.75) and inorganic ion transport (ratio 1.75) are elevated, while coenzyme,
nucleotide and amino acid metabolism and replication are depleted. Carbohydrate metabolism
is the largest specific category among the ancestral losses at 7.8% but does not differ
from background (ratio 1.18, q = 0.37).

**Independent recovery.** Both glucose-6-phosphate dehydrogenase and 6-phosphogluconate
dehydrogenase appear among the annotated branch losses, recovering the Step 6 result by a
different route.

**Carbohydrate enzymes among the ancestral losses.** Twelve of the 909 annotated ancestral
losses carry a carbohydrate-active enzyme assignment, four of them involving GH18, including
one architecture combining a lytic polysaccharide monooxygenase with chitin-binding modules.
Some chitin machinery was therefore already lost at the shared ancestor, which cuts across
Step 7, where GH18 is the family enriched in *Akkermansia*.

**Two limits on reading this.** The background is the set retained across both gut genera, not the genome, so
these are differences between lost and retained genes rather than enrichment against a neutral
expectation. And the two gut genera are sisters, sharing a branch that no single free-living
genus shares, so some shared absence is expected under any history: the shared ancestor is a
position on the tree, not a demonstration that gut life caused the loss. And unassigned function is itself elevated among
the ancestral losses (27.9% against 13.9%), which partly reflects that genes present only
in free-living Verrucomicrobiota are less well characterized.

**Not interpreted.** The pangenome partitioning is kept in the record for provenance, but
its family counts are not read as biological quantities: clustering identity strongly affects
them across a two-genus span. The orthogroup analysis, not the pangenome, carries the result.

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
| `scripts/74_extract_polarity_reps_B.py` | Representative proteins for the ancestral losses |
| `jobs/75_run_eggnog_polarity_B.sh` | Ancestral losses annotated functionally |
| `scripts/76_summarize_polarity_B.py` | Functional categories of the ancestral losses summarized |
| `scripts/77_extract_polarity_reps_any.py` | Representative proteins for any polarity category |
| `jobs/78_run_eggnog_polarity_shared.sh` | Retained orthogroups annotated as the background |
| `scripts/79_polarity_cog_enrichment.py` | Loss categories tested against the retained background |
| `scripts/80_fig_orthogroup_polarity.py` | Orthogroup occupancy figure |

Output: `figures/Figure_orthogroup_polarity.pdf`

### Step 5. The conserved mucin-degradation repertoire

Mucin-degrading enzyme families are scored across the candidate genus and four
lineage-matched Akkermansia groups under one uniform evidence filter. The
comparison is directed at a specific published claim, that mucin degradation
was gained at the last common ancestor of Akkermansia and constitutes the
innovation founding the genus. That claim was made using free-living Verrucomicrobiales as outgroups, all of
which fall outside the 292-tip node identified in Step 2, so the branch it describes skips the
entire sister lineage tested here. The annotation method matches: the same enzyme family
models at the same evidence threshold, with an additional coverage filter applied here, so
this comparison is at least as conservative as the one it addresses.

**Main result.** The mucin-associated repertoire characteristic of Akkermansia is largely
conserved in the candidate genus. Twelve families were retained for the panel by a 10%
prevalence floor, which is a selection criterion rather than a result. Within that panel the candidate genus carries every family, at 86.7 to
100% prevalence, median 94.8%, against group medians of 97.2 to 98.5% in the four Akkermansia
groups, with GH2, GH20, GH95 and GH109 at 100%. Four families reported absent from both the Akkermansia ancestor and
Akkermansia genomes are at or below 1.9% in both gut lineages here (GH85, GH98 and GH129 at
0%, GH101 at 1.9% in the candidate genus and 0% in Akkermansia), so the two datasets agree
where they should. Three of the four are carried by free-living relatives (GH85 9.6%, GH129
4.6%, GH98 3.9%), so their absence is specific to the gut lineages rather than general to the
family. All five families
reported as gained at the last Akkermansia common ancestor are present in the
sister genus: GH20 at 100%, GH29 at 92.4%, GH33 at 97.1%, GH35 at 96.2%, and
GH95 at 100%. The panel is not uniform inside Akkermansia either: GH27 falls to 9.9% in
mammal-associated genomes and 26.2% in the GTDB set, against 86.7% here, and GH123 falls to
56.3% in mammal-associated genomes. The parsimonious reading of the shared families is that
the repertoire was present in the common ancestor of the gut-associated Akkermansiaceae and
inherited. Presence in the sister genus makes inheritance the simpler explanation; it does not
exclude independent acquisition, and no ancestral-state reconstruction is claimed here.

| File | Purpose |
|---|---|
| `jobs/41_run_novel_genus_function.sh` | Functional panel submitted |
| `scripts/42_novel_genus_function.py` | Mucin family prevalence across groups |
| `scripts/43_build_mucin_matrix.py` | Mucin presence matrix assembled |
| `scripts/44_fig3_mucin_5group.py` | Mucin conservation heatmap |

Output: `figures/Figure_mucin_conservation.pdf`

### Step 6. Loss of the oxidative pentose phosphate pathway in Akkermansia

The three pathway genes are scored across the candidate genus, three separate
*Akkermansia* collections, and the free-living genera. Denominators come from database
directory counts rather than hit parsing, since a parse wobble would change a prevalence
even when the numerator is zero.

**Present in the candidate genus.** Glucose-6-phosphate dehydrogenase in 100 of 105,
6-phosphogluconate dehydrogenase in 101 of 105, the accessory subunit in 97 of 105.

**Present in the free-living genera.** Across 146 free-living Akkermansiaceae genomes,
glucose-6-phosphate dehydrogenase is in 142, 6-phosphogluconate dehydrogenase in 136, and the
accessory subunit in 122. Every genus carries all three: Haloferula 15 of 15, Rubritalea 11 of
11, Luteolibacter 93 to 111 of 113 depending on the gene, and both Roseibacillus lineages
positive.

**Absent from *Akkermansia*.** 0 of 94 amphibian, 0 of 172 family reference, 0 of 60
GTDB. 0 of 326 in total. Because the pathway is present in both the sister genus and the outgroups, this
reads directly as loss on the *Akkermansia* branch, with no model required.

**Detection is not limiting.** A positive control was run through the same search: a GH20
query from the candidate genus, using a gene family known to be present in *Akkermansia*. It
returns hits at 59.5% amino acid identity in amphibian *Akkermansia* (524 aligned residues,
E = 1e-189) and 62.8% in the GTDB set (529 residues, E = 1e-205), and recovers homologues in
the free-living genera down to 36.0% identity. The search therefore detects diverged
homologues of a gene that is present, so a zero for the pathway genes reflects absence rather
than a failed search.

**Physical linkage.** In the 95 candidate genomes carrying both, the dehydrogenase and
its accessory subunit share a contig in 91, with a median intergenic gap of 20 bp (12 to
60, all under 100) and co-orientation in all 91. 6-phosphogluconate dehydrogenase is unlinked. The four pairs on different
contigs are not evidence against linkage: contig breaks in metagenome-assembled genomes are
assembly artefacts, so the same-contig count is a floor.

**Not an assembly artifact.** The gut genomes are 98.0% MAGs at median completeness
96.4%; the free-living genomes that retain the pathway are 90.7% MAGs at 93.1%. The set
lacking the pathway is the better-assembled one.

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
| `scripts/84_id_ppp_trio.py` | Pathway query proteins identified by sequence motif and length |
| `scripts/85_ppp_freeliving_per_gene.py` | Pathway prevalence per gene across the free-living genera |
| `scripts/86_ppp_operon_tee.py` | Operon adjacency, intergenic gap and co-orientation written to file |

Output: `figures/Figure_ppp_loss.pdf`, `figures/Figure_ppp_operon.pdf`

### Step 7. Divergence in the carbohydrate-active enzyme repertoire

Family prevalence is censused across all annotated Verrucomicrobiota under one evidence
filter, after an earlier census was found to have applied different filters to two halves
of the same dataset. Families are compared across five lineage-matched groups, so any
difference is within Akkermansiaceae rather than between genera of different habitat, then
polarized against the free-living genera as in Step 4.

**Chitin-related families run in opposite directions.** The four comparison groups are
amphibian *Akkermansia* (n = 94), *Podarcis* wall lizard (137), mammal (71) and GTDB reference
(42). GH75 is in 80.0% of candidate genomes against 20.2%, 6.6%, 67.6% and 73.8% in those four.
GH18 runs the other way: 11.4% against 84.0%, 91.2%, 59.2% and 61.9%. GH46 is in 31 of 105
candidate genomes (29.5%) and in none of the 344 *Akkermansia* genomes across all four groups.
The control family CE11 is flat over the same groups at 88.6, 89.4, 89.1, 98.6 and 100%, so the
contrast is not an artefact of assembly or annotation depth.

**Repertoire size.** The candidate genus carries a median of 61.0 enzyme families per genome
(interquartile range 46.0 to 66.0, n = 105) against 48.0 in *Akkermansia* (42.0 to 53.0, n =
346). Removing the three candidate genomes that did not pass the chimerism screen leaves the
median and the interquartile range unchanged.

**Against the free-living genera.** GH75 is common throughout (61.1 to 100%) and GH18 is not
(0 to 40.7%). A further set including GH92, GH38, GH139, GH120, PL33, GH154, CBM91 and
GH141 is carried by 53 to 67% of candidate genomes and by the free-living genera, while
sitting between 0.3% and 1.4% across *Akkermansia*.

**Not reported.** Gain of families: two criteria on the same 294 families return different
counts and neither is adopted. The direction of the GH75 difference: reported as a
prevalence contrast, not a loss, because the formal reconstruction did not meet its
pre-registered criterion.

The five groups cover 344 of the 346 annotated *Akkermansia* genomes.

### A recent reduction inside the candidate genus

One species cluster, 18 genomes from newts only, is reduced in both size and GC.

**Family richness scales with genome size.** At 20.0 families per Mb against 20.4 in the
other 87 genomes, family density is indistinguishable (Mann-Whitney p = 0.84). The raw counts,
44 families against 64, follow from the smaller genome (2.19 Mb against 3.09) and should not be
read as selective family loss.

**Copy-number density does not.** 87 enzyme proteins against 171, which is 51% of the
comparison group against 71% of its genome size, and 39.4 against 54.5 proteins per Mb.
Reduction here has removed gene copies rather than gene families.

**Three depths.** That places three reductions in the same clade: 1,047 orthogroups lost at
the shared gut ancestor, 57 on the *Akkermansia* branch, and this cluster, whose median of 44
families per genome is already narrower than the 48 of *Akkermansia* itself.

**The nesting is expected, not diagnostic.** Every family absent from *Akkermansia* but
retained across the candidate genus is also absent from this cluster. But both sets sit in
the same narrow prevalence band, and a prevalence-matched null recovers a mean of 6.65 of
the 9 observed overlaps.

**Scope.** Enzyme families only. The orthogroup work used cluster-level representative
proteomes rather than individual genomes, so the same question cannot be asked there.

**Not chimerism.** A wider repertoire is what a chimeric bin would produce, so the chimerism
screen in Step 1 is the control for this section. Genomes that fail it carry fewer enzyme
families than genomes that pass, within each arm, and the candidate arm failed less often than
the comparison arm.

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
