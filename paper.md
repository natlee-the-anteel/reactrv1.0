---
title: 'reactr: a modular Snakemake pipeline for comparative gene-family identification and phylogenetic analysis'
tags:
  - Python
  - Snakemake
  - bioinformatics
  - comparative genomics
  - phylogenetics
  - gene family analysis
  - plant genomics
authors:
  - name: Nathan Lee
    affiliation: 1
    corresponding: true
affiliations:
  - name: Independent Researcher, United States
    index: 1
date: 27 July 2026
bibliography: paper.bib
---

# Summary

`reactr` is a Snakemake-based [@Koster:2012; @Molder:2021] comparative
genomics pipeline
that automates a common but tedious research task: given a gene or gene
family of interest in one species (a "base" genome), find and
characterize its homologs in a second, typically less-studied species (a
"target" genome), and place them in a phylogenetic context. Starting from
a set of query protein sequences, `reactr` retrieves and prepares genome
assemblies, performs homology search using both BLAST [@Camacho:2009] and
profile hidden Markov models against Pfam [@Eddy:2011; @Mistry:2021],
extracts and aligns coding sequences with MUSCLE [@Edgar:2004], discovers
shared motifs with MEME [@Bailey:2015], infers maximum-likelihood
phylogenies with IQ-TREE [@Minh:2020], and optionally continues on to
synteny analysis with MCScanX [@Wang:2012], subcellular localization
prediction, CRISPR guide design, and PCR primer design. The pipeline is
organized as two Snakemake workflows -- `LoadDatasets.smk` for genome and
annotation acquisition, and `MainPipeline.smk`, modularized into thirteen
ordered rule files (`00_init.smk` through `12_crispr.smk`) covering each
analysis stage -- with all software dependencies pinned in a single conda
environment specification to support reproducible re-execution.

# Statement of need

Identifying and characterizing a gene family in a newly sequenced or
under-studied genome is a routine but labor-intensive step in plant and
comparative genomics: it typically involves manually chaining together a
BLAST search, a domain scan, sequence extraction, alignment, and tree
building, each run by hand with ad hoc scripts gluing the steps together.
This is exactly the workflow underlying, for example, recent
genome-wide gene-family characterization studies in non-model plant
species such as the identification of the DELLA gene family in blueberry
[@Zhou:2024]. `reactr` targets researchers -- including students and
independent researchers without access to an established bioinformatics
core -- who need to repeat this style of analysis across different
species pairs and gene families without re-writing the glue code each
time. By expressing the full workflow as a declarative Snakemake DAG,
`reactr` makes each analysis re-runnable, resumable, and independent of
the specific query or target species, so that a new characterization
task is a matter of editing a configuration file rather than writing new
scripts.

# State of the field

Several existing tools address parts of this problem. `OrthoFinder`
[@Emms:2019] provides robust, whole-proteome orthology inference across
many species and is the appropriate choice when the goal is
genome-wide orthogroup assignment. `MCScanX` [@Wang:2012] is the
standard tool for synteny and collinearity detection and is used by
`reactr` as a component, not replaced by it. Manual BLAST- and
HMMER-driven identification, as used directly in individual
gene-family studies such as @Zhou:2024, is flexible but not reusable
across projects without re-implementing the same glue logic each time.
`reactr` occupies a different niche from these: it is not a
whole-genome orthology inference tool like `OrthoFinder`, nor a
single-purpose synteny tool like `MCScanX`, but a re-runnable pipeline
for the specific, common task of taking a small bait gene set,
identifying its homologs in one target genome, and carrying that result
through alignment, motif discovery, phylogenetics, and several optional
downstream analyses (synteny, localization prediction, CRISPR guide and
primer design) in a single reproducible run. The choice to build a new,
modular Snakemake pipeline rather than extend an existing tool was
motivated by the need to chain these specific, heterogeneous
downstream steps together behind one configuration file and one conda
environment, rather than to improve on any single step in isolation.

# Software design

`reactr` splits data acquisition from analysis into two independently
runnable workflows. `LoadDatasets.smk` handles genome and annotation
retrieval and caching, itself decomposed into nine sequentially numbered
rule files (`00` through `08`), so that the same downloaded assembly can
be reused across multiple analysis runs without re-downloading.
`MainPipeline.smk` is decomposed into a further thirteen sequentially
numbered rule files, each responsible for one stage of the
analysis (initialization, query filtering, HMMER search, alignment and
phylogenetics, domain analysis, proteome preparation, annotation,
mapping with GMAP [@Wu:2005], primer design, conservation/synteny analysis,
visualization, and CRISPR design); later stages are individually
skippable or replaceable, and several -- such as subcellular
localization prediction via DeepLoc2 [@Odum:2024], which requires a
separately licensed academic install -- are gated behind an explicit, documented
opt-in rather than bundled into the default run. This design was chosen
over a single monolithic script so that (a) partial pipeline runs are
possible and cheap to resume via Snakemake's dependency tracking, (b)
individual stages can be swapped (for example, substituting a different
phylogenetic inference method) without touching unrelated code, and (c)
optional, license-restricted, or slow stages do not block the core
identification-and-phylogeny workflow that most users need by default.
All tool versions are pinned in a single `environment.yaml`, and the two
dependencies not available via conda (FlashFry [@McKenna:2018], MCScanX)
are fetched
from their upstream GitHub releases by the setup script, so that the
full dependency graph is reproducible from a clean environment.

# Research impact statement

`reactr` has been used by the developer to independently reproduce the
gene-family identification and phylogenetic placement results of three
previously published, unrelated plant gene-family studies: the DELLA
gene family in blueberry (*Vaccinium darrowii*) [@Zhou:2024], the
phytochrome (PHY) gene family in rice, and the HOOKLESS gene family in
tomato. In each case, `reactr` was configured with the same
base-species bait sequences used in the original publication and run
against the same target genome assembly, and the resulting
maximum-likelihood phylogeny was compared against a Newick
reconstruction of the originally published tree using Robinson-Foulds
distance [@RobinsonFoulds:1981]. For the PHY and HOOKLESS comparisons,
one and three false or truncated isoforms respectively (plus one extra
paralog in the HOOKLESS case) were manually removed from `reactr`'s
raw output prior to comparison; this reflects the pipeline's default
sensitivity-favoring behavior during homology search, which does not
automatically collapse isoforms, rather than an error in the underlying
homology calls, and is consistent with the behavior of other
established gene-family tools such as OrthoFinder [@Emms:2019]. After
this cleaning step, all three comparisons returned a Robinson-Foulds
distance of 0 after harmonizing tip-label naming conventions between
the two trees, indicating exact topological agreement with the
independently published results. Full validation methodology and the
per-species comparison tables are included in the repository as
`REACTR_Validation_Summary_Table.pdf`. This validation was additionally
re-run end-to-end on a clean GitHub Actions Ubuntu runner, independent
of the developer's local machine, and returned identical results,
indicating that the pipeline's output does not depend on
machine-specific configuration.

# AI usage disclosure

Generative AI (Claude, Anthropic) was used substantially during the
development of `reactr` and during the drafting of this paper. Its use
included: debugging Snakemake rule failures and dependency errors during
pipeline development; assisting with the design and implementation of
individual rule files; and drafting and revising the text of this paper
based on the author's pipeline, configuration files, and validation
results. All code was run and its output inspected by the author; the
validation results reported in the Research Impact Statement above
(the DELLA, PHY, and HOOKLESS tree comparisons and their Robinson-Foulds
distances) were independently verified by the author using the T-REX
web server, not merely generated and reported without inspection. No
part of the underlying biological or phylogenetic results reported here
was fabricated or asserted by an AI tool without independent
verification against the pipeline's actual output.

# Acknowledgements

The author thanks Zhou et al. and the authors of the PHY and HOOKLESS
studies referenced above for their original published gene-family
characterizations, against which this pipeline's output was validated.

# References
