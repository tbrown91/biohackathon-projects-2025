---
title: 'BioHackEU25 report: Automatic workflow for benchmarking BUSCO genes for phylogenomics'
title_short: 'BioHackEU25 #03: BUSCO for phylogenomics'
tags:
  - Phylogenomics 
  - BUSCO
  - Snakemake
  - Galaxy
authors:
  - name: Julian Chiu
    affiliation: 1
    orcid: 0009-0008-2264-0056
  - name: Solenne Correard
    orcid: 0000-0002-0554-5443
    affiliation: 3
  - name: Thomas Crow
    orcid: 0009-0004-1195-0287
    affiliation: 4
  - name: Tom Brown
    orcid: 0000-0001-8293-4816
    affiliation: 2
  - name: Anestis Gkanogiannis
    orcid: 0000-0002-6441-0688
    affiliation: 5
  - name: George He
    orcid: 0009-0002-7417-8998
    affiliation: 5
  - name: Iker Irisarri
    orcid: 0000-0002-3628-1137
    affiliation: 5
  - name: Nikita Kulikov
    orcid: 0000-0002-2566-6197
    affiliation: 5
  - name: Nikita Kulikov
    orcid: 0000-0002-2566-6197
    affiliation: 5
  - name: Tereza Manousaki
    orcid: 0000-0001-7518-0542
    affiliation: 5
  - name: Sebastian Martin
    orcid: 0000-0003-3171-4420
    affiliation: 5
  - name: Rafał Miłodrowski
    orcid: 0009-0002-4040-1231
    affiliation: 5
  - name: Cleopatra Petrohilos
    orcid: 0000-0002-8675-2726
    affiliation: 5
  - name: Gareth Price
    orcid: 0000-0003-2439-8650
    affiliation: 5
  - name: Joy (Uditi) Shah
    affiliation: 5
  - name: Nefeli Kleopatra Venetsianou
    orcid: 0009-0004-2472-2056
    affiliation: 5
  - name: Robert M. Waterhouse
    orcid: 0000-0003-4199-9052
    affiliation: 5
  - name: Jia Zhang
    orcid: 0000-0003-1161-4836
    affiliation: 5


affiliations:
  - name: First Affiliation
    index: 1
  - name: ELIXIR Europe
    ror: 044rwnt51
    index: 2
date: 13 November 2025
cito-bibliography: paper.bib
event: BH25EU
biohackathon_name: "BioHackathon Europe 2025"
biohackathon_url:   "https://biohackathon-europe.org/"
biohackathon_location: "Bad Saarow, Germany, 2025"
group: Project 03
# URL to project git repo --- should contain the actual paper.md:
git_url: https://github.com/elixir-europe/biohackathon-projects-2025/tree/main/03-automatic-workflow-for-benchmarking/buscophy/paper/
# This is the short authors description that is used at the
# bottom of the generated paper (typically the first two authors):
authors_short: First Author \emph{et al.}
---


# Introduction

As part of the BioHackathon Europe 2025, we report on our project "Automatic workflow for benchmarking BUSCO genes for phylogenomics", where we aimed to design an end-to-end pipeline which identified gene sequences within given genomes, mapped orthologous genes to each other across species, performed a filtering of paralogs and then created a species tree from these filtered alignments.

Phylogenomics is a central aspect of biodiversity genomics, as it reveals the relationships among organisms and key evolutionary processes such as introgression and gene flow. Genome-scale datasets are increasingly a reality in phylogenomics due to the availability of genomes for an ever-growing number of species. BUSCO datasets (universal single-copy orthologs [@Tegenfeldt2024] ) have become standard in assessing genome assembly completeness and are fully integrated into the pipelines of large genome consortia such as ERGA. Due to their low-copy nature, BUSCO genes are also increasingly used in phylogenomics, from genome skimming data to high-quality chromosome-scale genomes. Yet, their phylogenetic performance has not been thoroughly explored. Preliminary analyses show that BUSCO genes can recover robust phylogenetic relationships, but their single-copy nature is challenged: most BUSCO genes display varying levels of paralogy when using biodiverse species sets, and failure to account for this can negatively affect phylogenetic reconstruction.

This  project aimed to build an automatic phylogenomics pipeline using the output of the BUSCO software. Contrary to existing pipelines, we aimed to explicitly resolve paralogy events, thereby resulting in larger and more informative datasets. This pipeline will be used to benchmark the phylogenetic performance of the newly defined BUSCO lineage datasets, identifying not only the prevalence and evolutionary depth of the various paralogs but also resolving them for improved phylogenetic utility of BUSCO genes. This project aimed to result in a fully-fledged FAIR-compliant phylogenomics pipeline based on BUSCO and an assessment of the phylogenetic performance of new BUSCO gene sets (version odb12).

# Methodology

During the BioHackathon, we aimed to adapt an existing Snakemake pipeline called [buscophy](https://gitlab.leibniz-lib.de/smartin/buscophy.git), which generated phylogenetic trees of species based on identified single-copy BUSCO genes. Our main aims were:

* Ensure the pipeline is fully FAIR and deployable on all computational systems

* Add an additional paralog filtering step, to allow the inclusion of multi-copy BUSCO genes in the phylogentic reconstructions

To address the first step, we modified existing conda yamls and containers to ensure that consistent versions were used, independent of which software management system was used, and that tools were readily available in Galaxy. The latter point required building wrappers for BUSCO version 6, pal2nal and amas within the Galaxy Australia platform.

The addition of the paralog filtering step required the construction of individual gene trees, instead of immediately concatentating alignments into a supermatrix, and then processing the alignments and gene trees via a separate rule. The filtered alignments would then be used to construct a supermatrix and the ultimate species tree (Fig. 1).

![Proposed pipeline structure. Highlighted in the box are the rules that diverged from the buscophy pipeline and needed to be input.](./workflow.png){ width=200px }

## Results

Table: Pull Requests detailing tools requiring Galaxy wrappers

| Tool | URL |
| -------- | -------- |
| pal2nal | https://github.com/galaxyproject/tools-iuc/pull/7444 |
| amas | https://github.com/galaxyproject/tools-iuc/pull/7443 |
| BUSCO | https://github.com/galaxyproject/tools-iuc/pull/7139 |

![Mammalian tree](./mammal_tree.png){ width=200px }

# Discussion

## Acknowledgements

## References
