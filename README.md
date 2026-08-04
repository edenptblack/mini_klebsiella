# Mini-analysis of prophage-encoded defence/anti-defence systems and AMR genes in *Klebsiella pneumoniae*

Small-scale pilot investigating the contents of prophage cargo across 10 complete *K. pneumoniae* genomes from Nanchang, China.

## Description

A pipeline running from NCBI accessions to analysis & interpretation:

* Genome selection
* Sequence typing
* Prophage detection
* Prophage annotation
* Defence/anti-defence system detection
* Defence system/prophage alignment
* Interpretation

Candidate genomes were drawn from complete RefSeq assemblies originating in Nanchang, China, filtered on assembly quality, then typed with mlst; ten were selected to span the available lineages. DefenseFinder was run on whole-genome proteins rather than extracted prophage proteins to investigate the distribution of systems between prophages and native bacterial DNA. Location was determined afterwards by intersecting each defence gene's coordinates with geNomad's prophage intervals.

## Findings

geNomad identified 60 prophages across the ten genomes, ranging from two to nine per genome. Pharokka assigned PHROG categories to 2,970 prophage coding sequences, of which 43.9% fell into unknown function, exceeding any structural category.  DefenseFinder identified 305 defence and anti-defence systems, of which 75 (24.6%) fell within a predicted prophage, though the proportion varied from none in the ST15 isolate to over half in one ST11 genome. Three AMR genes were detected, all on a single prophage in one genome: sul1, qacEΔ1 and aadA2. The first two constitute the 3' conserved segment characteristic of clinical class 1 integrons, so these represent one integron structure rather than three independent acquisitions; whether the prophage genuinely carries an integron or the predicted boundary has run into an adjacent element is unclear. No prophage-borne virulence factors were detected.

## Stack & methods

Python - pandas - bash - conda - NCBI Datasets - mlst - geNomad - Pharokka - DefenseFinder - AntiDefenseFinder - R - tidyverse - ggplot2 - pheatmap - RColorBrewer

## Reproducability

Environments are recorded as .yml files. Each primary tool (geNomad, Pharokka, DefenseFinder) requires its own separate environment as a precaution against dependency conflicts. Set "channel_priority: strict" with conda-forge above bioconda prior to installing.

Each primary tool requires a database; the commands for installing these databases are not recorded in the repository (may be added later), but can be found in the documentation for each tool. 

Execution order: candidate_selection.Rmd -> scripts (01-02) -> 10 candidates selected manually based on MLST scores and sample independence, written into "accession_selected.txt" -> scripts (03-06) -> final_analysis.Rmd

## Limitations

Using only 10 genomes, this analysis is purely descriptive with no statistical element and makes no broader inferences. The genomes were sourced from the same city (Nanchang, China) as a case study; while efforts were made to take isolates from different sources, the ST11 genomes cannot be assumed to be independent. Prophage boundaries are geNomad predictions not verified by att site identification. 