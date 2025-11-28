# Puukkosuo soil microbes 2025

**Note:** This repository contains analysis codes and scripts for the manuscript:

*"Soil microbiome structure and function reflect environmental gradients and reindeer exclusion in a northern peatland (2025)"*.

## Overview
This repository contains the analysis scripts for the workflow used to process,
analyze, and visualize the data in this study. 

**Status:** Repository under active development.  
Analysis scripts are being updated and cleaned for publication. The repository will be 
expanded with more scripts and documentation as this work progresses.

## Required software
Analyses for this workflow were performed on the CSC - IT Centre for Science Puhti supercomputer with:

- Operating system: Rocky Linux 8.10 (Green Obsidian)
- OS family: RHEL-compatible (rocky/rhel/centos/fedora)
- Job scheduler: Slurm version 24.05.8 on Puhti

The following software have been utilized in the workflow:
- FastQC v0.11.9
- MultiQC v1.19
- CutAdapt v4.6
- phyloFlash v3.4.2 (preformatted SILVA 138.1 database)
- SortMeRNA v4.3.6 (sensitive reference database v4.3)
- DIAMOND blastx v2.1.6
- anvi’o v8-dev (KEGG snapshot v2025-02-04)
- MEGAHIT v1.2.9
- Prodigal v2.6.3
- Bowtie2 v2.4.4
- SAMtools v1.18
- MetaBAT2 v2.17
- GTDB-Tk v2.4.0 (database release 220)
- FastANI v1.34
- CoverM v0.7.0 
- Minimap2 v2.28
- featureCounts (Subread v2.0.6)

R environments v4.3.2 and v4.4.1 (for MOFA on a personal computer) were used
with the main following packages (additional packages could have been used for visualization):
- KEGGREST v1.4.2
- taxonomizr v0.10.6
- phyloseq v1.46.0
- microbiome v1.24.0
- MOFA2 v1.14.0
- randomForest v4.7-1.1
- caret v6.0-94
- vegan v2.6-4
- lme4 v1.1-35.1
- lmerTest v3.1-3
- fgsea v1.28.0 

## Contact information
**Tommi Välikangas**
[tommi.valikangas@luke.fi](mailto:tommi.valikangas@luke.fi)

**Jenni Hultman**
[jenni.hultman@luke.fi](mailto:jenni.hultman@luke.fi)

