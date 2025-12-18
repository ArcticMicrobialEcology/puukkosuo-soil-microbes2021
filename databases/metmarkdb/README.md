# Greening Lab metabolic marker gene database

This directory contains the metabolic marker gene database used by the workflow.

## Contents

- `Funcgenes_51_Dec2021.fasta` – FASTA file of the curated concatenated metabolic marker gene sequences
- `all_lengths_51_dbs_subgroup_Jan2022.txt` – Annotation and metadata for the metabolic marker gene sequences
- `Some_manually_observed_eukaryotic_proteins.txt` – A list of some manually observed proteins in the database (information used when removing euakryotic protein from the database)

## Database construction
This database contains the Greening Lab metabolic marker gene database:

Leung PM, Greening C. Compiled Greening lab metabolic marker gene databases. 2021, DOI: 10.26180/14431208.v1

In this study, a preformatted concatenated version of the database was used, shared to the Jenni Hultman by Dr Pok Man Leung, Department of Microbiology, Biomedicine Discovery Institute, Monash University, Clayton, Victoria 3800, Australia.

## Usage
In this study, the taxonomizr R package was used to download the annotated National Center for Biotechnology Information (NCBI) protein database and remove all annotated eukaryotic proteins from the concatenated marker gene database. The NCBI taxonomy database was further supplemented with additional manually identified eukaryotic proteins (provided here).

The database is used in the following analysis steps:
- DIAMOND blastx searches
- Functional annotation

- `scripts/1_setup_sample_names_databases.sh`
- `scripts/read_based/8_mg_read_diamond_metmarkdb_combine_pairs_summarize.sh`
- `scripts/read_based/9_mt_read_diamond_metmarkdb_combine_pairs_summarize.sh`
- `scripts/contig_based/5_mg_contig_anvio_workflow.sh`
