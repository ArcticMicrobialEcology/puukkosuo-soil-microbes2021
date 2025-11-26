#!/bin/bash
#SBATCH --job-name=setup_metadata_databases
#SBATCH --error=1_setup_metadata_databases_%A_%a_err.txt
#SBATCH --output=1_setup_metadata_databases_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time 71:00:00
#SBATCH --account project_2005827

# set up  sample information and compile some installed databases needed by the workflow

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the metagenomics and metatranscriptomics raw data directories
mg_raw_data_dir="/scratch/project_2009164/2_OULANKA/1_RAW/metagenomics"
mt_raw_data_dir="/scratch/project_2009164/2_OULANKA/1_RAW/metatranscriptomics"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts"

# the directory for the external installed R packages, such as taxonomizr (required)
r_package_dir="$software_base_dir/project_rpackages_430"

# the location of the functional 51 metabolic marker gene database
met_mark_annot_dir="$database_base_dir/metmarkdb" 
metmark_database_file="$database_base_dir/metmarkdb/Funcgenes_51_Dec2021.fasta"

# the file for some manually observed additional eukaryotic proteins in the functional 51 metabolic marker gene database
met_mark_euk_prot_man_file="$database_base_dir/metmarkdb/Some_manually_observed_eukaryotic_proteins.txt"

# the downloaded ncbi protein to tax id sql database file
ncbi_database_file="$database_base_dir/ncbi_protein/accessionTaxa.sql"

# load the R environment
module load r-env/432 # version

# run a small R script to set up sample names, metadata, etc.
# first input parameter needs to be the metadata directory for the study. The metadata directory
# needs to contain the essential information files, such as the Oulanka_ACAP_study_site.xlsx, Puukkosuo_vegetatation_clusters.xlsx

echo "[$(date -Is)] Setting up sample names and such, modify the R_scripts/1_SETUP_sample_names_etc.R as necessary to change the initialization..."
Rscript "$r_script_dir/1a_SETUP_sample_names_etc.R" "$project_root/metadata" "$mg_raw_data_dir" "$mt_raw_data_dir"
echo "*******************************************************"
echo "*******************************************************"

# run a small R script to parse the marker gene database metadata.
# requires as input parameter the directory containing the database and the required information.
# e.g. all_lengths_51_dbs_subgroup_Jan2022.txt, Funcgenes_51_Dec2021.fasta, Some_manually_observed_eukaryotic_proteins.txt

echo "[$(date -Is)] Setting up the metabolic marker database information..."
Rscript "$r_script_dir/1b_SETUP_metabolic_marker_database.R" "$met_mark_annot_dir"
echo "*******************************************************"
echo "*******************************************************"

# run a small R script to remove the eukaryotic proteins from the metabolic marker gene database
echo "[$(date -Is)] Filtering the metabolic marker database for eukaryotic proteins and leaving only prokaryotic proteins..."
Rscript "$r_script_dir/1c_SETUP_remove_euk_metmarkdb.R" "$r_package_dir" "$metmark_database_file" "$ncbi_database_file" "$met_mark_euk_prot_man_file" "$met_mark_annot_dir"
echo "*******************************************************"
echo "*******************************************************"

# make a DIAMOND database for the marker genes
echo "[$(date -Is)] Making a DIAMOND database for the functional metabolic marker gene database..."
module purge
module load diamond/2.1.6 # check version if necessary

diamond makedb --in "$met_mark_annot_dir/MetMarkDB_Euk_Filtered.fasta" \
-d "$met_mark_annot_dir/MetMarkDB_prok" \
-p "$SLURM_CPUS_PER_TASK"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All set up and ready for analysis. Bye."


