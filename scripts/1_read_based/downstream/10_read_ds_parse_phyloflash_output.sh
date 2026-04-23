#!/bin/bash
#SBATCH --job-name=R_parse_phyloflash
#SBATCH --error=10_read_ds_parse_phyloflash_%A_err.txt
#SBATCH --output=10_read_ds_parse_phyloflash_%A_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --account=your_project_account

# a small script to combine and process the phyloflash runs for the different samples into a neat matrix and objects ready for downstream processing

# define the needed directories
# this is the project root directory under which all the subdirectories, results etc. are created.
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/read_based/downstream"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# the phyloflash result directory for metagenomics
mg_phylo_dir="$mg_main_dir/phyloflash"

# the phyloflash result directory for metatranscriptomics
mt_phylo_dir="$mt_main_dir/phyloflash"

# the metadata file with sample design
metadata_file="$project_root/metadata/Study_Metadata.RData"

# the directory containing the silva taxonomy rank mapping files
silva_taxonomy_dir="$database_base_dir/phyloflash_databases/138.1_silva_taxonomy_files"

# load the R environment
module load r-env/432 # version

# compose otu and taxonomy tables from phyloflash output - process both metagenomics and metatransciptomics
echo "[$(date -Is)] Compose OTU and taxonomy tables from phyloFlash output..."
Rscript "$r_script_dir/4_READ_parse_phyloflash_output.R" "$mg_phylo_dir/csv_files" "$mt_phylo_dir/csv_files" "$metadata_file" "$silva_taxonomy_dir"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done for now. Remember to double check everything! Bye bye."

