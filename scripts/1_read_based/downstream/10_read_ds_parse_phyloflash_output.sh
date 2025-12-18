#!/bin/bash
#SBATCH --job-name=R_parse_phyloflash
#SBATCH --error=10_read_ds_parse_phyloflash_DA_%A_%a_err.txt
#SBATCH --output=10_read_ds_parse_phyloflash_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time 2:00:00
#SBATCH --account project_2005827

# a small script to combine and process the phyloflash runs for the different samples into a neat matrix and objects ready for downstream processing

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

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

# load the R environment
module load r-env/432 # version

# compose otu and taxonomy tables from phyloflash output - process both metagenomics and metatransciptomics
echo "[$(date -Is)] Compose OTU and taxonomy tables from phyloFlash output..."
Rscript "$r_script_dir/4_READ_parse_phyloflash_output.R" "$mg_phylo_dir/csv_files" "$mt_phylo_dir/csv_files" "$metadata_file"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done for now. Remember to double check everything! Bye bye."

