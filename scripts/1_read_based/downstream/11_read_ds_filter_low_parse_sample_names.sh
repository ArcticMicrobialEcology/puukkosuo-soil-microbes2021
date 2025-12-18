#!/bin/bash
#SBATCH --job-name=R_filter_datas
#SBATCH --error=11_read_ds_filter_datas_%A_%a_err.txt
#SBATCH --output=11_read_ds_filter_datas_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time 2:00:00
#SBATCH --account project_2005827

# a script to filter all the composed datas for lowly expressed features similarly

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

# load the R environment
module load r-env/432 # version

# filter all the different functional datas and taxonomy datas for lowly expressed features, parse sample names and order according to metadata
echo "[$(date -Is)] Filter all the datasets for lowly expressed features, parse sample names and order according to metadata..."
Rscript "$r_script_dir/5_READ_filter_all_datas.R" "$mg_main_dir" "$mt_main_dir"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done for now. Remember to double check everything! Bye bye."


