#!/bin/bash
#SBATCH --job-name=R_perform_lmm_functional
#SBATCH --error=12_read_ds_lmm_func_%A_%a_err.txt
#SBATCH --output=12_read_ds_lmm_func_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time 48:00:00
#SBATCH --account project_2005827

# a script to perform differential abundance and expression analysis using linear mixed effects models for all the functional datas

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

# the result directory for the LMM analysis
lmm_res_dir="$project_root/downstream/lmm_analysis"

# create the results directory
mkdir -p "$lmm_res_dir"

# load the R environment
module load r-env/432 # version

# perform differential abundance and expression analysis using linear mixed effects models for the KEGG datas
echo "[$(date -Is)] Perform linear mixed effects model analysis for the KEGG KO and marker gene datasets..."
Rscript "$r_script_dir/6a_READ_perform_lmm_for_kegg_datas.R" "$mg_main_dir" "$mt_main_dir" "$lmm_res_dir"
echo "*******************************************************"
echo "*******************************************************"

# perform differential abundance and expression analysis using linear mixed effects models for the nitrogen and methane related genes
echo "[$(date -Is)] Perform linear mixed effects model analysis for the KEGG KO and marker gene datasets..."
Rscript "$r_script_dir/6b_READ_perform_lmm_for_metmarkdb_datas.R" "$mg_main_dir" "$mt_main_dir" "$lmm_res_dir"
echo "*******************************************************"
echo "*******************************************************"

echo "[$(date -Is)] All done for now. Remember to double check everything! Bye bye."



