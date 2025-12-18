#!/bin/bash
#SBATCH --job-name=R_summarize_kegg
#SBATCH --error=7b_mt_read_R_summarize_kegg_%A_%a_err.txt
#SBATCH --output=7b_mt_read_R_summarize_kegg_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time 48:00:00
#SBATCH --account project_2009164

# a script to summarize the combined KEGG aligment results for the metatranscriptomics data to the KO level

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/read_based"

# load the R environment
module load r-env/432 # version
echo "[$(date -Is)] Summarizing the combined KEGG alignment hits for the paired reads to KO group level..."
Rscript "$r_script_dir/2b_READ_summarize_KEGG_aligment.R" "$mt_main_dir/kegg_diamond/sample_details" "$mt_main_dir/kegg_diamond/combined_hits_gene_aggregated" "$mt_main_dir/kegg_diamond"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done. Bye bye"	
