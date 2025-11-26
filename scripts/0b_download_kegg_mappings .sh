#!/bin/bash
#SBATCH --job-name=download_kegg
#SBATCH --error=0b_download_kegg_%A_%a_err.txt
#SBATCH --output=0b_download_kegg_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time 72:00:00
#SBATCH --account project_2009164

# a script to update KEGG prokaryotic gene KO group mappings - not necessary to run, if some old version exists?

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts"

# this is the directory where the downloaded KEGG mappings will be saved, currently the metadata directory
metadata_dir="$project_root/metadata"

# update KEGG annotations
mkdir -p "$metadata_dir"
module load r-env/432
echo "[$(date -Is)] Downloading KEGG KO→gene mappings and annotations from KEGG using KEGG REST..."
Rscript "$r_script_dir/0b_SETUP_download_KEGG_prok_gene_KO_mappings.R" "$metadata_dir" 
echo "*******************************************************"
echo "*******************************************************"

echo "Everything done. See ya later alligator."