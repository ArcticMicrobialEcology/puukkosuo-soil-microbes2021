#!/bin/bash
#SBATCH --job-name=R_estimate_kegg_module_completeness
#SBATCH --error=14_read_ds_kegg_module_completeness_%A_%a_err.txt
#SBATCH --output=14_read_ds_kegg_module_completeness_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time 4:00:00
#SBATCH --account project_2005827

# a script to estimate KEGG module completeness for each sample seperately in the KEGG KO datasets 

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# this is the directory where all the required external R packages will be installed
r_package_dir="$software_base_dir/project_rpackages_430"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/read_based/downstream"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# the directory where the anvio KEGG database is installed
kegg_data_dir="$database_base_dir/KEGG" 

# the result directory for the LMM analysis
lmm_res_dir="$project_root/downstream/lmm_analysis"

# where the KO groups from the samples in the MG and MT data are saved as enzymes text within data and sample-specific subdirectories
base_dir="$project_root/downstream/kegg_modules"

# create the results directory
mkdir -p "$base_dir"

# load the R environment
module load r-env/432 # version
echo "[$(date -Is)] Exporting the KO groups present in samples in the KEGG KO datas..."
Rscript "$r_script_dir/9_READ_export_ko_groups_samples_for_anvio.R" "$mg_main_dir" "$mt_main_dir" "$base_dir"
echo "*******************************************************"
echo "*******************************************************"

# estimate module completeness in KEGG datas in metagenomics using an installed developer version of anvi'o
echo "[$(date -Is)] Estimating KEGG module completeness in the different samples in metagenomics data..."

export PATH="$software_base_dir/anvio_dev/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/sandbox:$PATH"
export PYTHONPATH="${PYTHONPATH:-}:$software_base_dir/anvio_git/anvio"

# some extra modules needed for this installation
module load gcc/11.3.0
module load biopythontools/11.3.0_3.10.6
module load hmmer/3.4

# go trough all the subdirectories
cd "$base_dir/kegg_mg_read"
DIRS=(*/)
for i in "${!DIRS[@]}"; do
	cd "${DIRS[$i]}"
	anvi-estimate-metabolism --enzymes-txt kegg_enzymes_for_anvio.txt --kegg-data-dir "$kegg_data_dir" --include-kos-not-in-kofam --exclude-dashed-reactions --output-modes modules,module_steps,module_paths --add-copy-number
	cd ..
done
echo "*******************************************************"
echo "*******************************************************"

# metatranscriptomics
echo "[$(date -Is)] Estimating KEGG module completeness in the different samples in metatranscriptomics data..."
cd "$base_dir/kegg_mt_read"
DIRS=(*/)
for i in "${!DIRS[@]}"; do
	cd "${DIRS[$i]}"
	anvi-estimate-metabolism --enzymes-txt kegg_enzymes_for_anvio.txt --kegg-data-dir "$kegg_data_dir" --include-kos-not-in-kofam --exclude-dashed-reactions --output-modes modules,module_steps,module_paths --add-copy-number
	cd ..
done
echo "*******************************************************"
echo "*******************************************************"

echo "[$(date -Is)] All done for now. Remember to double check everything! Bye bye."
