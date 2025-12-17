#!/bin/bash
#SBATCH --job-name=anvio_metabat2_mag_binning
#SBATCH --error=6_mg_mag_binning_%A_%a_err.txt
#SBATCH --output=6_mg_mag_binning_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time 12:00:00
#SBATCH --account project_2007998
#SBATCH --gres=nvme:50
#SBATCH --array 1-2

# a small script to bin metagenome assembled genomes (MAGs) from the contigs from the different assemblies using metabat2 and anvio 
# (metabat2 needs to be installed and accessible by anvio as a driver)

# define the needed directories
# this is the project root directory under which all the subdirectories, results etc. are created.
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="" # e.g. "/projappl/project_number" 

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# the list of treatment groups for which the coassembly is performed separately - prepared by the 1_SETUP script
all_conditions_file="$project_root/metadata/Chosen_Conditions.txt"

# define which coassembly is to be processed with the used array
CONDITION="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$all_conditions_file")"

# directories and variables for the project
main_name="$CONDITION"
main_dir="$mg_main_dir/mag_based/$CONDITION"

# a workflow to bin the discovered contigs in the coassemblies to metagenome assembled genomes or MAGs metabat2
echo "[$(date -Is)] Binning the discovered contigs in the coassemblies to metagenome assembled genomes or MAGs with Metabat2..."

# process the chosen conditions and the related contigs
echo "[$(date -Is)] Processing $CONDITION condition..."

# using an installed developer version of anvio
export PATH="$software_base_dir/anvio_dev/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/sandbox:$PATH"
export PYTHONPATH="${PYTHONPATH:-}:$software_base_dir/anvio_git/anvio"

# some extra modules needed for this installation
module load gcc/11.3.0
module load biopythontools/11.3.0_3.10.6
module load hmmer/3.4

cd "$main_dir"
mkdir -p "$main_dir/MAGs"

# run metabat2 to bin the contigs
echo "[$(date -Is)] Running anvi-cluster-contigs with Metabat2 to bin the contigs into metagenome assembled genomes..."
anvi-cluster-contigs \
	--contigs-db "$main_dir/${main_name}_contigs.db" \
	--profile-db SAMPLES-MERGED/PROFILE.db \
	--collection-name metabat2 \
	--num-threads "$SLURM_CPUS_PER_TASK" \
	--driver metabat2 \
	--seed 123 \
	--log-file MAGs/metabat2_log.txt \
	--just-do-it
echo "*******************************************************"
echo "*******************************************************"

# renaming & copying the bins for refinement
mkdir -p "$main_dir/MAGs/metabat2"
echo "[$(date -Is)] Copying the Metabat2 bins for refining with anvi-rename-bins into a new collection named metabat2_refined with the prefix metabat2_refined..."
anvi-rename-bins \
    --contigs-db "$main_dir/${main_name}_contigs.db" \
    --profile-db SAMPLES-MERGED/PROFILE.db \
    --collection-to-read metabat2 \
    --collection-to-write metabat2_refined \
	--prefix metabat2_refined \
	--report-file MAGs/metabat2/rename_bins.txt
echo "*******************************************************"
echo "*******************************************************"

# summarize the metabat2 bins
echo "[$(date -Is)] Summarizing the Metabat2 bins with anvi-summarize..."
anvi-summarize \
    --contigs-db "$main_dir/${main_name}_contigs.db" \
    --pan-or-profile-db SAMPLES-MERGED/PROFILE.db \
    --collection-name metabat2_refined \
    --output-dir MAGs/metabat2/summary
echo "*******************************************************"
echo "*******************************************************"

echo "All done. At this point it would be we wise to check your bins and manually refine them before running the next step."
echo "See you later alligator."

# perform manual refinement of the metabat2 bins with anvi-refine
# for example:

# anvi-refine -c Grazed_contigs.db -p SAMPLES-MERGED/PROFILE.db --server-only --port-number 8080 --collection-name metabat2_refined --bin-id metabat2_refined_Bin_00001
# anvi-refine -c UnGrazed_contigs.db -p SAMPLES-MERGED/PROFILE.db --server-only --port-number 8080 --collection-name metabat2_refined --bin-id metabat2_refined_Bin_00004
