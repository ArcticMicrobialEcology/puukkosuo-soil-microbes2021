#!/bin/bash
#SBATCH --job-name=mag_processing_gtdbktk
#SBATCH --error=7_mg_mag_processing_gtdbktk_%A_%a_err.txt
#SBATCH --output=7_mg_mag_processing_gtdbktk_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH --time 71:00:00
#SBATCH --account project_2007998
#SBATCH --gres=nvme:100
#SBATCH --array 1-2

# a small script to process the metabat2 binned and manually refined MAGs and taxononomically annotate them with gtdb-tk
# requires manual refinement of the discovered MAGs resulting from the previous script to have been performed and 
# stored in the anvio contigs db with the name metabat2_refined

# define the needed directories
# this is the project root directory under which all the subdirectories, results etc. are created.
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="" # e.g. "/projappl/project_number" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# the list of treatment groups for which the coassembly is performed separately - prepared by the 1_SETUP script
all_conditions_file="$project_root/metadata/Chosen_Conditions.txt"

# define which coassembly is to be processed with the used array
CONDITION="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$all_conditions_file")"

# directories and variables for the project
main_name="$CONDITION"
main_dir="$mg_main_dir/mag_based/$CONDITION"

# using an installed developer version of anvio
export PATH="$software_base_dir/anvio_dev/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/sandbox:$PATH"
export PYTHONPATH="${PYTHONPATH:-}:$software_base_dir/anvio_git/anvio"

# some extra modules needed for this installation
module load gcc/11.3.0
module load biopythontools/11.3.0_3.10.6
module load hmmer/3.4

# start processing
echo "Processing MAGs in the $CONDITION condition..."
cd "$main_dir"

# call MAGs; summarize the refined bins and MAGs - by using the MIMAG standard =50% complete and<10% redundant
# rename the bins also
echo "[$(date -Is)] Copying the refined Metabat2 bins with anvi-rename-bins and calling MAGs as at least 50% complete and max 10% redundant into a new collectin called metabat2_refined_mags..."
anvi-rename-bins \
    --contigs-db "$main_dir/${main_name}_contigs.db" \
    --profile-db SAMPLES-MERGED/PROFILE.db \
    --collection-to-read metabat2_refined \
    --collection-to-write metabat2_refined_mags \
	--prefix metabat2_refined_mags \
	--call-MAGs \
	--min-completion-for-MAG 50 \
	--max-redundancy-for-MAG 10 \
	--exclude-bins \
	--report-file MAGs/metabat2/rename_bins_refined_mags.txt

# summarize the metabat2 MAGs
echo "[$(date -Is)] Summarizing the Metabat2 MAGs bins with anvi-summarize..."
anvi-summarize \
	--contigs-db "$main_dir/${main_name}_contigs.db" \
	--pan-or-profile-db SAMPLES-MERGED/PROFILE.db \
	--collection-name metabat2_refined_mags \
	--output-dir MAGs/metabat2/summary_refined_mags
echo "Done. Summary of the refined MAGs are in the condition subdirectory: MAGs/metabat2/summary_refined_mags."
echo "*******************************************************"
echo "*******************************************************"

# copy the bin fastas into one directory for gtdb-tk
cd MAGs/metabat2/summary_refined_mags
mkdir -p fasta
cd bin_by_bin
BINS=(*/)
for i in "${!BINS[@]}"; do
    cd "${BINS[$i]}"
    cp *contigs.fa ../../fasta
    cd ..
done

# run GTDB-Tk to get the taxonomy for the refined MAGs
echo "[$(date -Is)] Getting taxonomy for the refined MAGs with GTDB-Tk..."

# make gtdbtk available in the path
export PATH="$software_base_dir/gtdbtk/bin:$PATH"

# export gtdbtk database into PATH as well
export GTDBTK_DATA_PATH="$database_base_dir/gtdbtk_data/release220"

cd "$main_dir"
cd MAGs/metabat2/summary_refined_mags
mkdir -p tmp
gtdbtk classify_wf -x fa --genome_dir fasta --out_dir gtdbtk --cpus "$SLURM_CPUS_PER_TASK" --tmpdir tmp --mash_db "$GTDBTK_DATA_PATH/mash_db"
rm -rf tmp
echo "*******************************************************"
echo "*******************************************************"
