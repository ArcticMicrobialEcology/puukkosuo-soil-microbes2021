#!/bin/bash
#SBATCH --job-name=mag_processing_dereplication_coverm
#SBATCH --error=8_mg_mag_dereplication_coverm_final_mags_%A_%a_err.txt
#SBATCH --output=8_mg_mag_dereplication_coverm_final_mags_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --time 48:00:00
#SBATCH --account project_2009164
#SBATCH --gres=nvme:50

# define the needed directories
# this is the project root directory under which all the subdirectories, results etc. are created.
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="" # e.g. "/projappl/project_number" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/mag_based"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# here are the trimmed fastq sequences files for the samples
fastq_dir="$mg_main_dir/cutadapt_trimmed_fastq" 

# the list of treatment groups for which the coassembly is performed separately - prepared by the 1_SETUP script
all_conditions_file="$project_root/metadata/Chosen_Conditions.txt"

# save the conditions for the coassemblies in a variable
mapfile -t CONDITIONS < "$all_conditions_file"

# save the internal genomes file in the MAG directory of the first condition, indexing starts from 0
main_dir="$mg_main_dir/mag_based/${CONDITIONS[0]}"

# the result directory for the contig-based analysis
contig_res_dir="$project_root/downstream/contig_based"

# dereplicate the MAGs from the different conditions
echo "[$(date -Is)] Processing the discovered Metabat2 and manually refined bins into dereplicated annotated MAGs..."

# create internal genomes text file for the refined MAGs from the different coasssemblies for Anvio
module load r-env/432 
Rscript "$r_script_dir/2_MAG_create_internal_genomes_refined_MAGs_for_anvio.R" "$main_dir/MAGs/metabat2" "$mg_main_dir/mag_based/${CONDITIONS[0]}/MAGs/metabat2/summary_refined_mags/bin_by_bin" "$mg_main_dir/mag_based/${CONDITIONS[1]}/MAGs/metabat2/summary_refined_mags/bin_by_bin"

# using an installed developer version of anvio for processing
export PATH="$software_base_dir/anvio_dev/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/sandbox:$PATH"
export PYTHONPATH="${PYTHONPATH:-}:$software_base_dir/anvio_git/anvio"

# some extra modules needed for this installation
module load gcc/11.3.0
module load biopythontools/11.3.0_3.10.6
module load hmmer/3.4

# make fastani available in the path
export PATH="$software_base_dir/fastani/bin:$PATH"

cd "$main_dir/MAGs/metabat2"
anvi-dereplicate-genomes --internal-genomes internal_genomes_refined_MAGs.txt \
                         --output-dir dereplicated_MAGs \
                         --program fastANI \
                         --similarity-threshold 0.99 \
                         --num-threads "$SLURM_CPUS_PER_TASK"
echo "*******************************************************"
echo "*******************************************************"

# process the dereplicated MAGs, get taxonomy, marker gene copies, etc.
echo "[$(date -Is)] Processing the dereplicated MAGs, getting taxonomy, marker gene copies, etc..."
mkdir -p "$mg_main_dir/mag_based/final_mags"
cd "$mg_main_dir/mag_based/final_mags"
Rscript "$r_script_dir/3_MAG_copy_dereplicated_MAGs_final_MAGs.R" "$main_dir/MAGs/metabat2/dereplicated_MAGs" "$mg_main_dir/mag_based/final_mags" "$mg_main_dir/mag_based"
echo "The final MAGs and related data can be found in the directory:" "$mg_main_dir/mag_based/final_mags"
echo "*******************************************************"
echo "*******************************************************"

# map the MAG genes to contig genes to use these for more accurate estimation of gene expression in MAGs
echo "[$(date -Is)] Mapping the MAG genes to the contig genes to assess gene expression in MAGs..."
Rscript "$r_script_dir/4_MAG_map_contig_genes_to_mag_genes.R" "$mg_main_dir/contig_based/reformat_report.txt" "$contig_res_dir" "$main_dir/MAGs/metabat2" "$main_dir/prodigal_genes" "$mg_main_dir/mag_based/${CONDITIONS[1]}/prodigal_genes"
echo "*******************************************************"
echo "*******************************************************"

# get the gene abundance and transcript expression datas summarized to marker gene level for the MAGs
echo "[$(date -Is)] Getting the gene transcript expression datas summarized to marker gene level for the MAGs..."
Rscript "$r_script_dir/5_MAG_get_mag_gene_transcript_count_datas.R" "$contig_res_dir" "$main_dir/MAGs/metabat2" "$mg_main_dir/mag_based/final_mags" "$mg_main_dir/kegg_diamond/Matrices_For_Downstream.RData" "$mt_main_dir/kegg_diamond/Matrices_For_Downstream.RData"
echo "*******************************************************"
echo "*******************************************************"

# finally get MAG abundance with CoverM
echo "[$(date -Is)] Get the relative abundance of MAGs in different samples with CoverM..."

# copy MAG fastas into one directory..."
mkdir mag_fasta
BINS=(cluster_*/)
for i in "${!BINS[@]}"; do
		cd "${BINS[$i]}"
		fasta=($(ls *.fa))
		old_name="${BINS[$i]}"
		new_name="${old_name%?}"
		cp "$fasta" "../mag_fasta/$new_name.fa"
		cd ..
done

cd "$mg_main_dir/mag_based/final_mags"
mkdir coverm

# make coverm available in the path
export PATH="$software_base_dir/coverm/bin:$PATH"
coverm genome \
  -1 "${fastq_dir}"/*_R1*.fastq \
  -2 "${fastq_dir}"/*_R1*.fastq \
  --genome-fasta-files mag_fasta/*.fa \
  --methods relative_abundance mean covered_fraction tpm \
  --min-read-aligned-percent 80 \
  --threads "$SLURM_CPUS_PER_TASK" \
  --output-file coverm/coverM_MG_final_metrics.tsv \
  --output-format dense
echo "*******************************************************"
echo "*******************************************************"


echo "All done, thank you. Please double-check everything..."
echo "*******************************************************"
echo "*******************************************************"
