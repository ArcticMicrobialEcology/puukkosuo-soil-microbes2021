#!/bin/bash
#SBATCH --job-name=diamond_metmarkdb
#SBATCH --error=9_mt_read_diamond_metmarkdb_%A_%a_err.txt
#SBATCH --output=9_mt_read_diamond_metmarkdb_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time 71:00:00
#SBATCH --account project_2005827

# a script to align the metatranscriptomics reads to the metabolic marker gene database, combine the results from the paired reads and summarize to marker gene level

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# here are the trimmed fastq files filtered for ribosomal RNA for the samples
fastq_dir="$mt_main_dir/rrna_filtered_trimmed_fastq" 

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/read_based"

# this database is 51 functional marker gene DIAMOND database (Greening Lab) -with eukaryotic proteins filtered- and related info
# prepared  by the 1_setup script
met_mark_diamond_db="$database_base_dir/metmarkdb/MetMarkDB_prok.dmnd"
met_markdb_database_dir="$database_base_dir/metmarkdb"

# sample name file for the metagenomics samples - prepared by the 1_SETUP script
mt_sample_name_file="$project_root/metadata/MT_Sample_Names.txt"

mapfile -t SAMPLES < "$mt_sample_name_file"

# make the result directories
mkdir -p "$mt_main_dir/metmarkdb_diamond/raw"
cd "$mt_main_dir/metmarkdb_diamond/raw"

# align for both read pairs using DIAMOND
echo "[$(date -Is)] Running DIAMOND and aligning the metatranscriptomic reads to the 51 functional metabolic marker gene prokaryotic database..."
module load diamond/2.1.6 # check version if necessary

# all samples in one slurm job
for j in "${!SAMPLES[@]}"; do
	
	echo "[$(date -Is)] Running DIAMOND and processing the sample: ${SAMPLES[$j]}"
	# R1
	diamond blastx --db "$met_mark_diamond_db" \
		 --query "$fastq_dir/${SAMPLES[$j]}_fwd.fq" \
		 --out "${SAMPLES[$j]}_R1_metmarkdb.txt" \
		 --query-cover 80 --max-target-seqs 1 --threads "$SLURM_CPUS_PER_TASK" \
		 --outfmt 6 qseqid stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp slen

	# R2
	diamond blastx --db "$met_mark_diamond_db" \
		 --query "$fastq_dir/${SAMPLES[$j]}_rev.fq" \
		 --out "${SAMPLES[$j]}_R2_metmarkdb.txt" \
		 --query-cover 80 --max-target-seqs 1 --threads "$SLURM_CPUS_PER_TASK" \
		 --outfmt 6 qseqid stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp slen
		 
	echo "*******************************************************"
	echo "*******************************************************"
done

# combine the hits from the paired reads using a custom R script
module purge

# load the R environment
module load r-env/432

echo "[$(date -Is)] Combining the MetMarkDB alignment hits for the paired reads and aggregating them to gene level..."
Rscript "$r_script_dir/3_READ_combine_MetMarkDB_alignment_paired_reads_summarize.R" "$mt_main_dir/kegg_diamond/sample_details" "$mt_main_dir/metmarkdb_diamond/raw" "$mt_main_dir/metmarkdb_diamond" "0"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done. Bye bye"

