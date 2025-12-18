#!/bin/bash
#SBATCH --job-name=diamond_kegg
#SBATCH --error=6a_mg_read_diamond_kegg_%A_%a_err.txt
#SBATCH --output=6a_mg_read_diamond_kegg_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time 71:00:00
#SBATCH --account project_2005827
#SBATCH --array 1-36

# a script to align the metagenomic reads to the KEGG prokaryotic database and combine the results from the paired reads

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# here are the trimmed fastq sequences files for the samples
fastq_dir="$mg_main_dir/cutadapt_trimmed_fastq" 

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/read_based"

# this database is the KEGG prokaryotic gene amino acid sequences translated into an DIAMOND database - Acquired from Jenni Hultman
# copied by the setup script 0a
kegg_diamond_db="$database_base_dir/KEGG_PROK/PROKARYOTES.dmnd"

# this is the mapping file between KEGG genes and KO numbers -  needs to preexist - downloaded and prepared by the setup script 0b (or copied from somewhere)
kegg_gene_ko_mapping_file="$project_root/metadata/2025_02_28_KO_Mapping_Genes_Prokaryotes.txt" 

# sample name file for the metagenomics samples - prepared by the 1_SETUP script
mg_sample_name_file="$project_root/metadata/MG_Sample_Names.txt"

# define the sample to be processed, run as an slurm array job
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$mg_sample_name_file")"

# make the result directories
mkdir -p "$mg_main_dir/kegg_diamond/raw"
cd "$mg_main_dir/kegg_diamond/raw"

# align for both read pairs using DIAMOND
module load diamond/2.1.6

echo "[$(date -Is)] Running DIAMOND and aligning the metagenomic reads to the KEGG prokaryotic database..."
#R1
diamond blastx --db "$kegg_diamond_db" \
	--query "$fastq_dir/${SAMPLE}_R1_trimmed.fastq" \
	--out "${SAMPLE}_kegg_R1.txt" \
	--query-cover 80 --max-target-seqs 1 --threads "$SLURM_CPUS_PER_TASK" \
	--outfmt 6 qseqid stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp slen

#R2
diamond blastx --db "$kegg_diamond_db" \
	--query "$fastq_dir/${SAMPLE}_R2_trimmed.fastq" \
	--out "${SAMPLE}_kegg_R2.txt" \
	--query-cover 80 --max-target-seqs 1 --threads "$SLURM_CPUS_PER_TASK" \
	--outfmt 6 qseqid stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp slen

echo "*******************************************************"
echo "*******************************************************"	

# combine the hits from the paired reads using a custom R script
module purge

# load the R environment
module load r-env/432
echo "Combining the KEGG alignment hits for the paired reads and aggregating them to gene level..."
Rscript "$r_script_dir/2a_READ_combine_KEGG_alignment_paired_reads.R" "${SAMPLE}_kegg_R1.txt" "${SAMPLE}_kegg_R2.txt" "$kegg_gene_ko_mapping_file" "0" "$mg_main_dir/kegg_diamond"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done. Bye bye"
	
