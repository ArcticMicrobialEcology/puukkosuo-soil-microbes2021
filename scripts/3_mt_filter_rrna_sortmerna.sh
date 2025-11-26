#!/bin/bash
#SBATCH --job-name=sortmerna
#SBATCH --error=3_mt_sortmerna_%A_%a_err.txt
#SBATCH --output=3_mt_sortmerna_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=320G
#SBATCH --time 71:00:00
#SBATCH --account project_2009164
#SBATCH --array 1-36

# a script to trim the metatranscriptomic reads for ribosomal RNA (rRNA)

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# sample name file for the metatranscriptomics samples - prepared by the 1_SETUP script
mt_sample_name_file="$project_root/metadata/MT_Sample_Names.txt"

# the location of sortmerna install
sortmerna_dir="$software_base_dir/sortmerna/bin"

# the location of the installed sortmerna database to be used
sortmerna_db="$database_base_dir/rRNA_databases_v4/smr_v4.3_sensitive_db.fasta" 

echo "[$(date -Is)] Running SortMeRNA and filtering the metatranscriptomics reads for ribosomalic RNA..."
# define the sample to be processed, run as an slurm array job
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$mt_sample_name_file")"

# make sortmerna available
export "PATH=$sortmerna_dir:$PATH"

# create the result and working directories
mkdir -p "$mt_main_dir/rrna_filtered_trimmed_fastq"
mkdir -p "$mt_main_dir/sortmerna_workdir"

cd "$mt_main_dir/sortmerna_workdir"
mkdir "$SAMPLE"

# run sortmerna
sortmerna -ref "$sortmerna_db" \
-reads "$mt_main_dir/cutadapt_trimmed_fastq/${SAMPLE}_R1_trimmed.fastq" \
-reads "$mt_main_dir/cutadapt_trimmed_fastq/${SAMPLE}_R2_trimmed.fastq" \
-other "$mt_main_dir/rrna_filtered_trimmed_fastq/$SAMPLE" \
--workdir "$mt_main_dir/sortmerna_workdir/$SAMPLE" \
--fastx --paired_in --out2 --threads "$SLURM_CPUS_PER_TASK"
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done. Bye
	

