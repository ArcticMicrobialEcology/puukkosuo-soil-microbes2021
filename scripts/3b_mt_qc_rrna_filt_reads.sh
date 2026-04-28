#!/bin/bash
#SBATCH --job-name=fastqc_cutadapt
#SBATCH --error=3b_qc_trim_rrna_filt_%A_err.txt
#SBATCH --output=3b_qc_trim_rrna_filt_%A_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --account=your_project_account

# a script to QC th rRNA filtered metatranscriptomics sequencing files

# define the needed directories
# this is the project root directory under which all the subdirectories, results etc. are created.
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# the directory that contains the rRNA filtered metatranscriptomis data
mt_raw_data_dir="$project_root/metatranscriptomics/rrna_filtered_trimmed_fastq"

# sample name file for the metatranscriptomics samples - prepared by the 1_SETUP script
mt_sample_name_file="$project_root/metadata/MT_Sample_Names.txt"

# perform fastqc for the fastq files
module load biokit/11.3.0 # fastqc 0.11.9 is used

echo "[$(date -Is)] Running Fastqc and Multiqc for the metatranscriptomics rRNA filtered fastq files..."

# create the directory
mkdir -p "$mt_main_dir/fastqc/cutadapt_trimmed/rrna_filtered" 
cd "$mt_raw_data_dir"
for j in * 
do
	fastqc --noextract -t "$SLURM_CPUS_PER_TASK" -o "$mt_main_dir/fastqc/cutadapt_trimmed/rrna_filtered/" "$j"
done

# perform multiqc for the fastqc results
module purge
module load multiqc/1.19
cd "$mt_main_dir/fastqc/cutadapt_trimmed/rrna_filtered"
multiqc .
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] Metatranscriptomics samples processed. Everything done. Bye for now."
