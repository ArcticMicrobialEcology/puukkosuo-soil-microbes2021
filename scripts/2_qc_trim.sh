#!/bin/bash
#SBATCH --job-name=fastqc_cutadapt
#SBATCH --error=2_qc_trim_%A_%a_err.txt
#SBATCH --output=2_qc_trim_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time 71:00:00
#SBATCH --account project_2009164

# a script to QC and trim the raw sequencing files for metagenomics and metatranscriptomics

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"
# create the directory
mkdir -p "$mg_main_dir" 

# the diretory that contains the raw metagenomics fastq files
mg_raw_data_dir="/scratch/project_2009164/2_OULANKA/1_RAW/metagenomics"

# sample name file for the metagenomics samples - prepared by the 1_SETUP script
mg_sample_name_file="$project_root/metadata/MG_Sample_Names.txt"

# perform fastqc for the raw data fastq files
module load biokit/11.3.0 # contains various programs, fastqc 0.11.9 is used
echo "[$(date -Is)] Running Fastqc and Multiqc for the metagenomics raw fastq files..."
mkdir -p "$mg_main_dir/fastqc" # create the directory
cd "$mg_raw_data_dir"
for j in * 
do
	fastqc --noextract -t "$SLURM_CPUS_PER_TASK" -o "$mg_main_dir/fastqc/" "$j"
done

# perform multiqc for the fastqc results
module purge
module load multiqc/1.19
cd "$mg_main_dir/fastqc"
multiqc .
echo "*******************************************************"
echo "*******************************************************"

# trim the raw fastq files with cutadapt
echo "[$(date -Is)] Running Cutadapt and trimming the metagenomics raw fastq files..."
# create the directory
mkdir -p "$mg_main_dir/cutadapt_trimmed_fastq" 

module purge
module load cutadapt/4.6

# map the samples into a variable
mapfile -t SAMPLES < "$mg_sample_name_file"

# trim each sample
for j in "${!SAMPLES[@]}"; do
	cutadapt "$mg_raw_data_dir/${SAMPLES[$j]}_R1_001.fastq.gz" \
			 "$mg_raw_data_dir/${SAMPLES[$j]}_R2_001.fastq.gz" \
			 -o "$mg_main_dir/cutadapt_trimmed_fastq/${SAMPLES[$j]}_R1_trimmed.fastq" \
			 -p "$mg_main_dir/cutadapt_trimmed_fastq/${SAMPLES[$j]}_R2_trimmed.fastq" \
			 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC \
			 -A CTGTCTCTTATACACATCTGACGCTGCCGACGA \
			 -m 50 \
			 -j "$SLURM_CPUS_PER_TASK" \
			 --nextseq-trim 20 > "$mg_main_dir/cutadapt_trimmed_fastq/${SAMPLES[$j]}_trim.log"
done
echo "*******************************************************"
echo "*******************************************************"

# perform fastqc for the cutadapt filtered files
module purge
module load biokit/11.3.0 # fastqc 0.11.9 is used

echo "[$(date -Is)] Running Fastqc and Multiqc for the trimmed metagenomics fastq files..."
# create the directory
mkdir -p "$mg_main_dir/fastqc/cutadapt_trimmed" 

cd "$mg_main_dir/cutadapt_trimmed_fastq"
for j in *.fastq 
do
	fastqc --noextract -t "$SLURM_CPUS_PER_TASK" -o "$mg_main_dir/fastqc/cutadapt_trimmed/" "$j"
done

# perform multiqc for the fastqc results
module purge
module load multiqc/1.19
cd "$mg_main_dir/fastqc/cutadapt_trimmed"
multiqc .
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] Metagenomics samples processed. Now processing metatranscriptomics samples."

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"
# create the directory
mkdir -p "$mt_main_dir" 

# the diretory that contains the raw metatranscriptomics fastq files
mt_raw_data_dir=/scratch/project_2009164/2_OULANKA/1_RAW/metatranscriptomics

# sample name file for the metatranscriptomics samples - prepared by the 1_SETUP script
mt_sample_name_file="$project_root/metadata/MT_Sample_Names.txt"

# perform fastqc for the raw data fastq files
module load biokit/11.3.0 # fastqc 0.11.9 is used

echo "[$(date -Is)] Running Fastqc and Multiqc for the metatranscriptomics raw fastq files..."
# create the directory
mkdir -p "$mt_main_dir/fastqc" 
cd "$mt_raw_data_dir"
for j in * 
do
	fastqc --noextract -t "$SLURM_CPUS_PER_TASK" -o "$mt_main_dir/fastqc/" "$j"
done

# perform multiqc for the fastqc results
module purge
module load multiqc/1.19
cd "$mt_main_dir/fastqc"
multiqc .
echo "*******************************************************"
echo "*******************************************************"

# trim the raw fastq files with cutadapt
echo "[$(date -Is)] Running Cutadapt and trimming the metatranscriptomics raw fastq files..."
# create the directory
mkdir -p "$mt_main_dir/cutadapt_trimmed_fastq" 

module purge
module load cutadapt/4.6

# map the samples into a variable
mapfile -t SAMPLES < "$mt_sample_name_file"

for j in "${!SAMPLES[@]}"; do
	cutadapt "$mt_raw_data_dir/${SAMPLES[$j]}_R1_001.fastq.gz" \
			 "$mt_raw_data_dir/${SAMPLES[$j]}_R2_001.fastq.gz" \
			 -o "$mt_main_dir/cutadapt_trimmed_fastq/${SAMPLES[$j]}_R1_trimmed.fastq" \
			 -p "$mt_main_dir/cutadapt_trimmed_fastq/${SAMPLES[$j]}_R2_trimmed.fastq" \
			 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
			 -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
			 -m 50 \
			 -j "$SLURM_CPUS_PER_TASK" \
			 --nextseq-trim 20 > "$mt_main_dir/cutadapt_trimmed_fastq/${SAMPLES[$j]}_trim.log"
done
echo "*******************************************************"
echo "*******************************************************"

# perform fastqc for the cutadapt filtered files
module purge
module load biokit/11.3.0 # fastqc 0.11.9 is used

echo "[$(date -Is)] Running Fastqc and Multiqc for the trimmed metatranscriptomics fastq files..."
# create the directory
mkdir -p "$mt_main_dir/fastqc/cutadapt_trimmed" 

cd "$mt_main_dir/cutadapt_trimmed_fastq"
for j in *.fastq 
do
	fastqc --noextract -t "$SLURM_CPUS_PER_TASK" -o "$mt_main_dir/fastqc/cutadapt_trimmed/" "$j"
done

# perform multiqc for the fastqc results
module purge
module load multiqc/1.19
cd "$mt_main_dir/fastqc/cutadapt_trimmed"
multiqc .
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] Metatranscriptomics samples processed. Everything done. Bye for now."
