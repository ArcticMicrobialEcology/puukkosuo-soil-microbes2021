#!/bin/bash
#SBATCH --job-name=megahit_assembly
#SBATCH --error=4_mg_contig_megahit_coassembly_%A_%a_err.txt
#SBATCH --output=4_mg_contig_megahit_coassembly_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=300G
#SBATCH --time 71:00:00
#SBATCH --account project_2005827

# a script to perform assembly for the metagenomics reads into contigs using megahit
# perform coassembly separately for the samples in the grazed upslope and ungrazed downslope conditions / levels of the exclusion treatment

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# the list of treatment groups for which the coassembly is performed separately - prepared by the 1_SETUP script
all_conditions_file="$project_root/metadata/Chosen_Conditions.txt"
mapfile -t CONDITIONS < "$all_conditions_file"

# create the needed directories and result directories
mkdir -p "$mg_main_dir/cutadapt_trimmed_fastq/concatenated_fastq"
mkdir -p "$mg_main_dir/megahit_assemblies"

# use the trimmed fastq files
cd "$mg_main_dir/cutadapt_trimmed_fastq"

# load megahit
module load biokit/11.3.0 # megahit v1.2.9 is used

# loop through all the chosen conditions and make coassemblies
echo "[$(date -Is)] Concatenating the trimmed metagenomics fastq files related to the chosen conditions and composing coassemblies from them using Megahit..."
for j in "${!CONDITIONS[@]}"; do

	# concantenate sample fastqs into single fastq file for the condition
	echo "[$(date -Is)] Concatenating metagenomics fastq files related to the ${CONDITIONS[$j]} condition..."
	
	# sample name file for the metagenomics samples in the condition of interest - prepare
	sample_name_file="$project_root/metadata/MG_Sample_Names_${CONDITIONS[$j]}.txt"
	
	FILES_R1=()
	while IFS=$'\n' read -r line; do
	  FILES_R1+=( "${line}_R1_trimmed.fastq" )
	done < "$sample_name_file"
	
	FILES_R2=()
	while IFS=$'\n' read -r line; do
	  FILES_R2+=( "${line}_R2_trimmed.fastq" )
	done < "$sample_name_file"
	
	cat "${FILES_R1[@]}" > "$mg_main_dir/cutadapt_trimmed_fastq/concatenated_fastq/${CONDITIONS[$j]}_samples_R1.fastq"
	cat "${FILES_R2[@]}" > "$mg_main_dir/cutadapt_trimmed_fastq/concatenated_fastq/${CONDITIONS[$j]}_samples_R2.fastq"
	
	# run megahit and compose coassemblies
	echo "[$(date -Is)] Running Megahit for the samples related to the ${CONDITIONS[$j]} condition and coassembling reads into contigs..."
	
	megahit -1 "$mg_main_dir/cutadapt_trimmed_fastq/concatenated_fastq/${CONDITIONS[$j]}_samples_R1.fastq" \
        -2 "$mg_main_dir/cutadapt_trimmed_fastq/concatenated_fastq/${CONDITIONS[$j]}_samples_R2.fastq" \
        --out-dir "$mg_main_dir/megahit_assemblies/${CONDITIONS[$j]}" \
        --min-contig-len 1000 \
        --k-min 57 \
        --k-max 157 \
        --k-step 10 \
        --num-cpu-threads "$SLURM_CPUS_PER_TASK" &> "$mg_main_dir/megahit_assemblies/${CONDITIONS[$j]}_megahit_log.txt"
	
	echo "*******************************************************"
	echo "*******************************************************"
done
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All finished. Bye"
