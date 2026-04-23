#!/bin/bash
#SBATCH --job-name=contig_mapping_stats
#SBATCH --error=7_contig_mapping_stats_%A_err.txt
#SBATCH --output=7_contig_mapping_stats_%A_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --account=project_2009164

# define the needed directories
# this is the project root directory under which all the subdirectories, results etc. are created.
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# here are the trimmed fastq sequences files for the samples
fastq_dir="$mg_main_dir/cutadapt_trimmed_fastq" 

# the list of treatment groups for which the coassembly is performed separately - prepared by the 1_SETUP script
all_conditions_file="$project_root/metadata/Chosen_Conditions.txt"
mapfile -t CONDITIONS < "$all_conditions_file"

# create the output directory
mkdir -p "$mg_main_dir/megahit_assemblies/bowtie_evaluation"

# run bowtie to see how the samples in the different coassemblies map to the generated contigs in those conditions
module load bowtie2/2.4.4

# loop through all the chosen conditions and map samples to the assembled contigs
echo "[$(date -Is)] Running bowtie and mapping the samples in different conditions to the assembled contigs..."

for j in "${!CONDITIONS[@]}"; do

	echo "[$(date -Is)] Mapping samples related to the ${CONDITIONS[$j]} condition..."
	
	# sample name file for the metagenomics samples in the condition of interest - prepare
	sample_name_file="$project_root/metadata/MG_Sample_Names_${CONDITIONS[$j]}.txt"
	
	# map the original reads to the coassembled contig
	mapfile -t SAMPLES < "$sample_name_file"
	
	# create the output directory
	mkdir -p "$mg_main_dir/megahit_assemblies/bowtie_evaluation/${CONDITIONS[$j]}"
	
	# index the reference genome to be used from the coassembled contig
	bowtie2-build "$mg_main_dir/megahit_assemblies/${CONDITIONS[$j]}/final.contigs.fa" "$mg_main_dir/megahit_assemblies/bowtie_evaluation/${CONDITIONS[$j]}" --threads "$SLURM_CPUS_PER_TASK"

	for i in "${!SAMPLES[@]}"; do
	   bowtie2 \
			-1 "$fastq_dir/${SAMPLES[$i]}_R1_trimmed.fastq" \
			-2 "$fastq_dir/${SAMPLES[$i]}_R2_trimmed.fastq" \
			-x "$mg_main_dir/megahit_assemblies/bowtie_evaluation/${CONDITIONS[$j]}" \
			-S "$mg_main_dir/megahit_assemblies/bowtie_evaluation/${CONDITIONS[$j]}/${SAMPLES[$i]}.sam" \
			--threads "$SLURM_CPUS_PER_TASK" \
			--no-unal &> "$mg_main_dir/megahit_assemblies/bowtie_evaluation/${CONDITIONS[$j]}_${SAMPLES[$i]}_bowtie_log.txt"
	done 
	
	echo "*******************************************************"
	echo "*******************************************************"
done

echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done for now. See you later."
