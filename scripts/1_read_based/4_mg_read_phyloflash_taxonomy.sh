#!/bin/bash
#SBATCH --job-name=phyloflash
#SBATCH --error=4_mg_read_phyloflash_%A_%a_err.txt
#SBATCH --output=4_mg_read_phyloflash_%A_%a_out.txt
#SBATCH --partition small
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 300G
#SBATCH --time 71:00:00
#SBATCH --account project_2009164

# a script to determine the taxonomy in the metagenomics samples using phyloflash

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

# the location of phyloflash install, version v3.4.2 used
phyloflash_dir="$software_base_dir/phyloflash/bin" 

# # the location of the installed phyloflash database to be used
phyloflash_db="$database_base_dir/phyloflash_databases/138.1" 

# sample name file for the metagenomics samples - prepared by the 1_SETUP script
mg_sample_name_file="$project_root/metadata/MG_Sample_Names.txt"

# map the samples into a variable
mapfile -t SAMPLES < "$mg_sample_name_file"

# as phyloflash is quite fast for the MG samples, run as a single job for MG
echo "[$(date -Is)] Running phyloFlash and determining the taxonomic composition..."

# make phyloflash available
export "PATH=$phyloflash_dir:$PATH"

# make the result directories
mkdir -p "$mg_main_dir/phyloflash"
cd "$mg_main_dir/phyloflash"
mkdir -p "csv_files"

# analyze the samples
for j in "${!SAMPLES[@]}"; do

	mkdir -p "${SAMPLES[$j]}"
	cd "${SAMPLES[$j]}"
	
	echo "[$(date -Is)] Running phyloFlash for sample ${SAMPLES[$j]}..."
	
	# run phyloFlash using similar parameters to the recommended -almost_everything parameter
	# -almost everything corresponds to  -poscov -treemap -zip -log. We don't need treemap or spades assembly which takes a long time. 
	# Could be simplified further as only the rRNA Silva database taxonomy information is used? 
	
	phyloFlash.pl -dbhome "$phyloflash_db" -lib "${SAMPLES[$j]}" \
	-CPUs "$SLURM_CPUS_PER_TASK" -read1 "$fastq_dir/${SAMPLES[$j]}_R1_trimmed.fastq" \
	-read2 "$fastq_dir/${SAMPLES[$j]}_R2_trimmed.fastq" \
	-skip_spades -zip -log -poscov
	
	# unpack and copy the relevant output phyloFlash files for downstream processing
	tar xzf "${SAMPLES[$j]}.phyloFlash.tar.gz" --wildcards --no-anchored "*NTU*.csv"
	cp "*NTU*.csv" "$mg_main_dir/phyloflash/csv_files"
	
	cd ..
	echo "*******************************************************"
	echo "*******************************************************"
done
echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done. Bye bye"