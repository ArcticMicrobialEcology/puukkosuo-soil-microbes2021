#!/bin/bash
#SBATCH --job-name=phyloflash
#SBATCH --error=5_mt_read_phyloflash_%A_%a_err.txt
#SBATCH --output=5_mt_read_phyloflash_%A_%a_out.txt
#SBATCH --partition hugemem
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 30
#SBATCH --mem 800G
#SBATCH --time 71:00:00
#SBATCH --account project_2009164
#SBATCH --array 1-36

# a script to determine the taxonomy in the metagenomics samples using phyloflash
# running phyloflash for the metatranscriptomics samples with lots of rRNA may take a long time and require a lot of memory

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# here are the trimmed fastq sequences files for the samples
fastq_dir="$mt_main_dir/cutadapt_trimmed_fastq" 

# the location of phyloflash install, version v3.4.2 used
phyloflash_dir="$software_base_dir/phyloflash/bin" 

# # the location of the installed phyloflash database to be used
phyloflash_db="$database_base_dir/phyloflash_databases/138.1" 

# sample name file for the metagenomics samples - prepared by the 1_SETUP script
mt_sample_name_file="$project_root/metadata/MT_Sample_Names.txt"

echo "[$(date -Is)] Running phyloFlash and determining the taxonomic composition..."

# define the sample to be processed, run as an slurm array job
SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$mt_sample_name_file")"

# make phyloflash available
export "PATH=$phyloflash_dir:$PATH"

# make the result directories
mkdir -p "$mt_main_dir/phyloflash"
cd "$mt_main_dir/phyloflash"

mkdir -p "csv_files"
mkdir -p "$SAMPLE"
cd "$SAMPLE"

# run phyloFlash using similar parameters to the recommended -almost_everything parameter
# -almost everything corresponds to  -poscov -treemap -zip -log. We don't need treemap or spades assembly which takes a long time. 
phyloFlash.pl -dbhome "$phyloflash_db" \
-lib "$SAMPLE" -CPUs "$SLURM_CPUS_PER_TASK" \
-read1 "$fastq_dir/${SAMPLE}_R1_trimmed.fastq" \
-read2 "$fastq_dir/${SAMPLE}_R2_trimmed.fastq" \
-skip_spades -zip -log -poscov

# unpack and copy the relevant output phyloFlash files for downstream processing
tar xzf "$SAMPLE.phyloFlash.tar.gz"  --wildcards --no-anchored "*NTU*.csv"
cp "*NTU*.csv" "$mt_main_dir/phyloflash/csv_files"

echo "*******************************************************"
echo "*******************************************************"
echo "[$(date -Is)] All done. Bye bye"
