#!/bin/bash
#SBATCH --job-name=anvio_contig_workflow
#SBATCH --error=5_mg_contig_anvio_complete_%A_%a_err.txt
#SBATCH --output=5_mg_contig_anvio_complete_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=320G
#SBATCH --time 72:00:00
#SBATCH --account project_2009164
#SBATCH --gres=nvme:200

# a complete workflow to process the contigs assembled from metagenomics data using anvio
# this workflow concatenates the contigs from different coassemblies into one and processes them all together

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# the database directory where the needed databases are installed
database_base_dir="$project_root/databases"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts/contig_based"

# the main folder, where all the metagenomics analysis steps and results will be stored
mg_main_dir="$project_root/metagenomics"

# here are the trimmed fastq sequences files for the samples
fastq_dir="$mg_main_dir/cutadapt_trimmed_fastq" 

# directories for various software to be used
feature_counts_dir="$software_base_dir/subread/bin" # the location of featureCounts install

# directories and variables for the project
main_name="contig_based"
main_dir="$mg_main_dir/$main_name"
anvio_project_name="Oulanka_$main_name"

# the list of treatment groups for which the coassembly is performed separately - prepared by the 1_SETUP script
all_conditions_file="$project_root/metadata/Chosen_Conditions.txt"

# sample names
all_sample_name_file="$project_root/metadata/MG_Sample_Names.txt"

# to be (created) gene subdirectory
gene_dir="prodigal_genes"

# directories and variables related to metatranscriptomics data
# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# sample name file for the metatranscriptomics samples - prepared by the 1_SETUP script
mt_sample_name_file="$project_root/metadata/MT_Sample_Names.txt"

# here are the trimmed fastq files filtered for ribosomal RNA for the samples
mt_fastq_dir="$mt_main_dir/rrna_filtered_trimmed_fastq" 

# this database is 51 functional marker gene DIAMOND database (Greening Lab) -with eukaryotic proteins filtered- and related info
# prepared  by the 1_setup script
met_mark_diamond_db="$database_base_dir/metmarkdb/MetMarkDB_prok.dmnd"
met_mark_annot_file="$database_base_dir/metmarkdb/Gene_metadata.RData"

# directory for anvio database(s)
scg_gtdb_data_dir="$database_base_dir/SCG_GTDB"

# using an installed developer version of anvio
export PATH="$software_base_dir/anvio_dev/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/sandbox:$PATH"
export PYTHONPATH="${PYTHONPATH:-}:$software_base_dir/anvio_git/anvio"

# some extra modules needed for this installation
module load gcc/11.3.0
module load biopythontools/11.3.0_3.10.6
module load hmmer/3.4

mkdir -p "$main_dir"
mkdir -p "$main_dir/bam_files"
cd "$main_dir"

# first, reformat the assembly contig fasta for the different coassemblies / conditions
echo "[$(date -Is)] Running anvi-script-reformat-fasta, reformatting the coassembly contig fastas..."

# format the assembled contigs fasta from both conditions 
mapfile -t CONDITIONS < "$all_conditions_file"

for j in "${!CONDITIONS[@]}"; do
	echo "Processing coassembly from condition: ${CONDITIONS[$j]}"
	
	contig_fasta="$mg_main_dir/megahit_assemblies/${CONDITIONS[$j]}/final.contigs.fa"
	
	anvi-script-reformat-fasta \
		--min-len 2000 \
		--simplify-names \
		-o "$main_dir/${CONDITIONS[$j]}_contigs.fasta" \
		--report-file "$main_dir/${CONDITIONS[$j]}_reformat_report.txt" \
		"$contig_fasta"
	echo "*******************************************************"
	echo "*******************************************************"
done 

echo "Done reformatting the coassembly contig fastas..."
echo "*******************************************************"
echo "*******************************************************"

# concatenate the contigs from the different coasseblies into one
# do this because we want to process all the contigs simultaneuosly downstream. To achieve this, we need to get gene counts for the
# called contig genes by Anvio - that is we want to align our metagenomic reads againts the called genes. We want each gene to be aligned only once. Since
# different coassemblies may contain very similar or the same genes, if we were to align the reads separately for the coassemblies and then count 
# gene hits, we would most likely use one read more than once (which is something we don't want to do). In the case of multiple identical genes where
# a read might align to, typically it is randomly assigned to one of these genes. However, in our case this will not produce a problem, as we are
# composing the data to functional unit level (e.g. KO, marker gene) and the same genes (or very similar) from the different coassemblies will have
# the same functionality and then each read is correctly counted only once. If then again we were to align the reads separately for the coassemblies,
# summarize to functional unit level and then try to combine the functional datas for the different coassemblies (e.g. by summing) - we would most likely have problems
# as we would have counted many reads for both coassemblies and would then thus easily overestimate the abundances for many functional units.

echo "[$(date -Is)] Concatenating contig databases into one..."

# add coassembly name for each fasta header - this way we can keep track of the contigs more easily, from which condition they came from.
while read line
do
	sed '/^>/s/$/_'${line}'/' "$main_dir/${line}_contigs.fasta" > "$main_dir/${line}_headermod_contigs.fasta"
	
done < "$all_conditions_file"

# concatenate the modified fasta files into one
> "$main_dir/${main_name}_contigs_concatenated.fasta"
while read -r line; do
    cat "$main_dir/${line}_headermod_contigs.fasta" >> "$main_dir/${main_name}_contigs_concatenated.fasta"
done < "$all_conditions_file"
echo "*******************************************************"
echo "*******************************************************"

# reformat the assembly contig fasta - to let anvio create unique names for the contigs
# enforce still the same limits for the contigs (should already be the case after running anvi-script-reformat-fasta separately for the different coasseblies 
contig_fasta="$main_dir/${main_name}_contigs_concatenated.fasta"
echo "[$(date -Is)] Running anvi-script-reformat-fasta, reformatting the assembly contig fasta..."
anvi-script-reformat-fasta \
    --min-len 2000 \
    --simplify-names \
    -o "$main_dir/${main_name}_contigs.fasta" \
    --report-file "$main_dir/reformat_report.txt" \
    "$contig_fasta"
echo "*******************************************************"
echo "*******************************************************"

# generate the contigs DB
echo "[$(date -Is)] Running anvi-gen-contigs-database, generating the contigs DB..."
anvi-gen-contigs-database \
    -f "$main_dir/${main_name}_contigs.fasta" \
    -o "$main_dir/${main_name}_contigs.db" \
    -n "$anvio_project_name" \
    -T "$SLURM_CPUS_PER_TASK"
echo "*******************************************************"
echo "*******************************************************"

echo "[$(date -Is)] Set up some missing Anvi'o databases using DIAMOND..."
module load diamond/2.1.6
# setup the single-copy-core-gene taxonomy data
echo "Setting up the scg taxonomy database..."
anvi-setup-scg-taxonomy \
	-T "$SLURM_CPUS_PER_TASK" \
	--scgs-taxonomy-data-dir "$scg_gtdb_data_dir"
echo "*******************************************************"
echo "*******************************************************"

# annotate the single-copy-core-genes using all the default sources in anvio
echo "[$(date -Is)] Running anvi-run-hmms and annotating the contigs and contig genes with all the anvio default HMM sources..."
anvi-run-hmms \
    -c "$main_dir/${main_name}_contigs.db" \
    -T "$SLURM_CPUS_PER_TASK"
echo "*******************************************************"
echo "*******************************************************"

# annotate single-copy core genes with taxonomy
echo "[$(date -Is)] Running anvi-run-scg-taxonomy, annotating single-copy core genes..."
anvi-run-scg-taxonomy \
    -c "$main_dir/${main_name}_contigs.db" \
	--scgs-taxonomy-data-dir "$scg_gtdb_data_dir" \
    -T "$SLURM_CPUS_PER_TASK"
echo "*******************************************************"
echo "*******************************************************"

# export the called genes by prodigal in anvio
mkdir -p "$main_dir/$gene_dir"
echo "[$(date -Is)] Running anvi-export-gene-calls, exporting called prodigal genes with AA sequences..."
anvi-export-gene-calls \
	-c "$main_dir/${main_name}_contigs.db" \
	--gene-caller prodigal \
	-o "$main_dir/$gene_dir/prodigal_genes.txt"	
echo "*******************************************************"
echo "*******************************************************"

# run again without sequence reporting for downstream scripts
anvi-export-gene-calls \
	-c "$main_dir/${main_name}_contigs.db" \
	--gene-caller prodigal \
	--skip-sequence-reporting \
	-o "$main_dir/$gene_dir/prodigal_genes_no_aaseq.txt"	
echo "*******************************************************"
echo "*******************************************************"

# export gene sequences into fasta
echo "[$(date -Is)] Running anvi-get-sequences-for-gene-calls, exporting ALL called gene sequences into fasta..."
anvi-get-sequences-for-gene-calls \
	-c "$main_dir/${main_name}_contigs.db" \
	-o "$main_dir/$gene_dir/gene_sequences_all.fasta"	
echo "*******************************************************"
echo "*******************************************************"

# export gene sequences as GFF
echo "[$(date -Is)] Running anvi-get-sequences-for-gene-calls, exporting ALL called gene sequences into GFF..."
anvi-get-sequences-for-gene-calls \
	-c "$main_dir/${main_name}_contigs.db" \
	--export-gff3 \
	-o "$main_dir/$gene_dir/gene_sequences_all.gff"	
echo "*******************************************************"
echo "*******************************************************"

# use DIAMOND to blast exported gene sequences against the functional metabolic marker 51 database
echo "[$(date -Is)] Running DIAMOND to blast the exported gene sequences against the functional metabolic marker 51 database..."
module load diamond/2.1.6 
diamond blastx --db "$met_mark_diamond_db" \
	--query "$main_dir/$gene_dir/gene_sequences_all.fasta" \
	--out "$main_dir/$gene_dir/prodigal_gene_metmarkdb_hits.txt" \
	--query-cover 80 --max-target-seqs 1 --threads "$SLURM_CPUS_PER_TASK" \
	--outfmt 6 qseqid stitle pident length qstart qend sstart send evalue bitscore qcovhsp scovhsp slen
echo "*******************************************************"
echo "*******************************************************"

# run an R script to parse the DIAMOND alignment results into the correct format for importing into anvio
module load r-env/432 
Rscript "$r_script_dir/2_CONTIG_parse_metmarkdb.R" "$main_dir/$gene_dir" "$met_mark_annot_file" 0

# import the metabolic marker gene hits into the anvio contigs database
echo "[$(date -Is)] Running anvi-import-functions, importing the metabolic marker gene hits into the anvio contigs database..."
anvi-import-functions \
	-c "$main_dir/${main_name}_contigs.db" \
	-i "$main_dir/$gene_dir/MetMarkDB_for_AnvioFunc.txt"
echo "*******************************************************"
echo "*******************************************************"

# run bowtie 2 and align the samples to the contigs
echo "[$(date -Is)] Running Bowtie2 and aligning the metagenome samples to the contigs..."
module load bowtie2/2.4.4
module load samtools/1.18

# index the reference genome to be used from the coassembled contigs
bowtie2-build --threads "$SLURM_CPUS_PER_TASK" "$main_dir/${main_name}_contigs.fasta" "$main_dir/$main_name"

# align the samples to the coassembled contigs
mapfile -t SAMPLES < "$all_sample_name_file"

for j in "${!SAMPLES[@]}"; do
   bowtie2 \
        -1 "$fastq_dir/${SAMPLES[$j]}_R1_trimmed.fastq" \
        -2 "$fastq_dir/${SAMPLES[$j]}_R2_trimmed.fastq" \
        -x "$main_dir/$main_name" \
        -S "$main_dir/${SAMPLES[$j]}.sam" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --no-unal

    samtools view -@ "$SLURM_CPUS_PER_TASK" -F 4 -bS "$main_dir/${SAMPLES[$j]}.sam" |\
        samtools sort -@ "$SLURM_CPUS_PER_TASK" > "$main_dir/bam_files/${SAMPLES[$j]}.bam"
    
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$main_dir/bam_files/${SAMPLES[$j]}.bam"
    rm 	"$main_dir/${SAMPLES[$j]}.sam"
done 
echo "*******************************************************"
echo "*******************************************************"

# create sample profiles in anvio
echo "[$(date -Is)] Running anvi-profile and creating sample profiles in anvio..."

# sample name (defined by the -S parameter) needs to be corrected for anvio. It does not accept sample names 
# starting with numbers (add sample_ as prefix for each sample) or sample names including special characters 
# (in this case - which is replaced by _ ).

for j in "${!SAMPLES[@]}"; do 
    anvi-profile \
        -i "$main_dir/bam_files/${SAMPLES[$j]}.bam" \
        -c "$main_dir/${main_name}_contigs.db" \
        -S sample_"${SAMPLES[$j]//[-]/_}" \
		--min-contig-length 2000 \
        -o "$main_dir/${SAMPLES[$j]}_PROFILE" \
        -T "$SLURM_CPUS_PER_TASK"
done
echo "*******************************************************"
echo "*******************************************************"

# merge sample profiles for anvio
echo "[$(date -Is)] Running anvi-merge and merging sample profiles for anvio..."
anvi-merge \
    -o "$main_dir/SAMPLES-MERGED" \
    -c "$main_dir/${main_name}_contigs.db" \
    "$main_dir"/*_PROFILE/PROFILE.db 
echo "*******************************************************"
echo "*******************************************************"	

# export the metabolic marker gene annotations for downstream analysis
echo "[$(date -Is)] Exporting metabolic marker gene database annotations..."
anvi-export-functions \
-c "$main_dir/${main_name}_contigs.db" \
-o "$main_dir/$gene_dir/exported_functions_metmarkdb.txt" \
--annotation-sources MetMarkDB
echo "*******************************************************"
echo "*******************************************************"	

# estimate taxonomy for the contigs with anvio
echo "[$(date -Is)] Running anvi-estimate-scg-taxonomy, estimating taxonomy for the contigs and exporting those from anvio..."
anvi-estimate-scg-taxonomy \
	-c "$main_dir/${main_name}_contigs.db" \
	--metagenome-mode \
	-o "$main_dir/anvio_scg_estimated_contigs_taxonomy.txt" 
echo "*******************************************************"
echo "*******************************************************"	

# finally, get gene counts to prodigal genes for the different samples with featureCounts
echo "[$(date -Is)] Get counts for all the genes in the contigs for each sample with featureCounts..."
export PATH="$feature_counts_dir:$PATH"
cd "$main_dir"
featureCounts -t CDS -g ID -p -O \
	-a "$main_dir/$gene_dir/gene_sequences_all.gff" \
	-o "$main_dir/$gene_dir/gene_counts_all_samples.txt" \
	bam_files/*.bam
echo "*******************************************************"
echo "*******************************************************"	

# this is the end of the metagenomics contig worklow
echo "[$(date -Is)] This concludes the metagenomics contigs workflow"

# align also metatranscriptomics sample reads to the contigs and produce count matrices
echo "Now aligning metatranscriptomics samples mRNA reads to the contigs and providing count information similar to the metagenomics samples"

# run bowtie 2 and align the metatranscriptome samples to the contigs
echo "[$(date -Is)] Running Bowtie2 and aligning the metatranscriptome samples to the contigs..."
mkdir -p "$main_dir/bam_files_metat"
module purge
module load bowtie2/2.4.4
module load samtools/1.18

# map the rRNA filtered metatranscriptomic reads to the coassembled contigs - the interest is mainly in the metabolic marker genes
mapfile -t SAMPLES < "$mt_sample_name_file"

for j in "${!SAMPLES[@]}"; do
   bowtie2 \
        -1 "$mt_fastq_dir/${SAMPLES[$j]}_fwd.fq" \
        -2 "$mt_fastq_dir/${SAMPLES[$j]}_rev.fq" \
        -x "$main_dir/$main_name" \
        -S "$main_dir/${SAMPLES[$j]}.sam" \
        --threads "$SLURM_CPUS_PER_TASK" \
        --no-unal

    samtools view -@ "$SLURM_CPUS_PER_TASK" -F 4 -bS "$main_dir/${SAMPLES[$j]}.sam" |\
        samtools sort -@ "$SLURM_CPUS_PER_TASK" > "$main_dir/bam_files_metat/${SAMPLES[$j]}.bam"
    
    samtools index -@ "$SLURM_CPUS_PER_TASK" "$main_dir/bam_files_metat/${SAMPLES[$j]}.bam"
    rm 	"$main_dir/${SAMPLES[$j]}.sam"
done 
echo "*******************************************************"
echo "*******************************************************"

# finally finally, get transcript counts to prodigal genes for the metatranscriptomic samples with featureCounts
echo "[$(date -Is)] Get counts for all the genes in the contigs for each sample in the metatranscriptomics data with featureCounts..."
module purge
export PATH="$feature_counts_dir:$PATH"
cd "$main_dir"
featureCounts -t CDS -g ID -p -O \
	-a "$main_dir/$gene_dir/gene_sequences_all.gff" \
	-o "$main_dir/$gene_dir/transcript_counts_all_samples.txt" \
	bam_files_metat/*.bam
echo "*******************************************************"
echo "*******************************************************"	
echo "This is it for now. Goodbye."
