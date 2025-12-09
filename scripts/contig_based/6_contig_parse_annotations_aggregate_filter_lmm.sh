#!/bin/bash
#SBATCH --job-name=process_contig_data
#SBATCH --error=6_contig_parse_count_aggregate_filter_lmm_%A_%a_err.txt
#SBATCH --output=6_contig_parse_count_aggregate_filter_lmm_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time 24:00:00
#SBATCH --account project_2005827

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

# the main folder, where all the metatranscriptomics analysis steps and results will be stored
mt_main_dir="$project_root/metatranscriptomics"

# other required directories and information
contig_main_dir="$mg_main_dir/contig_based"
gene_dir="$contig_main_dir/prodigal_genes"
anvio_project_name=Oulanka_contig_based

# preprocessed metadata file for the study
meta_file="$project_root/metadata/Study_Metadata.RData"

# the result directory for the contig-based analysis
contig_res_dir="$project_root/downstream/contig_based"

# the result directory for the LMM analysis
lmm_res_dir="$contig_res_dir/lmm_analysis"

# create the results directories
mkdir -p "$contig_res_dir"
mkdir -p "$lmm_res_dir"

# load used modules 
module load r-env/432 # check version and change if necessary

# parse all the functional annotations, gene and transcript counts and link the annotations to the counts
echo "[$(date -Is)] Parsing all the functional annotations, gene and transcript counts and linking the annotations to the counts..."
Rscript "$r_script_dir/3_CONTIG_parse_annotations_get_gene_transcript_counts.R" "$gene_dir" "$contig_res_dir" "$anvio_project_name"
echo "*******************************************************"
echo "*******************************************************"

# aggregate / compose the functional metabolic matrices
echo "[$(date -Is)] Aggregating gene and transcript counts and composing metabolic functional matrices for downstream..."
Rscript "$r_script_dir/4_CONTIG_compose_func_matrices.R" "$mg_main_dir/kegg_diamond/Matrices_For_Downstream.RData" "$mt_main_dir/kegg_diamond/Matrices_For_Downstream.RData" "$contig_res_dir"
echo "*******************************************************"
echo "*******************************************************"

# filter all the composed matrices
echo "[$(date -Is)] Filter all the datasets for lowly expressed features, parse sample names and order according to metadata..."
Rscript "$r_script_dir/5_CONTIG_filter_metmark_datas.R" "$contig_res_dir/metagenomics" "$contig_res_dir/metatranscriptomics" "$meta_file"
echo "*******************************************************"
echo "*******************************************************"

# perform differential abundance and expression analysis using linear mixed effects models for the nitrogen and methane related genes
# use the same script as for the read-based analysis
echo "[$(date -Is)] Perform linear mixed effects model analysis for the marker gene datasets..."
Rscript "$scripts_dir/R_scripts/read_based/downstream/6b_READ_perform_lmm_for_metmarkdb_datas.R" "$contig_res_dir/metagenomics" "$contig_res_dir/metatranscriptomics" "$lmm_res_dir"
echo "*******************************************************"
echo "*******************************************************"


echo "All done for now. Remember to double check everything! Bye bye."



