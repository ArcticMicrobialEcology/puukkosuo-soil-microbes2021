#!/bin/bash
#SBATCH --job-name=setup_databases
#SBATCH --error=0a_setup_databases_%A_%a_err.txt
#SBATCH --output=0a_setup_databases_%A_%a_out.txt
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=256G
#SBATCH --time 71:00:00
#SBATCH --account project_2009164
#SBATCH --gres=nvme:200

# a script to install some external databases required to run the metagenomics / metatranscriptomics workflow on the CSC Puhti supercomputer. 
# these scripts will use several software installed in the Puhti supercomputer. See a separate list of the installed software and version numbers.
# Also some metadata and database files are required to exist in the project directories. See the separate list for this. 
echo "Installing some external databases required to run the metagenomics / metatranscriptomics workflow..."
echo "This will most likely take some time to get everything set up..."
echo "[INFO] Starting database setup at $(date -Is) on $(hostname)"

# define the needed directories
# this is the project root directory
project_root="/scratch/project_2009164/2_OULANKA/Tommi/final"

# this is the directory where all the external needed software not found already in Puhti have been installed to
software_base_dir="/projappl/project_2009164" 

# this is where the required databases will be installed
database_base_dir="$project_root/databases"

# this is the directory where all the required external R packages will be installed
r_package_dir="$software_base_dir/project_rpackages_430"

# the directories for the scripts and the R scripts
scripts_dir="$project_root/scripts"
r_script_dir="$scripts_dir/R_scripts"

# this the directory where the NCBI protein to taxonomy mapping sql database will be downloaded by the taxonomizr R package
ncbi_database_dir="$database_base_dir/ncbi_protein"

# this is the file, where the KEGG DIAMOND prokaryote database is located and copied from
# this database is received from Jenni Hultman
# this database is the KEGG prokaryotic gene amino acid sequences translated into an DIAMOND database.
kegg_prok_db="/scratch/project_2009164/KEGG_DB/PROKARYOTES.dmnd"

# this is the directory where the anvio KEGG databases are downloaded to
kegg_data_dir="$database_base_dir/KEGG"

# set up Anvi'o databases, using a developer version of anvi'o v8-dev, pulled from anvi'o git and installed 3.3.2025
echo "Setting up databases with Anvio..."

export PATH="$software_base_dir/anvio_dev/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/bin:$PATH"
export PATH="$software_base_dir/anvio_git/anvio/sandbox:$PATH"
export PYTHONPATH="${PYTHONPATH:-}:$software_base_dir/anvio_git/anvio"

# some extra modules needed for this installation
module load gcc/11.3.0
module load biopythontools/11.3.0_3.10.6
module load hmmer/3.4

# download and setup KEGG modules and related information using anvio
# echo "Using the latest snapshot in Anvi'o dev, should be new enough. Check date from output."
echo "Downloading and parsing KEGG metabolic module data using Anvio..."
mkdir -p "$database_base_dir" "$kegg_data_dir"
cd "$database_base_dir"
anvi-setup-kegg-data -T "$SLURM_CPUS_PER_TASK" --kegg-data-dir KEGG
echo "*******************************************************"
echo "*******************************************************"
module purge

# use the taxonomizr R package to install the NCBI protein to taxon ID mapping database - currently installed version date 22_5_2024
module load r-env/432
echo "Using the taxonomizr R package to get the NCBI taxon id mappings..."
mkdir -p "$ncbi_database_dir"
Rscript "$r_script_dir/0a_SETUP_install_ncbi_protein_database.R" "$r_package_dir" "$ncbi_database_dir" "$LOCAL_SCRATCH" 

# copying the KEGG prokaryotic database into the database directory
echo "Copying the DIAMOND KEGG Prokaryotic amino acid sequence database..."
cd "$database_base_dir"
mkdir -p KEGG_PROK
cp "$kegg_prok_db" "$database_base_dir/KEGG_PROK"
echo "*******************************************************"
echo "*******************************************************"
echo "All done. Bye"
