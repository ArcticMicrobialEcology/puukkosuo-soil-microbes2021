args <- commandArgs(trailingOnly = TRUE)
# the first parameter needs to be the directory for the custom R package installations where the taxonomizr package is installed

# the second parameter needs to be downloaded metabolic marker gene database fasta

# the third parameter needs to be the ncbi database file, where the ncbi protein to taxon mapping database has been downloaded and installed

# the fourth parameter needs to be the file for some manually observed additional eukaryotic proteins in the database

# the fifth parameter needs to be the output directory where the filtered fasta file and the relevant taxonomy are saved

# e.g.
# args <- character(5)
# args[1] <- "/projappl/project_2009164/project_rpackages_430"
# args[2] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/databases/metmarkdb/Funcgenes_51_Dec2021.fasta"
# args[3] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/databases/ncbi_protein/accessionTaxa.sql"
# args[4] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/databases/metmarkdb/Some_manually_observed_eukaryotic_proteins.txt"
# args[5] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/databases/metmarkdb"

# load libraries
.libPaths(c(as.character(args[1]), .libPaths()))
library(taxonomizr)

# read-in the metabolic marker gene fasta file
fasta_file <- readLines(as.character(args[2]))

# pick protein headers
header_locs <- grep(">", fasta_file)
header_texts <- fasta_file[header_locs]

# parse protein names
prot_names <- gsub(">", "", header_texts)
prot_names <- strsplit(x = prot_names, split = "[-]")
prot_names <- unlist(lapply(prot_names, function(x) x[2]))
prot_names <- strsplit(x = prot_names, split = " ")
prot_names <- unlist(lapply(prot_names, function(x) x[1]))

# asssign NCBI taxonomy identifier to the proteins - this will take a while
taxaId <- accessionToTaxa(c(prot_names), sqlFile = as.character(args[3]))

# get taxonomy for the retrieved IDs
prot_taxa <- getTaxonomy(taxaId, sqlFile = as.character(args[3]))
prot_taxa <- data.frame(prot_taxa, stringsAsFactors = F, check.names = F)
prot_taxa$taxid <- taxaId
prot_taxa$protein <- prot_names 
prot_taxa$header <- gsub(">", "", header_texts)

# set the superkingdom as Eukaryota for the manually detected Eukaryotic proteins
man_eukaryotes <- read.csv(file = as.character(args[4]), header = F, sep = "\t", quote = "")
man_eukaryotes <- as.character(man_eukaryotes[,1])
locs_euk <- which(prot_taxa$header%in%man_eukaryotes)
if(length(locs_euk)>0){
  prot_taxa$superkingdom[locs_euk] <- "Eukaryota"
}

# remove Eukaryotic taxa
all_euk_locs <- which(prot_taxa$superkingdom%in%"Eukaryota")

# remove these from the fasta
all_rem_rows <- numeric(0)
for(i in 1:length(all_euk_locs)){
  rem_rows <- c(header_locs[all_euk_locs[i]]:(header_locs[all_euk_locs[i]+1]-1)) 
  all_rem_rows <- c(all_rem_rows, rem_rows)
}

# separate the fasta database into non-eukaryotic and eukaryotic proteins
fasta_non_euk <- fasta_file[-all_rem_rows]
fasta_euk <- fasta_file[all_rem_rows]

# save the fasta
setwd(as.character(args[5]))
writeLines(fasta_non_euk, "MetMarkDB_Euk_Filtered.fasta")
writeLines(fasta_euk, "MetMarkDB_Eukaryotes.fasta")

# save the eukaryotic and non eukaryotic proteins wit their taxa
prot_taxa_non_euk <- prot_taxa[-all_euk_locs,]
prot_taxa_euk <- prot_taxa[all_euk_locs,]

write.table(x = prot_taxa_non_euk, file = "MetMarkDB_Prokaryotic_Proteins.txt", quote = F, row.names = T, col.names = T, sep = "\t")
write.table(x = prot_taxa_euk, file = "MetMarkDB_Eukaryotic_Proteins.txt", quote = F, row.names = T, col.names = T, sep = "\t")

save(prot_taxa_non_euk, file="MetMarkDB_Prokaryotic_Proteins.RData")
save(prot_taxa_euk, file="MetMarkDB_Eukaryotic_Proteins.RData")

# print out session info
print("SessionInfo:")
sessionInfo()

