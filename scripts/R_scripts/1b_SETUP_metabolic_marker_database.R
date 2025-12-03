args <- commandArgs(trailingOnly = TRUE)
# requires as input parameter the directory containing the database and the required information.

# a script to setup the metabolic marker gene DB.
# modify as needed

# The functional 51 marker genes protein database and associated metadata has been downloaded from google drive
# https://drive.google.com/drive/folders/1IeDF4iUbTtPhuaBOyQkrAlXO67EbkJ-A
# Contact person Dr Pok Man Leung (Bob)
# Department of Microbiology, Biomedicine Discovery Institute, Monash University, Clayton, Victoria 3800, Australia

# The database:https://bridges.monash.edu/collections/_/5230745
# Leung, Pok Man; Greening, Chris (2020). Greening lab metabolic marker gene databases. 
# Monash University. Collection. https://doi.org/10.26180/c.5230745

# A marine microbiome study using this database
# https://doi.org/10.1038/s41564-023-01322-0

# the directory containing the database and the required information, such as the following file.
setwd(args[1]) 

# metadata about the database entries
master_lengths <- read.delim("all_lengths_51_dbs_subgroup_Jan2022.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE, strip.white = TRUE)
names(master_lengths) <- c("Sequence", "Length_aa", "Gene", "Length_kb", "Subgroup", "Definition1", "Definition2")

# parse gene names
master_lengths$Gene=gsub("-", "_", master_lengths$Gene)

# get unique definitions for the genes
unique_genes <- unique(master_lengths$Gene)
unique_genes <- unique_genes[-which(unique_genes%in%"Copper_containing_nitrite_reductase_NirK_Jan2021")]
gene_metadata <- matrix(nrow = length(unique_genes), ncol=2)

# parse definitions for the genes
for(i in 1:length(unique_genes)){
  locs <- which(master_lengths$Gene%in%unique_genes[i])[1] # all have the same definitions
  gene_metadata[i,1] <- paste(unique(master_lengths$Definition1[locs]), collapse=";")
  gene_metadata[i,2] <- paste(unique(master_lengths$Definition2[locs]), collapse=";")
}
gene_metadata <- data.frame(gene_metadata, stringsAsFactors = F, check.names = F)
rownames(gene_metadata) <- unique_genes
colnames(gene_metadata) <- c("Definition1", "Definition2")

# save
save(gene_metadata, file = "Gene_metadata.RData")

# print out session info
print("SessionInfo:")
sessionInfo()
