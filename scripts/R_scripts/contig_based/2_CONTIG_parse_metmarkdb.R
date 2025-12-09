args <- commandArgs(trailingOnly = TRUE)

# the first input parameter needs to be the gene directory for the progigal genes
# the second input parameter needs to be prepared functional marker gene information RData file
# the third input parameter needs to be if the MetMarkDB hits should be filtered (1) or not (0) with the defined default limits.
# currently the default limits: e-value 10^-6 and bitscore 50, change if needed

# e.g.
# args <- character(3)
# args[1] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/metagenomics/contig_based/prodigal_genes"
# args[2] <- "/scratch/project_2009164/2_OULANKA/Tommi/final/databases/metmarkdb/Gene_metadata.RData"
# args[3] <- "0"

# a small script to parse the results the DIAMOND mapping of the Anvio called prodigal genes to the metabolic marker gene database into
# a suitable format for Anvio

# read in  the DIAMOND blastx aligment results - this format is fixed and follows the one defined in the shell script for the DIAMOND search run.
setwd(args[1])
genes_func_align <-read.csv("prodigal_gene_metmarkdb_hits.txt", sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F, 
                            colClasses = c("character","character","numeric", "integer", "integer", "integer", 
                                           "integer", "integer","numeric","numeric", "numeric", "numeric", "integer"))
colnames(genes_func_align) <- c("qseqid", "stitle", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovhsp", "scovhsp", "slen")


# what was imported?
print(paste("Read-in results file with", nrow(genes_func_align), "hits"))

# an extra option to filter the alignment hits with more stringent criteria
if(as.numeric(args[3]) == 1){
  print("Filtering hits with evalue <= 10^-6 and bitscore >= 50")
  genes_func_align_filtered <- genes_func_align[which(genes_func_align$evalue<= 10^(-6) & genes_func_align$bitscore>50),]
  
  print(paste("hits remaining after filtering:", nrow(genes_func_align_filtered)))
}else{
  genes_func_align_filtered <- genes_func_align
}

# remove possible ribosomal RNA hits, included in the anvi-get-sequences-for-gene-calls - shouldn't be too many if any.
# but no in the exported amino acid sequences
prodigal_genes <- read.csv("prodigal_genes_no_aaseq.txt", sep = "\t", header = F, quote = "", fill = F, stringsAsFactors = F)
colnames(prodigal_genes) <- prodigal_genes[1,]
prodigal_genes <- prodigal_genes[-1,]

keep <- which(genes_func_align_filtered$qseqid%in%prodigal_genes$gene_callers_id)
genes_func_align_filtered <- genes_func_align_filtered[keep,]

# parse into the correct format anvio
# format should be:
#  gene_callers_id	source	accession	function	e_value

# examples of genes_func_align_filtered$stitle  
# NuoF-WP_129628804.1 NADH-quinone oxidoreductase subunit NuoF [Candidatus Oscillochloris fontis]
# FdhA-WP_036701914.1 formate dehydrogenase subunit alpha [Paracoccus sanguinis]
# FdhA-WP_085053520.1 formate dehydrogenase subunit alpha [Nitrospirae bacterium HCH-1]
# FdhA-WP_085053520.1 formate dehydrogenase subunit alpha [Nitrospirae bacterium HCH-1]
# NuoF-WP_120169180.1 4Fe-4S binding protein, partial [Thermohalobacter berrensis]
# NrfA-WP_015351334.1 - Myxococcus stipitatus
# SdhA_FrdA-Candidatus Protochlamydia amoebophila (Group 1b)

# these need to be processed separately from the rest of the marker genes as they are annotated differently in the DB
locs_sdha <- grep("\\bSdhA_FrdA\\b", genes_func_align_filtered$stitle) 

# separate based on space " "
func_genes_desc <- lapply(genes_func_align_filtered$stitle, function(x) strsplit(x = as.character(x[1]), split = " "))

# separate the first entry from func_genes_desc based on "-" get metabolic marker gene names for the hits
func_genes <- unlist(lapply(func_genes_desc, function(x) strsplit(x[[1]], "-")[[1]][1]))

# separate the second entry from func_genes_desc based on "-" get protein annotations for the hits
func_genes_proteins <- unlist(lapply(func_genes_desc, function(x) strsplit(x[[1]], "-")[[1]][2]))

# process SdhA_FrdA separately, different type of naming convention in the database - no protein annotation
func_genes_proteins[locs_sdha] <- "unknown"

# parse the protein in also the protein descriptions
func_genes_desc2 <- unlist(lapply(func_genes_desc, function(x) {
  if(x[[1]][2]=="-"){
    paste(x[[1]][3:length(x[[1]])], collapse = " ")
  }else{
    paste(x[[1]][2:length(x[[1]])], collapse = " ")
  }
}))

# parse SdhA_FrdA still separately
func_genes_desc3 <- unlist(lapply(genes_func_align_filtered$stitle, function(x) strsplit(x, split = "-")[[1]][2]))

# add SdhA_FrdA to the same description with others
func_genes_desc2[locs_sdha] <- func_genes_desc3[locs_sdha]

# make a data frame for Anvio
anvio_functions <- data.frame(matrix(nrow = nrow(genes_func_align_filtered), ncol=5))
colnames(anvio_functions) <- c("gene_callers_id",	"source",	"accession", "function",	"e_value")
anvio_functions[,1] <- genes_func_align_filtered$qseqid
anvio_functions[,2] <- rep("MetMarkDB", nrow(genes_func_align_filtered))
anvio_functions[,3] <- paste(func_genes, func_genes_proteins, sep = "_")
anvio_functions[,4] <- func_genes_desc2  
anvio_functions[,5] <- genes_func_align_filtered$evalue

# read-in marker gene metadata, needs to be precompiled already
load(args[2])
gene_metadata$Definition <- paste(gene_metadata$Definition1, gene_metadata$Definition2, sep = ";")

# parse in the marker gene names to the protein annot
temp_genes <- unlist(lapply(anvio_functions$accession, function(x) strsplit(x = x, split = "_")[[1]][1]))
locs_sdha <- grep("\\bSdhA\\b", temp_genes) # correct these
temp_genes[locs_sdha] <-"SdhA_FrdA"
all(temp_genes%in%rownames(gene_metadata)) #TRUE

# add still gene name to the protein description - easier parsing downstream
anvio_functions$`function` <- paste(temp_genes, anvio_functions$`function`, sep = "; ")

# some example rows of what the parsed input data.frame for anvio will look like
# gene_callers_id	source	accession	function	e_value
# 49	MetMarkDB	NuoF_WP_129628804.1	NuoF; NADH-quinone oxidoreductase subunit NuoF [Candidatus Oscillochloris fontis]	2.72e-114
# 182	MetMarkDB	FdhA_WP_036701914.1	FdhA; formate dehydrogenase subunit alpha [Paracoccus sanguinis]	6.21e-33
# 329	MetMarkDB	FdhA_WP_085053520.1	FdhA; formate dehydrogenase subunit alpha [Nitrospirae bacterium HCH-1]	1.88e-28
# 330	MetMarkDB	FdhA_WP_085053520.1	FdhA; formate dehydrogenase subunit alpha [Nitrospirae bacterium HCH-1]	3.28e-158
# 331	MetMarkDB	NuoF_WP_120169180.1	NuoF; 4Fe-4S binding protein, partial [Thermohalobacter berrensis]	6.15e-172
# 424	MetMarkDB	NrfA_WP_015351334.1	NrfA; Myxococcus stipitatus	2.87e-233
# 440	MetMarkDB	SdhA_FrdA_unknown	SdhA_FrdA; Candidatus Protochlamydia amoebophila (Group 1b)	8.68e-09

# save the diamond hits formatted to be inputted into anvio
write.table(x = anvio_functions, file = "MetMarkDB_for_AnvioFunc.txt", quote = F, col.names = T, row.names = F, sep = "\t")

# print out session info
print("SessionInfo:")
sessionInfo()