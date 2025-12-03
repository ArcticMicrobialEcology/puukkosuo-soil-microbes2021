args <- commandArgs(trailingOnly = TRUE)
# requires the output directory for the downloaded KEGG mappings as the first parameter

# a simple script to retrieve KEGG gene to KO group mappings through the KEGG REST API

# load libraries
library(KEGGREST)

# get all databases
all_databases <- listDatabases()

# get all organisms
orgs <- keggList("organism")

# only include prokaryotes
proks <- orgs[grep("Prokaryotes", orgs[,4]),]

# get all KO identifiers
kos <- keggList("ko")

# make a timestamp of date
timestamp <- gsub("-", "_", as.character(Sys.Date()))

# save the downloaded mappings into this directory
setwd(as.character(args[1]))

# define output file - now hardcoded, could also be as an input parameter?
outfile <- "KO_Mapping_Genes_Prokaryotes.txt"

# add timestamp
outfile <- paste(timestamp, outfile, sep = "_")

# link KO identifiers to genes in groups of 100
n_loops <- ceiling((length(kos)/100))

# initialize index
ind <- 1

# open output file connection
outConnection <- file(outfile, "a")

print("Starting the download of KO mappings to genes")
for(i in 1:n_loops){
  # get the linking to genes
  if(i<n_loops){ # not the last loop
    query_res <- tryCatch({
      keggLink("genes", c(names(kos)[ind:(ind+99)]))
    },error = function(e) {
      NA
    })
    expected_length <- 100
  }else{
    query_res <- tryCatch({
      keggLink("genes", c(names(kos)[ind:length(kos)])) #Last part
    },error = function(e) {
      NA
    })
    expected_length <- length(ind:length(kos))
  }
  
  if(length(query_res)==1){
    if(is.na(query_res)){
      print(paste("There was a problem retrieving anything at iteration:", i))
    }
  }else if(length(unique(names(query_res)))!=expected_length){
    print(paste("Gene annotations for all the KOs were not retrieved at iteration:", i))
  }
  
  # pick only the species name
  query_parsed <- unlist(lapply(strsplit(query_res, ":"), function(x) x[[1]][1]))
  
  # remove non-prokaryotes
  rems <- which(!query_parsed%in%proks[,2])
  
  if(length(rems)>0){
    query_filt <- query_res[-rems]
  }else{
    query_filt <- query_res
  }
  
  # write to file
  if(length(query_filt)>0){
    textToWrite <- paste(query_filt, gsub("ko:", "",names(query_filt)), sep = "\t")
    writeLines(c(textToWrite), outConnection)
  }
  
  # print out progress
  print(paste((100 * round(i/n_loops,3)), "% of gene mappings downloaded", sep = ""))
  
  # increase the index
  ind <- ind+100
  
  # sleep a small while as not to flood the server
  Sys.sleep(0.2)
}
# close output file connection
close(outConnection)

print("Done downloading the KO mappings to genes")

# parse into matrix
ko_matrix <- matrix(nrow = length(kos), ncol=2)
ko_matrix[,1] <- names(kos)
ko_matrix[,2] <- kos

# save the downloaded mappings - filename hard coded for now
write.table(x = ko_matrix, file = paste(timestamp, "_KO_Name_Mapping_REST.txt", sep = ""), row.names = F, col.names = F, sep = "\t", quote = F)

# print out session info
print("SessionInfo:")
sessionInfo()

# save sessionInfo
s_info <- sessionInfo()
save(s_info, file = gsub(".txt", ".RData", outfile))
