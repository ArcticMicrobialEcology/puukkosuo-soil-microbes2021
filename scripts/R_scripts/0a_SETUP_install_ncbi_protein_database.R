args <- commandArgs(trailingOnly = TRUE)
# the first parameter needs to be the directory for custom R packages- where the taxonomizr package is installed
# the second parameter needs to be the directory where the NCBI protein database will be downloaded and saved
# the third parameter needs to be the local temp directory defined by variable $LOCAL_SCRATCH in the batch script
# please use > 100 GB of memory and > 100 GB of local disc space for this script

.libPaths(c(as.character(args[1]), .libPaths()))
library("taxonomizr")

# download the NCBI taxonomy 
setwd(as.character(args[2]))
prepareDatabase(types=c('prot'), sqlFile = 'accessionTaxa.sql', tmpDir = as.character(args[3]), vocal = TRUE, extraSqlCommand=c("pragma temp_store = 2;"))

# print out session info
print("SessionInfo:")
sessionInfo()
