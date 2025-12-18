# a small script to prepare some environmental variables for use for the later plotting and analysis scripts.

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(readxl)

# read-in the gas flux data
setwd(paste(project_root, "/metadata/environmental_data", sep = ""))

gas_fluxes1 <- data.frame(read_excel("2021-09-12_Oulanka_flux.xlsx", sheet = 4), stringsAsFactors = F, check.names = F)
gas_fluxes2 <- data.frame(read_excel("2021-09-20_Oulanka_flux.xlsx", sheet = 4), stringsAsFactors = F, check.names = F)

# remove empty rows from some tables
gas_fluxes1 <- gas_fluxes1[-which(gas_fluxes1$`CHAMBER ID [#]`==0),]

# load metadata
load(paste(project_root, "/metadata/Study_Metadata.RData", sep = ""))

# organize all methane flux data similarly - according to sample plot number
rownames(gas_fluxes1) <- gas_fluxes1$`CHAMBER ID [#]`
rownames(gas_fluxes2) <- gas_fluxes2$`CHAMBER ID [#]`

# name sample metadata also according to the plot number
rownames(metadata) <- metadata$Plot

# organize the methjane flux data according to the metadata
gas_fluxes1 <- gas_fluxes1[rownames(metadata),]
gas_fluxes2 <- gas_fluxes2[rownames(metadata),]

# combine the methane fluxes in the same table
methane_fluxes <- data.frame(gas_fluxes1$`CH4 FLUX [mg m-2 h-1]`,
                             gas_fluxes2$`CH4 FLUX [mg m-2 h-1]`)

# name the fluxes according to the collection date and plot number
colnames(methane_fluxes) <- c("12_9_2021", "20_9_2021")
rownames(methane_fluxes) <- rownames(metadata)

# also read-in pore water data from Summer 2021
pore_water <- data.frame(read_excel("Huokosvesi_Puukkosuo_kesä2021.xlsx"), stringsAsFactors = F, check.names = F)

# for pore water, the closest timepoint for which data is available to the microbial sampling time 
# is september, syyskuu. Use this data
pore_water_september <- pore_water[which(pore_water$Kuukausi=="syyskuu"),]

# organize the pore water data similarly to the gas flux data, 
# according to the sample plot numberin in sample metadata
rownames(pore_water_september) <- pore_water_september$RuutuID
pore_water_september <- pore_water_september[rownames(metadata),]

# set the relevant variables as numeric
colnames(pore_water_september)
for(i in 5:11){pore_water_september[,i] <- as.numeric(pore_water_september[,i])}
pore_water_nitrogen <- pore_water_september[,c(6,10)]

# save the parsed tables for later exploration & analysis with relevant microbial data
ses_info <- sessionInfo()
save(methane_fluxes, pore_water, pore_water_september, pore_water_nitrogen, metadata, ses_info, file = "Meth_fluxes_porewater.RData")

# print out session info
print("SessionInfo:")
sessionInfo()
