# generate a small mapping from Silva phylum names to GTDB phylum names for mapping

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(tibble)
silva_gtdb_map <- tribble(
  ~SILVA,                   ~GTDB,
  "Crenarchaeota",          "Thermoproteota",
  "Halobacterota",          "Halobacteriota",
  "Nanoarchaeota",          "Nanoarchaeota",
  "Aenigmarchaeota",        "Aenigmarchaeota",
  "Altiarchaeota",          "Altiarchaeota",
  "Asgardarchaeota",        "Asgardarchaeota",
  "Euryarchaeota",          "Halobacteriota/Methanobacteriota/Thermoplasmatota",
  "Hadarchaeota",           "Hadarchaeota",
  "Hydrothermarchaeota",    "Hydrothermarchaeota",
  "Iainarchaeota",          "Iainarchaeota",
  "Micrarchaeota",          "Micrarchaeota",
  "Thermoplasmatota",       "Thermoplasmatota",
  "Acidobacteriota",        "Acidobacteriota",
  "Actinobacteriota",       "Actinomycetota",
  "Bacteroidota",           "Bacteroidota",
  "Chloroflexi",            "Chloroflexota",
  "Cyanobacteria",          "Cyanobacteriota",
  "Desulfobacterota",       "Desulfobacterota (A-E)",
  "Entotheonellaeota",      "Tectomicrobia",
  "Firmicutes",             "Firmicutes (A-Z)",
  "Gemmatimonadota",        "Gemmatimonadota",
  "Latescibacterota",       "Latescibacterota",
  "MBNT15",                 "MBNT15",
  "Methylomirabilota",      "Methylomirabilota",
  "Myxococcota",            "Myxococcota",
  "NB1-j",                  "NB1-j",
  "Nitrospirota",           "Nitrospirota",
  "Patescibacteria",        "Patescibacteria",
  "Planctomycetota",        "Planctomycetota",
  "Proteobacteria",         "Pseudomonadota",
  "Spirochaetota",          "Spirochaetota",
  "Sva0485",                "Sva0485",
  "Verrucomicrobiota",      "Verrucomicrobiota",
  "Zixibacteria",           "Zixibacteria",
  "10bav-F6",               "10bav-F6",
  "Abditibacteriota",       "Abditibacteriota",
  "Acetothermia",           "Acetothermiota",
  "Aerophobota",            "Aerophobota",
  "AncK6",                  "AncK6",
  "Aquificota",             "Aquificota",
  "Armatimonadota",         "Armatimonadota",
  "Bdellovibrionota",       "Bdellovibrionota",
  "BHI80-139",              "BHI80-139",
  "Caldatribacteriota",     "Caldatribacteriota",
  "Caldisericota",          "Caldisericota",
  "Calditrichota",          "Calditrichota",
  "Campylobacterota",       "Campylobacterota",
  "Cloacimonadota",         "Cloacimonadota",
  "Dadabacteria",           "Dadabacteria",
  "Deferribacterota",       "Deferribacterota",
  "Deferrisomatota",        "Deferrisomatota",
  "Deinococcota",           "Deinococcota",
  "Dependentiae",           "Dependentiae",
  "Desantisbacteria",       "Desantisbacteria",
  "Dictyoglomota",          "Dictyoglomota",
  "DTB120",                 "DTB120",
  "Edwardsbacteria",        "Edwardsbacteria",
  "Elusimicrobiota",        "Elusimicrobiota",
  "FCPU426",                "FCPU426",
  "Fermentibacterota",      "Fermentibacterota",
  "Fibrobacterota",         "Fibrobacterota",
  "Firestonebacteria",      "Firestonebacteria",
  "Fusobacteriota",         "Fusobacteriota",
  "FW113",                  "FW113",
  "GAL15",                  "GAL15",
  "GN01",                   "GN01",
  "Halanaerobiaeota",       "Halanaerobiaeota",
  "Hydrogenedentes",        "Hydrogenedentes",
  "LCP-89",                 "LCP-89",
  "Margulisbacteria",       "Margulisbacteria",
  "Marinimicrobia (SAR406 clade)", "Marinimicrobia",
  "MAT-CR-M4-B07",          "MAT-CR-M4-B07",
  "Modulibacteria",         "Modulibacteria",
  "Nitrospinota",           "Nitrospinota",
  "NKB15",                  "NKB15",
  "PAUC34f",                "PAUC34f",
  "Poribacteria",           "Poribacteria",
  "RCP2-54",                "RCP2-54",
  "Rs-K70 termite group",   "Rs-K70",
  "SAR324 clade(Marine group B)", "SAR324",
  "Schekmanbacteria",       "Schekmanbacteria",
  "Sumerlaeota",            "Sumerlaeota",
  "Synergistota",           "Synergistota",
  "TA06",                   "TA06",
  "Thermotogota",           "Thermotogota",
  "TX1A-33",                "TX1A-33",
  "WOR-1",                  "WOR-1",
  "WPS-2",                  "Candidatus Eremiobacterota",
  "WS1",                    "WS1",
  "WS2",                    "WS2",
  "WS4",                    "WS4",
  "unknown",                "unclassified"
)

# transform into data.frame for easier use later
silva_gtdb_map <- data.frame(silva_gtdb_map, stringsAsFactors = F, check.names = F)


# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure3", sep = ""))
save(silva_gtdb_map, file = "Silva_GTDB_phylum_mapping.RData")

# print out session info
print("SessionInfo:")
sessionInfo()