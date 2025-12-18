# a small script to generate some more colors to be used for the MAG plotting

# define some needed directories
project_root="" # e.g. "/scratch/project_number/puukkosuo"

# load libraries
library(colorspace)

# convert any mix of color names/hex to hex strings
to_hex <- function(x) {
  is_hex <- grepl("^#", x)
  out <- x
  if (any(!is_hex)) {
    rgb_mat <- grDevices::col2rgb(x[!is_hex])
    out[!is_hex] <- grDevices::rgb(t(rgb_mat) / 255)
  }
  out
}

# greedy farthest-point palette picker using colorspace only
pick_distinct_palette <- function(existing_cols, k = 20, seed = 1) {
  set.seed(seed)
  
  # 1) Existing colors -> HEX -> LAB
  existing_hex <- to_hex(existing_cols)
  existing_lab <- coords(as(hex2RGB(existing_hex), "LAB"))  # matrix (n x 3)
  
  # 2) Build a large candidate pool of vivid qualitative colors (already HEX)
  cand_hex <- unique(c(
    qualitative_hcl(800, h = c(0, 360), c = 70, l = 60),
    qualitative_hcl(800, h = c(0, 360), c = 80, l = 55),
    qualitative_hcl(800, h = c(0, 360), c = 60, l = 70)
  ))
  
  # Filter to mid L / decent C in HCL (polarLUV): keep vivid, avoid too light/dark
  cand_hcl <- coords(as(hex2RGB(cand_hex), "polarLUV"))  # columns: L, C, H
  keep <- cand_hcl[,"C"] >= 35 & cand_hcl[,"L"] >= 35 & cand_hcl[,"L"] <= 80
  cand_hex <- cand_hex[keep]
  
  # Remove any color already present
  cand_hex <- setdiff(cand_hex, existing_hex)
  
  # Precompute LAB for candidates
  cand_lab <- coords(as(hex2RGB(cand_hex), "LAB"))
  
  # 3) Greedy selection: each step pick the candidate with the largest
  #    minimum LAB distance to the already-selected (incl. existing)
  selected_hex <- character(0)
  selected_lab <- existing_lab
  
  for (i in seq_len(k)) {
    # compute min distance to current selected set
    # (vectorized: for each candidate, compute distances to all selected, take min)
    dmin <- apply(cand_lab, 1, function(row) {
      min(sqrt(rowSums((selected_lab - matrix(row, nrow = nrow(selected_lab), ncol = 3, byrow = TRUE))^2)))
    })
    
    pick <- which.max(dmin)
    
    # append picked color
    selected_hex <- c(selected_hex, cand_hex[pick])
    selected_lab <- rbind(selected_lab, cand_lab[pick, , drop = FALSE])
    
    # remove from candidates
    cand_hex <- cand_hex[-pick]
    cand_lab <- cand_lab[-pick, , drop = FALSE]
  }
  
  selected_hex
}

# colors defined for figure 3
plot_colors <- c(
  "purple","cornflowerblue","mistyrose2","darkgrey","deeppink","brown","violet",
  "firebrick1","chocolate3","forestgreen","darkslateblue","darkolivegreen",
  "chocolate1","cadetblue1","darkgoldenrod3","blue","darkolivegreen3",
  "yellow","burlywood2","black","green"
)

# pick more colors that would be as disctint from the previous colors as possible
new20 <- pick_distinct_palette(plot_colors, k = 20, seed = 1)

# change directory into the plotting directory - needs to exist
setwd(paste(project_root, "/downstream/manuscript_figures/figure6", sep = ""))

# save these colors
save(new20, file="new_colors.RData")
