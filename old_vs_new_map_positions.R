library(snpStats)
library(tidyverse)

# plink name
sheep_plink_name <- "data/SNP_chip/Plates_1to87_QC3"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)


# old positions
map_old <- as_tibble(full_sample$map)

# new positions
map_new <- read_csv("data/SNP_chip/new_mapping_positions/OvineSNP50_B.csv.rambo_pos_20190513.csv")

not_in_new_map <- map_old$snp.name[!(map_old$snp.name %in% map_new$entry)]

test <- map_old %>% 
    filter(snp.name %in% not_in_new_map)

str_match("lol", "lo")
library(stringdist)
which(str_detect(map_new$entry, "OAR2_126564819.1"))
out <- stringdist(map_new$entry, "s18609.1")
min(out)
hist(out)
map_new$entry[out == 1]

# 
# plink name
sheep_plink_name <- "data/SNP_chip/Plates_1-2_HD_QC2"
sheep_plink_name <- "data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram"

# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample2 <- read.plink(sheep_bed, sheep_bim, sheep_fam)

map_new_hd <- read_csv("data/SNP_chip/new_mapping_positions/SheepHD_AgResearch_Cons_15041608_A.csv.rambo_pos_20190513.csv")
map_old_hd <- full_sample2$map

not_in_new_map_hd <- map_old_hd$snp.name[!(map_old_hd$snp.name %in% map_new_hd$entry)]
