# check oar31 50k and HD chips before merging.
# We want to check:
# (1) SNP at same positions have the same names
# (2) SNP with the same names have the same positions

library(tidyverse)
library(snpStats)

# plink name
sheep_plink_name <- "data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
hd_chip <- read.plink(sheep_bed, sheep_bim, sheep_fam)
map_hd <- hd_chip$map

# plink name
sheep_plink_name <- "data/SNP_chip/oar31_mapping/Plates_1to87_QC3"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
ld_chip <- read.plink(sheep_bed, sheep_bim, sheep_fam)
map_ld <- ld_chip$map
map_ld %>% filter(chromosome %in% 1:26) %>% as_tibble()
# overlap between chips
sum(map_hd$snp.name %in% map_ld$snp.name)

# check snps with different names but same position
same_pos <- map_ld %>% 
                filter(chromosome %in% 1:27) %>% 
                inner_join(map_hd, by = c("chromosome", "position")) %>% 
                mutate(same_name = ifelse(snp.name.x == snp.name.y, TRUE, FALSE))

# get names from ld chip which differ on hd chip
same_pos %>% 
    filter(!same_name) %>% 
    .$snp.name.x -> diff_names_ld

same_pos %>% 
    filter(!same_name) %>% 
    .$snp.name.y -> diff_names_hd

# check snps with same names but different position
same_name <- map_ld %>% 
    filter(!is.na(chromosome)) %>% 
    filter(chromosome %in% 1:27) %>% 
    inner_join(map_hd, by = c("snp.name")) %>% 
    mutate(same_pos = ifelse(position.x == position.y, TRUE, FALSE))
same_name %>% filter(!same_pos)


# filter SNPs on LD chip and save as file
map_ld %>% 
    # filter everything unmapped
    filter(chromosome %in% 1:27) %>% 
    as_tibble() %>% 
    filter(!(snp.name %in% diff_names_ld)) %>% 
    .$snp.name %>% 
    write_lines("data/oar31_snps_inc_for_merge_ld.txt")

map_hd %>%  # filter everything unmapped
    filter(chromosome %in% 1:27) %>% 
    as_tibble() %>% 
    filter(!(snp.name %in% diff_names_hd)) %>% 
    .$snp.name %>% 
    write_lines("data/oar31_snps_inc_for_merge_hd.txt")


# make snp filtered bed files
# HD
system(paste0("/usr/local/bin/plink --bfile data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2 --sheep ",      
              "--extract data/oar31_snps_inc_for_merge_hd.txt ",
              "--make-bed --out data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2_snpfile_for_merge"))
# LD          
system(paste0("/usr/local/bin/plink --bfile data/SNP_chip/oar31_mapping/Plates_1to87_QC3 --sheep ",      
              "--extract data/oar31_snps_inc_for_merge_ld.txt ",
              "--make-bed --out data/SNP_chip/oar31_mapping/Plates_1to87_QC3_snpfile_for_merge"))

# merge chips based on oar 3.1 mapping
system(paste0("/usr/local/bin/plink --bfile data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2_snpfile_for_merge --sheep ",      # --out output/ROH/roh_nofilt
              "--bmerge data/SNP_chip/oar31_mapping/Plates_1to87_QC3_snpfile_for_merge.bed data/SNP_chip/oar31_mapping/Plates_1to87_QC3_snpfile_for_merge.bim data/SNP_chip/oar31_mapping/Plates_1to87_QC3_snpfile_for_merge.fam ",
              "--make-bed --out data/SNP_chip/oar31_mapping/merged_sheep_geno_oar31 ",
              "--indiv-sort 0 --merge-equal-pos "))

# merged data
sheep_plink_name <- "data/SNP_chip/oar31_mapping/merged_sheep_geno_oar31"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
merged_genos <- read.plink(sheep_bed, sheep_bim, sheep_fam)
# all
dim(merged_genos$genotypes)
# without X
merged_genos$map %>% filter(chromosome %in% 1:26) %>% as_tibble()
