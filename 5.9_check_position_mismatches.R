# checking that the mapped positions are the same for both snp chips
library(tidyverse)

hd_snps <- read_delim("data/sheep_genome/aligned/final_probe_alignments/snp_set_hd.txt", delim = " ") %>% 
    mutate(chr = str_sub(chr, 7, 8)) %>% 
    mutate(chr = as.numeric(chr) - 71) %>% 
    # one snp is a duplicate but with difference alleles, remove it
    filter(snp != "oar3_OAR22_12486013")

ld_snps <- read_delim("data/sheep_genome/aligned/final_probe_alignments/snp_set_ld.txt", delim = " ")%>% 
    mutate(chr = str_sub(chr, 7, 8)) %>% 
    mutate(chr = as.numeric(chr) - 71)

# check that snps with the same names have the same position
overlap <- ld_snps$snp[ld_snps$snp %in% hd_snps$snp]

# which SNPs have different positions in the two sets?
snp_diff <- hd_snps %>% 
    filter(snp %in% overlap) %>% 
    left_join(ld_snps, by = "snp") %>% 
    mutate(pos_equal = ifelse((pos.x - pos.y == 0), 1, 0)) %>% 
    mutate(pos_diff = pos.x - pos.y) %>% 
    filter(pos_equal == 0) %>% 
    select(snp, flag.x, flag.y, chr.x, chr.y, pos.x, pos.y, mapq.x, mapq.y, pos_diff) %>% 
    mutate(keep = ifelse(abs(pos_diff) < 10, 1, 0)) %>% 
    mutate(new_pos_ld = pos.y + pos_diff) 

snp_new_pos <- snp_diff %>% filter(keep == 1) %>% select(snp, pos.x)

discard_snps <- snp_diff %>% filter(keep == 0) %>% .$snp
change_snps <- snp_diff %>% filter(keep == 1) %>% .$snp

# create new snp lists

ld_snps_new <- ld_snps  %>% 
    filter(!(snp %in% discard_snps)) %>% 
    # change positions to match with hd 
    left_join(snp_new_pos, by = "snp") %>% 
    mutate(pos = ifelse(is.na(pos.x), pos, pos.x)) %>% 
    dplyr::select(-pos.x)
 
hd_snps_new <- hd_snps %>% 
    filter(!(snp %in% discard_snps))

# double check that snps with same name have the same position now!
# check that snps with the same names have the same position
overlap <- ld_snps_new$snp[ld_snps_new$snp %in% hd_snps_new$snp]

# which SNPs have different positions in the two sets?
snp_diff2 <- hd_snps_new %>% 
    filter(snp %in% overlap) %>% 
    left_join(ld_snps_new, by = "snp") %>% 
    mutate(pos_equal = ifelse((pos.x - pos.y == 0), 1, 0)) %>% 
    mutate(pos_diff = pos.x - pos.y) %>% 
    # filter unequal positions
    filter(pos_equal == 0) 

ld_snps_new %>% write_delim("data/sheep_genome/aligned/final_probe_alignments/final_snp_sets/snp_set_ld_ram.txt", " ")
hd_snps_new %>% write_delim("data/sheep_genome/aligned/final_probe_alignments/final_snp_sets/snp_set_hd_ram.txt", " ")

# check whether probes are the same
# probes_hd <- read_csv("data/sheep_genome/manifest_files/SheepHD_AgResearch_Cons_15041608_A.csv", 
#                       skip = 7) %>% 
#              rename(snp = Name, seq_hd = AlleleA_ProbeSeq) %>% 
#              select(snp, seq_hd)
# probes_ld <- read_csv("data/sheep_genome/manifest_files/ovinesnp50-genome-assembly-oar-v3-1.csv", 
#                       skip = 7) %>% 
#              rename(snp = Name, seq_ld = AlleleA_ProbeSeq) %>% 
#              select(snp, seq_ld)
# 
# snp_diff_seq <- snp_diff %>% 
#     left_join(probes_hd, by = "snp") %>% 
#     left_join(probes_ld, by = "snp") %>% 
#     mutate(check_seq = ifelse(seq_hd == seq_ld, 1, 0))

# many of them seem to deviate very litte (-+ 3 BP)
# some are very different though, do they have the same probes?
