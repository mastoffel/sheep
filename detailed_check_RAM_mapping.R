# throuroughly check whether there are mistakes in the mapping to Ram.
# or whether it's assembly errors
library(tidyverse)
library(snpStats)
source("../sheep_ID/theme_simple.R")
# check HD chip

# old mapping
sheep_plink_name <- "data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2"
# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample_old <- read.plink(sheep_bed, sheep_bim, sheep_fam)
map_old <- full_sample_old$map %>% 
                dplyr::rename(chr_old = chromosome,
                       pos_old = position) %>% 
                as_tibble()

# new mapping
sheep_plink_name <- "data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram"
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
map_new <- full_sample$map %>% 
                dplyr::rename(chr_new = chromosome,
                              pos_new = position) %>% 
                as_tibble()

# new mapping by rudi
map_rudi <- read_csv("data/sheep_genome/aligned_by_rudi/SheepHD_AgResearch_Cons_15041608_A.csv.rambo_pos_20190513.csv") %>% 
    dplyr::rename(snp.name = entry,
                  pos_rudi = pos,
                  chr_rudi = chrom)

wrongly_mapped <- map_new %>% 
    left_join(map_old, by = "snp.name") %>% 
    filter(!is.na(chr_old) & !is.na(chr_new)) %>% 
    filter(chr_new != chr_old) %>% 
    write_delim("../sheep_ID/data/wrongly_mapped_hd.txt")

map_full <- map_new %>% 
            left_join(map_old, by = "snp.name") %>% 
            left_join(map_rudi, by = "snp.name") %>% 
            filter(!is.na(chr_old) & !is.na(chr_rudi) & !is.na(chr_new)) %>% 
            mutate(chr_diff = case_when(
                                        (chr_new != chr_old) & (chr_old == chr_rudi) ~ "only new (n=61)",
                                        (chr_new != chr_old) & (chr_old != chr_rudi) ~ "new and cons (n=3868)",
                                        (chr_new == chr_old) & (chr_old != chr_rudi) ~ "only cons (n=27)",
                                        (chr_new == chr_old) ~ "no (prob. correct)",
                                        )) %>% 
            mutate(chr_diff = factor(chr_diff, levels = c("no (prob. correct)", "new and cons (n=3868)", "only new (n=61)", "only cons (n=27)")))
            


table(map_full$chr_diff)
p <- ggplot(map_full, aes(pos_old, pos_new, color = chr_diff)) +
    geom_point(size = 1, alpha = 1) +
    #theme_simple(axis_lines = TRUE, grid_lines = TRUE, base_size = 10) +
    #scale_color_brewer(type = "qual", palette = "Set2") +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_color_discrete(name = "Mapped to different Chr.?") +
    facet_wrap(~chr_new) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("HD chip positions on oar3.1 (old) plotted against ramb v1.0 (new)",
            subtitle = "And compared to the new mapping position by the consortium (cons)") 

ggsave("figs/hd_new_vs_old_vs_consnew.jpg", p, width = 11, height = 7)


# LD chip

# plink name
sheep_plink_name <- "data/SNP_chip/oar31_mapping/Plates_1to87_QC3"

# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")

full_sample_old <- read.plink(sheep_bed, sheep_bim, sheep_fam)
map_old <- full_sample_old$map %>% 
    dplyr::rename(chr_old = chromosome,
                  pos_old = position) %>% 
    as_tibble()

# plink name
sheep_plink_name <- "data/SNP_chip/ramb_mapping/Plates_1to87_QC4_ram"

# read merged plink data
sheep_bed <- paste0(sheep_plink_name, ".bed")
sheep_bim <- paste0(sheep_plink_name, ".bim")
sheep_fam <- paste0(sheep_plink_name, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
map_new <- full_sample$map %>% 
    dplyr::rename(chr_new = chromosome,
                  pos_new = position) %>% 
    as_tibble()

wrongly_mapped <- map_new %>% 
    left_join(map_old, by = "snp.name") %>% 
    filter(!is.na(chr_old) & !is.na(chr_new)) %>% 
    filter(chr_new != chr_old) %>% 
    write_delim("../sheep_ID/data/wrongly_mapped_ld.txt")

# rudi alignment
map_rudi <- read_csv("data/sheep_genome/aligned_by_rudi/OvineSNP50_B.csv.rambo_pos_20190513.csv") %>% 
    dplyr::rename(snp.name = entry,
                  pos_rudi = pos,
                  chr_rudi = chrom)

map_full <- map_new %>% 
    left_join(map_old, by = "snp.name") %>% 
    left_join(map_rudi, by = "snp.name") %>% 
    filter(!is.na(chr_old) & !is.na(chr_rudi) & !is.na(chr_new)) %>% 
    mutate(chr_diff = case_when(
        (chr_new != chr_old) & (chr_old == chr_rudi) ~ "only new (n=3)",
        (chr_new != chr_old) & (chr_old != chr_rudi) ~ "new and cons (n=326)",
        (chr_new == chr_old) & (chr_old != chr_rudi) ~ "only cons (n=2)",
        (chr_new == chr_old) ~ "no (prob. correct)",
    )) %>% 
    mutate(chr_diff = factor(chr_diff, levels = c("no (prob. correct)", "new and cons (n=326)", "only new (n=3)", "only cons (n=2)")))



table(map_full$chr_diff)
p <- ggplot(map_full, aes(pos_old, pos_new, color = chr_diff)) +
    geom_point(size = 1, alpha = 1) +
    #theme_simple(axis_lines = TRUE, grid_lines = TRUE, base_size = 10) +
    #scale_color_brewer(type = "qual", palette = "Set2") +
    #scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
    scale_color_discrete(name = "Mapped to different Chr.?") +
    facet_wrap(~chr_new) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("LD chip positions on oar3.1 (old) plotted against ramb v1.0 (new)",
            subtitle = "And compared to the new mapping position by the consortium (cons)") 
p

ggsave("figs/ld_new_vs_old_vs_consnew.jpg", p, width = 11, height = 7)


# check imputed data


