# Probe sequences to fasta for alignment against new reference
library(tidyverse)
library(seqinr)
flanks_hd <- read_csv("data/sheep_genome/flanking_seqs/SheepHD_AgResearch_Cons_15041608_A.csv", skip = 7) %>% 
                filter(!is.na(AlleleA_ProbeSeq))
flanks_ld <- read_csv("data/sheep_genome/flanking_seqs/ovinesnp50-genome-assembly-oar-v3-1.csv", skip = 7) %>% 
                filter(!is.na(AlleleA_ProbeSeq))

# some checks
sum(duplicated(flanks_hd$Name))
sum(duplicated(flanks_ld$Name))

# are all probes 50 long?
any(str_count(flanks_hd$AlleleA_ProbeSeq) != 50)
any(str_count(flanks_ld$AlleleA_ProbeSeq) != 50)

# write_to_fasta
flanks_hd_fasta <- rep(NA, 2*nrow(flanks_hd))
flanks_hd_fasta[seq(from = 1, to = length(flanks_hd_fasta), by =2)] <- paste0(">",flanks_hd$Name)
flanks_hd_fasta[seq(from = 2, to = length(flanks_hd_fasta), by =2)] <- flanks_hd$AlleleA_ProbeSeq
write_lines(flanks_hd_fasta, "data/sheep_genome/sheep_hd_flanks.fasta")

# write_to_fasta
flanks_ld_fasta <- rep(NA, 2*nrow(flanks_ld))
flanks_ld_fasta[seq(from = 1, to = length(flanks_ld_fasta), by =2)] <- paste0(">",flanks_ld$Name)
flanks_ld_fasta[seq(from = 2, to = length(flanks_ld_fasta), by =2)] <- flanks_ld$AlleleA_ProbeSeq
write_lines(flanks_ld_fasta, "data/sheep_genome/sheep_ld_flanks.fasta")


# get HD SNP chip data
sheep_bed <- "data/SNP_chip/Plates_1-2_HD_QC2.bed"
sheep_bim <- "data/SNP_chip/Plates_1-2_HD_QC2.bim"
sheep_fam <- "data/SNP_chip/Plates_1-2_HD_QC2.fam"
# 
sample_hd <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_hd <- sample_hd$map %>% 
    #filter(chromosome %in% 1:26) %>% 
    .$snp.name

sum(snps_hd %in% flanks_snps)
