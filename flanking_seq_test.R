# check flanking sequences
library(seqinr)
library(tidyverse)
library(data.table)
library(snpStats)

flanks <- read.fasta("data/sheep_genome/20190913_Flanking_Sequences/20190912_HD_SNPs.fasta")
snps_flank <- map_chr(flanks, function(x) attr(x, "Annot")) %>% 
                str_remove("> ") %>% 
                str_remove("__Seq1") %>% 
                str_remove("__Seq2") 
snps_flank_hd <- snps_flank[!duplicated(snps_flank)]


flanks <- read_csv("data/sheep_genome/flanking_seqs/ovinesnp50-genome-assembly-oar-v3-1.csv", skip = 7)         

# get low density SNP chip data and take only SNPs from chromosome 1
sheep_bed <- "data/SNP_chip/Plates_1-2_HD_QC2.bed"
sheep_bim <- "data/SNP_chip/Plates_1-2_HD_QC2.bim"
sheep_fam <- "data/SNP_chip/Plates_1-2_HD_QC2.fam"
# 
sample_hd <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_hd <- sample_hd$map$snp.name

# ld
# get low density SNP chip data and take only SNPs from chromosome 1
sheep_bed <- "data/SNP_chip/Plates_1to87_QC3.bed"
sheep_bim <- "data/SNP_chip/Plates_1to87_QC3.bim"
sheep_fam <- "data/SNP_chip/Plates_1to87_QC3.fam"
# 
sample_ld <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_ld <- sample_ld$map$snp.name

all_chip <- sample_ld$map

sum(is.na(flanks$SNP))

full <- sample_ld$map %>% 
    filter(!(is.na(chromosome) | is.na(position))) %>% 
    left_join(flanks, by = c("chromosome", "position")) 


sum(duplicated(full$snp.name))

full[duplicated(full$snp.name), ]

# are all hd snps in the flanking file?
sum((snps_hd %in% snps_flank))

# which are not in flanking file
not_flank <- sample_hd$map[!(snps_hd %in% snps_flank), ]
sum(is.na(not_flank$chromosome))

snps_hd
sum((snps_ld %in% snps_flank))

out <- table(snps)

sum(out == 2)
hist(out)

str(flanks)
flanks[[1]]

?rbindlist

rbindlist(flanks[1:2])

test <- flanks[[1]]
