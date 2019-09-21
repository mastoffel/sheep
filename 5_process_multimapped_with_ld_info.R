# Process analysis of multimapped SNP probes, to decide which of them to
# use, following 3 steps:
# (1) Remove multiple mapping probes which do not map to a chromosome
# -> if only one mapped probe for a given SNP is left, retain 
# (2) Find SNP probes which map to multiple positions but have the same
# flanking SNPs -> keep one of the mapping positions as it won't matter
# for imputation
# (3) For all other multimapping probes, examine the LD with the surrounding
# ten SNPs on each side. If LD is high (> 0.3) in only one position and
# low (< 0.15) in all others, retain the position where LD is high.

library(tidyverse)
library(snpStats)

# get HD SNP chip data
sheep_bed <- "data/SNP_chip/Plates_1-2_HD_QC2.bed"
sheep_bim <- "data/SNP_chip/Plates_1-2_HD_QC2.bim"
sheep_fam <- "data/SNP_chip/Plates_1-2_HD_QC2.fam"
# 
sample_hd <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_hd <- sample_hd$map %>% 
    #filter(chromosome %in% 1:26) %>% 
    .$snp.name

# simple files
# load file with all mappings (only one row also for multiple alignments)
all_mapped_hd <- read_delim("data/ld_mapping/all_mapped_hd_chip.txt", delim = " ")
# load file with ALL multiple mappings
multimap_hd <- read_delim("data/ld_mapping/multimapped_hd_chip.txt", delim = " ")
# load file with LD infos
lds <- read_delim("data/ld_mapping/10_Multimapped_SNPs_MeanLD.txt", delim = "\t")

# remove all non-chromosome mappings and filter snps where only one alignment is left
snp_filt0.9 <- multimap_hd %>% 
    filter(str_detect(chr, "CM")) %>% 
    group_by(snp) %>% 
    tally() %>% 
    filter(n==1)
snp_filt1 <- multimap_hd %>% 
    filter(str_detect(chr, "CM")) %>% 
    filter(snp %in% snp_filt0.9$snp)
    
# find all probes where left and right snps (and thus MeanLD) 
# are the same, and where it hence doesnt matter which one to take.
# susie filtered out non-chromosome mappings in her file already. 
snp_filt2 <- lds %>% 
    group_by(snp) %>% 
    mutate(all_equal = ifelse(any(MeanLD > 0) & length(table(MeanLD)) == 1, 1, 0)) %>% 
    filter(all_equal == 1) %>% 
    slice(1) %>% 
    mutate(all_equal = NULL)

# filter with the help of LD
snp_filt3 <- lds %>% 
    filter(!(snp %in% snp_filt2$snp)) %>% 
    group_by(snp) %>% 
    mutate(highest_ld = ifelse(
        # check whether highest LD among multiple alignments is high and whether all others are small
        (max(MeanLD) > 0.3) & all(MeanLD[-which.max(MeanLD)] < 0.15), 1, 0
    )) %>% 
    filter(highest_ld == 1) %>% 
    filter(MeanLD == max(MeanLD))

snps_filts <- bind_rows(snp_filt1, snp_filt2, snp_filt3)

# put everthing together
# check where alignment info starts
system("grep -n @PG data/sheep_genome/aligned/sheep_hd_filt.sam")
sheep_hd_sam_org <- read_delim("data/sheep_genome/aligned/sheep_hd_filt.sam", skip = 2641, 
                               delim = "\t", 
                               col_names = c("snp", "flag", "chr", "pos", "mapq",
                                             "cigar", "rnext", "pnext", "tlen",
                                             "seq", "qual", "edit_dist", "mismatch",
                                             "alignment_score", "subopt_alignment_score",
                                             "alt_hits")) %>% 
    filter(flag %in% c(0,16))


final_snp_set <- sheep_hd_sam_org %>% 
    # selecting all multiple hits
    filter(is.na(alt_hits) | !(snp %in% snps_filts$snp)) %>% 
    filter(snp %in% snps_hd)

