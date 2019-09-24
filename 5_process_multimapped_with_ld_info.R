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
library("GenomicAlignments")

# get HD SNP chip data
sheep_bed <- "data/SNP_chip/Plates_1-2_HD_QC2.bed"
sheep_bim <- "data/SNP_chip/Plates_1-2_HD_QC2.bim"
sheep_fam <- "data/SNP_chip/Plates_1-2_HD_QC2.fam"
# 
sample_hd <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_hd <- sample_hd$map %>% 
     #filter(chromosome %in% 1:26) %>% 
    .$snp.name

# get LD SNP chip data
sheep_bed <- "data/SNP_chip/Plates_1to87_QC3.bed"
sheep_bim <- "data/SNP_chip/Plates_1to87_QC3.bim"
sheep_fam <- "data/SNP_chip/Plates_1to87_QC3.fam"
#
sample_ld <- read.plink(sheep_bed, sheep_bim, sheep_fam)
snps_ld <- sample_ld$map %>%
    #filter(chromosome %in% 1:26) %>%
    .$snp.name

# SNP chip type
chip <- "ld" # "hd" "ld"

if (chip == "ld") {
    snps_chip <- snps_ld
    sample_chip <- sample_ld
} else {
    snps_chip <- snps_hd
    sample_chip <- sample_hd
}

# simple files
# load file with all mappings (only one row also for multiple alignments)
all_mapped <- read_delim(paste0("data/sheep_genome/aligned/all_mapped_", chip ,"_chip.txt"), delim = " ")
# load file with ALL multiple mappings
multimap <- read_delim(paste0("data/sheep_genome/aligned/multimapped_", chip, "_chip.txt"), delim = " ") 
# load file with LD infos >>>> check that its the right file <<<<<<<<<<<<<<<
lds <- read_delim("data/sheep_genome/susie_ld_mapping_output/10_Multimapped_SNPs_MeanLD_50K.txt", delim = "\t")

# 2 duplicated rows here
#dup_rows <- duplicated(data.table::as.data.table(lds))
#lds <- lds[!dup_rows, ]

# (1) remove all non-chromosome mappings and get snps where only one alignment is left
snp_filt1 <- multimap %>% 
    #distinct() %>% 
    filter(str_detect(chr, "^CM")) %>% 
    group_by(snp) %>% 
    filter(n() == 1) %>% 
    filter(!(snp == "oar3_OAR15_14252885")) # this snp makes problems

# (2) find all probes where left and right snps 
# are the same, and where it hence doesnt matter which one to take (so take one of them)
# this filter only collapses the similar alignments, so a probe might still have
# multiple alignments
# susie filtered out non-chromosome mappings in her file already. 
lds2 <- lds %>% 
    # remove 2 rows which were doubled
    #distinct() %>% 
    filter(str_detect(chr, "^CM")) %>% 
    # add flag and cigar fields
    left_join(multimap, by = c("snp", "chr", "pos")) %>% 
    filter(!is.na(cigar)) %>% 
    # and filter out the snps which are 
    mutate(flanking_snps = paste0(LeftSNP, "_", RightSNP)) %>% 
    group_by(flanking_snps) %>%
    # arrange by position to always slice the same bit (the smaller pos)
    dplyr::arrange(chr, pos, MeanLD) %>% 
    dplyr::slice(1) 

# snps with only one mapping left after collapsing
snp_filt2 <- lds2 %>% 
    ungroup() %>% 
    group_by(snp) %>% 
    # get snp
    filter(n() == 1)

# how to filter LD in the rest?
ld_pattern <- lds2 %>% 
    filter(!(snp %in% snp_filt2$snp)) %>% 
    arrange(snp) %>% 
    group_by(snp) %>% 
    # find difference between highest and second highest vector
    summarise(max_ld = max(MeanLD), second_ld = sort(MeanLD, TRUE)[2]) %>% 
    mutate(diff = max_ld - second_ld)
hist(ld_pattern$diff, breaks = 100)
abline(v = 0.1)

# filter the rest with the help of LD
snp_filt3 <- lds2 %>% 
    filter(!(snp %in% snp_filt2$snp)) %>% 
   # arrange(snp)
    group_by(snp) %>% 
    mutate(highest_ld = ifelse(
        # check whether highest LD among multiple alignments is high and whether all others are small (hd?)
        #(max(MeanLD) > 0.3) & all(MeanLD[-which.max(MeanLD)] < 0.15), 1, 0
        # if highest LD is at least r2 = 0.1 higher than others, take that
        (max(MeanLD) - sort(MeanLD, decreasing = TRUE)[2]) >= 0.1, 1, 0)) %>%  #  MeanLD[-which.max(MeanLD)]
    filter(highest_ld == 1) %>%
    filter(MeanLD == max(MeanLD))

# sanity check
check_snps <- snp_filt3 %>% 
    arrange(snp)
sum(duplicated(check_snps$snp))

# bind all snps that we want to keep together
snps_filts <- bind_rows(snp_filt1, snp_filt2, snp_filt3) 
sum(duplicated(snps_filts$snp))

# put everthing together
# check where alignment info starts
num_to_skip <- system(paste0("grep -n @PG data/sheep_genome/aligned/sheep_", chip, "_filt.sam"), intern = TRUE) %>% 
           str_split(":") %>% 
           .[[1]] %>% 
           .[1] %>% 
           as.numeric()

sheep_sam_org <- read_delim(paste0("data/sheep_genome/aligned/sheep_", chip, "_filt.sam"), skip = num_to_skip, 
                               delim = "\t", 
                               col_names = c("snp", "flag", "chr", "pos", "mapq",
                                             "cigar", "rnext", "pnext", "tlen",
                                             "seq", "qual", "edit_dist", "mismatch",
                                             "alignment_score", "subopt_alignment_score",
                                             "alt_hits")) %>% 
                    filter(flag %in% c(0,16)) #%>% 
                    #filter(mapq >= 10)


# add mapq values and select relevant columns
snps_filts_final <- snps_filts %>% 
    left_join(sheep_sam_org[c("snp", "mapq")], by = "snp") %>% 
    dplyr::select(snp, chr, pos, cigar, flag, mapq) %>% 
    # set mapq value high
    mutate(mapq = 60)

# sheep_sam_org contains all mapped snps, and multimapped snps appear 
# only once with one of the mapped positions
final_snp_set <- sheep_sam_org %>% 
    dplyr::select(snp, chr, pos, cigar, flag, mapq) %>% 
    # remove snps which we have filtered 
    filter(!(snp %in% snps_filts_final$snp)) %>% 
    # and replace with correct positions
    bind_rows(snps_filts_final) %>% 
    # filter everything not on chromosomes
    filter(str_detect(chr, "^CM")) %>% 
    # filter everything with low mapq value and no multiple mappings
    filter(mapq > 30) %>% 
    # if aligned to forward strand, add length of alignment to reference, 
    # if aligned to reverse strand, substract 1
    mutate(pos = ifelse(flag == 0, pos+extractAlignmentRangesOnReference(cigar)[[1]]@width, pos-1)) %>% 
    # retain only snps on the chip
    filter(snp %in% snps_chip) 

# adjustment to sex chromosome snp position (from plink error)
if (chip == "hd") {
    final_snp_set[final_snp_set$snp == "oar3_OARX_58642295", "pos"] <- final_snp_set[final_snp_set$snp == "oar3_OARX_58642295", "pos"] + 1
    
}

write_delim(final_snp_set, paste0("data/sheep_genome/aligned/final_probe_alignments/snp_set_", chip, ".txt"),
            delim = " ")

# sanity check 
sum(duplicated(final_snp_set$snp))

#~~~~~~~~~~~~~~~~~~~~ plot old against new ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
snps_old <- sample_chip$map %>% 
                dplyr::rename(snp = snp.name,
                       chr = chromosome,
                       pos_old = position)

snps_new <- final_snp_set %>% 
                    mutate(chr = str_sub(chr, 7, 8)) %>% 
                    mutate(chr = as.numeric(chr) - 71) %>% 
                    dplyr::rename(pos_new = pos)

snps_new[snps_new$chr == -71, ]

snps_oldnew <- snps_new %>% 
                    left_join(snps_old, by = c("snp", "chr")) %>% 
                    dplyr::select(snp, chr, pos_new, pos_old) %>% 
                    mutate(pos_old = ifelse(is.na(pos_old), -1, pos_old),
                           pos_new = ifelse(is.na(pos_new), -1, pos_new))

library(viridis)

my_cols <- c("#ffd16a",
            "#8c0fcd",
            "#30ba00",
            "#486dff",
            "#aec700",
            "#e07cff",
            "#bbff71",
            "#4f006e",
            "#00e788",
            "#f30032",
            "#02f0d6",
            "#ca3400",
            "#4de2ff",
            "#e1a800",
            "#0f001e",
            "#dcfeff",
            "#63001d",
            "#2b6900",
            "#ff80ad",
            "#00281e",
            "#ff9156",
            "#019cbb",
            "#812d00",
            "#007866",
            "#750048",
            "#ffcab2",
            "#363400")
snps_oldnew %>% 
    mutate(chr = as.factor(chr)) %>% 
    sample_frac(1) %>% 
ggplot(aes(pos_new, pos_old, fill = chr)) + geom_point(shape = 21, size = 2, alpha = 1, lwd = 0.1) +
    scale_fill_manual(values = my_cols) +
    theme_bw() -> newold_pos

ggsave(paste0("figs/", chip, "_snps_new_vs_old.jpg"), height = 5, width = 7)

# plot new against new from Rudi #
rudi_aligned <- read_csv("data/sheep_genome/aligned_by_rudi/SheepHD_AgResearch_Cons_15041608_A.csv.rambo_pos_20190513.csv") %>% 
    dplyr::rename(chr = chrom, snp = entry, pos_rudi = pos)
rudi_aligned <- read_csv("data/sheep_genome/aligned_by_rudi/OvineSNP50_B.csv.rambo_pos_20190513.csv") %>% 
    dplyr::rename(chr = chrom, snp = entry, pos_rudi = pos)

# snps_new %>% 
#     left_join(rudi_aligned, by = c("snp", "chr")) %>% 
#     mutate(pos_rudi = ifelse(is.na(pos_rudi), -1, pos_rudi)) %>% 
#     ggplot(aes(pos_new, pos_rudi)) + geom_point() -> p1
# 
# ggsave(paste0("figs/", chip, "_snps_new_vs_rudi.jpg"), height = 5)

# alignment logic:
# aligned on forward strand: SNP position = pos + 50
# aligned on reverse strand: SNP position = pos - 1

snps_new %>% 
    left_join(rudi_aligned, by = c("snp", "chr")) %>% 
    mutate(pos_rudi = ifelse(is.na(pos_rudi), -1, pos_rudi)) %>% 
    #mutate(pos_new = pos_new + 50) %>% 
    # if mapped to forward strand, add 50 bp to get to snp
    mutate(pos_new = ifelse(flag == 0, pos_new + 50, pos_new)) %>% 
    ggplot(aes(pos_new, pos_rudi)) + geom_point(size = 0.1) -> p1
    # ggplot(aes(pos_new, pos_rudi)) + geom_point(size = 1) +
    # scale_x_continuous(limits = c(10000000, 10005000)) +
    # scale_y_continuous(limits =  c(10000000, 10005000)) +
    # geom_line(data = tibble(a=1e7:10005000, b=1e7:10005000), mapping = aes(a,b))

ggsave(paste0("figs/", chip, "_snps_new_vs_rudi.jpg"), height = 5)
