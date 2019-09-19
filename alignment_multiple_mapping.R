library(Rsamtools)
library(snpStats)
library(tidyverse)
library(processx)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~ load HD chip alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# check where alignment info starts
system("grep -n @PG data/sheep_genome/aligned/sheep_hd_filt.sam")

sheep_hd_sam_org <- read_delim("data/sheep_genome/aligned/sheep_hd_filt.sam", skip = 2641, 
                           delim = "\t", 
                           col_names = c("snp", "flag", "chr", "pos", "mapq",
                                         "cigar", "rnext", "pnext", "tlen",
                                         "seq", "qual", "edit_dist", "mismatch",
                                         "alignment_score", "subopt_alignment_score",
                                         "alt_hits"))

# we want to have the multiple alignments as tidy data in multiple rows
sheep_hd_sam <- sheep_hd_sam_org %>% 
    filter(!is.na(alt_hits)) %>% 
    mutate(alt_hits = str_remove(alt_hits, "XA:Z:")) %>% 
    mutate(alt_hits_list = str_split(alt_hits, ";"))

sheep_hd_sam %>% 
    dplyr::select(snp, chr, pos, cigar, alt_hits_list) %>% 
    pmap(function(snp, chr, pos, cigar, alt_hits_list) {
        first_row <- tibble(snp = snp, chr = chr, pos = pos, cigar = cigar)
        other_rows <- map_df(alt_hits_list, function(x) {
            if (x == "") return()
            splitted <- unlist(str_split(x, ","))
            tibble(snp = snp,
                   chr = splitted[1], 
                   pos = splitted[2], 
                   cigar = splitted[3])
        })
        out <- rbind(first_row, other_rows)
    }
    ) -> all_mults_df

all_mults_df_proc <- all_mults_df %>% 
    bind_rows() %>% 
    mutate(pos = str_remove(pos, "\\-"),
           pos = str_remove(pos, "\\+"))

all_mults_df_proc %>% 
    write_delim(path = "data/multimapped_hd_chip.txt")

sheep_hd_sam_org %>% 
    dplyr::select(snp, chr, pos) %>% 
    write_delim(path = "data/all_mapped_hd_chip.txt")







sheep_hd_bam <- BamFile("data/sheep_genome/aligned/sheep_hd_filt.bam")
#yieldSize(sheep_hd_bam) <- 604000
seqinfo(sheep_hd_bam)
seqs <- as.data.frame(scanBam(sheep_hd_bam)[[1]]$seq)
names(seqs) <- "seq"

scanBam(sheep_hd_bam, param = ScanBamParam(what=scanBamWhat()))
aln_df <- scanBam(sheep_hd_bam) %>% 
    .[[1]] %>% 
    .[1:11] %>% 
    bind_cols() %>% 
    mutate(seq = seqs$seq) %>% 
    # filter not uniquely mapped reads
    #filter(mapq < 30) %>% 
    # which were mapped properly anyway
    filter(flag %in% c(0,16)) %>% 
    group_by(rname) %>% 
    arrange(rname, pos)

#write_delim(aln_df_low, "mult_maps_flanks.txt", " ")

aln_df_multimapped <- scanBam(sheep_hd_bam) %>% 
    .[[1]] %>% 
    .[1:11] %>% 
    bind_cols() %>% 
    mutate(seq = seqs$seq) %>% 
    # filter not uniquely mapped reads
    filter(mapq <= 30) %>% 
    # which were mapped properly anyway
    filter(flag %in% c(0,16)) %>% 
    arrange(rname, pos)

# get optional positions for multimapped probes from sam (indicated by XA)
system("grep XA data/sheep_genome/aligned/sheep_hd_filt.sam > data/sheep_genome/aligned/multimapped_hd.txt")
multimapped_sam <- read_delim("data/sheep_genome/aligned/multimapped_hd.txt", "\t", col_names = FALSE)
multimapped_sam$X16
#qwidth = read width
# table(aln_df$qwidth)
# sort(table(aln_df$cigar), decreasing = TRUE)[1:500]
# plot(sort(table(aln_df$cigar), decreasing = TRUE)[2:1000])

snps_aln <- aln_df_low %>% 
    # filter(cigar == "100M") %>% 
    # mutate(qname = str_remove(qname, "_Seq[12]")) %>% 
    #mutate(complete_align = ifelse(cigar == "100M", 1, 0)) %>% 
    #group_by(qname) %>% 
    #tally()
    #mutate(complete_align_mean = mean(complete_align)) %>% 
    # check sequences where at least one is perfectly aligned
    #filter(complete_align >= 0.5) %>% 
    .$qname

sum(snps_hd %in% snps_aln)

#~~~~~~~~~~~~~~~~~ load second HD chip alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

sheep_hd_bam2 <- BamFile("data/sheep_genome/aligned/sheep_hd_ext_oneseq_filt.bam")

aln_df2 <- scanBam(sheep_hd_bam2) %>% 
    .[[1]] %>% 
    .[1:11] %>% 
    bind_cols() %>% 
    mutate(snps_name = str_replace(qname, "_Seq[12]", "")) %>% 
    # filter not uniquely mapped reads
    filter(mapq <= 30) %>% 
    # which were mapped properly anyway
    filter(flag %in% c(0,16))

new_snps <- aln_df_low$qname[aln_df_low$qname %in% unique(aln_df2$snps_name)]
all_aligned <- c(aln_df$qname, new_snps)

sum(snps_hd %in% all_aligned)


#~~~~~~~~~~~~~~~~~~ load ld chip alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sheep_ld_bam <- BamFile("data/sheep_genome/aligned/sheep_ld_filt.bam")
#yieldSize(sheep_ld_bam) <- 604000
seqinfo(sheep_ld_bam)
seqs <- as.data.frame(scanBam(sheep_ld_bam)[[1]]$seq)
names(seqs) <- "seq"

aln_df <- scanBam(sheep_ld_bam) %>% 
    .[[1]] %>% 
    .[1:11] %>% 
    bind_cols() %>% 
    mutate(seq = seqs$seq) %>% 
    # filter not uniquely mapped reads
    #filter(mapq == 0) %>% 
    # filter good quality (probably uniquely mapped) reads
    filter(mapq < 30) %>% 
    #write_delim("data/low_qual_alignments.txt", delim = "\t")
    # 0 mean that no flags are set, i.e. 0 means unpaired read, 
    # successfully mapped to reference and to forward strand (what we want)
    filter(flag %in% c(0,16))
# # flag 256 are secondary (and prob. incomplete) alignments
# filter(flag != 256) %>% 
# # reverse strand alignment and secondary alignment
# filter(flag != 272) %>% 
# # reverse strand alignment
# filter(flag == 16)


#qwidth = read width
# table(aln_df$qwidth)
# sort(table(aln_df$cigar), decreasing = TRUE)[1:500]
# plot(sort(table(aln_df$cigar), decreasing = TRUE)[2:1000])

snps_aln <- aln_df %>% 
    # filter(cigar == "100M") %>% 
    # mutate(qname = str_remove(qname, "_Seq[12]")) %>% 
    #mutate(complete_align = ifelse(cigar == "100M", 1, 0)) %>% 
    #group_by(qname) %>% 
    #tally()
    #mutate(complete_align_mean = mean(complete_align)) %>% 
    # check sequences where at least one is perfectly aligned
    #filter(complete_align >= 0.5) %>% 
    .$qname

sum(snps_ld %in% snps_aln)



