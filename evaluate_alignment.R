library(Rsamtools)
library(snpStats)
library(tidyverse)

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
     filter(chromosome %in% 1:26) %>% 
    .$snp.name

#~~~~~~~~~~~~~~~~~~~~~~~~~ load HD chip alignment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
sheep_hd_bam <- BamFile("data/sheep_genome/aligned/sheep_hd_filt.bam")
#yieldSize(sheep_hd_bam) <- 604000
seqinfo(sheep_hd_bam)
seqs <- as.data.frame(scanBam(sheep_hd_bam)[[1]]$seq)
names(seqs) <- "seq"

aln_df <- scanBam(sheep_hd_bam) %>% 
        .[[1]] %>% 
        .[1:11] %>% 
        bind_cols() %>% 
        mutate(seq = seqs$seq) %>% 
        # filter not uniquely mapped reads
        #filter(mapq == 0) %>% 
        # filter good quality (probably uniquely mapped) reads
        filter(mapq >= 30) %>% 
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


# are any reads mapped to more than one position?
# a few a mapped multiple times
check_dup <- aln_df %>% 
    group_by(qname) %>% 
    tally() %>% 
    filter(n >= 2) %>% 
    .$qname

aln_df_dup <- aln_df %>% 
    filter(qname %in% check_dup)
aln_df_dup 

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

sum(snps_hd %in% snps_aln)


# how do the new alignments compare to the old ones?
aln_df_comp <- aln_df %>% 
    filter(str_detect(rname, "CM")) %>% 
    mutate(chromosome = as.numeric(str_sub(rname, 7,8))-71) 

# old alignments
aln_df_old <- read_delim("data/sheep_genome/susie_extracted_flanks/20190912_HD_SNP_Positions.txt", delim = "\t") %>% 
    rename(qname = Name)

aln_full <- aln_df_comp %>% 
                left_join(aln_df_old) %>% 
                rename(old_position = Position,
                       new_position = pos) %>% 
                filter(old_position < 1e+09)

ggplot(aln_full, aes(old_position, new_position)) + geom_point() + theme_classic()

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
    filter(mapq >= 30) %>% 
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


# are any reads mapped to more than one position?
# a few a mapped multiple times
check_dup <- aln_df %>% 
    group_by(qname) %>% 
    tally() %>% 
    filter(n >= 2) %>% 
    .$qname

aln_df_dup <- aln_df %>% 
    filter(qname %in% check_dup)
aln_df_dup 

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



