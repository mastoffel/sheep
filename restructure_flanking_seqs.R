library(tidyverse)

# extract seq1s
hd_flanks <- readLines("data/sheep_genome/flanking_seqs/20190912_HD_SNPs.fasta")

hd_flanks_seq1 <- which(str_detect(hd_flanks, "Seq1")) 
hd_flanks_seq1 <- sort(c(hd_flanks_seq1, hd_flanks_seq1 + 1))

flanks_df <- tibble(snps = hd_flanks[seq(1, length(hd_flanks), by = 2)],
                    seq =  hd_flanks[seq(2, length(hd_flanks), by = 2)]) %>% 
             mutate(snps = str_remove(snps, "__Seq[12]")) %>% 
             mutate(snps = str_remove(snps, "> ")) %>% 
             group_by(snps) %>% 
             summarise(seqs = paste0(seq, collapse="")) %>% 
             mutate(snps = paste0(">",snps))

fasta_file <- as.character(rep(NA, nrow(flanks_df) * 2))
fasta_file[seq(1, length(fasta_file), by = 2)] <- flanks_df$snps
fasta_file[seq(2, length(fasta_file), by = 2)] <- flanks_df$seqs

write_lines(fasta_file, "data/sheep_genome/flanking_seqs/20190912_HD_SNPs_modified.fasta")


# some exploration
# all sequences are only there once
flanks_df <- tibble(snps = hd_flanks[seq(1, length(hd_flanks), by = 2)],
                    seq =  hd_flanks[seq(2, length(hd_flanks), by = 2)]) %>% 
    mutate(snps = str_remove(snps, "__Seq[12]")) %>% 
    mutate(snps = str_remove(snps, "> ")) %>% 
    group_by(snps) %>% 
    tally()
    
table(flanks_df$n)
