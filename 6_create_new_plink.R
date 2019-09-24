library(tidyverse)

hd_snps <- read_delim("data/sheep_genome/aligned/final_probe_alignments/final_snp_sets/snp_set_hd_ram.txt", delim = " ") 
table(hd_snps$chr, useNA = "always")

ld_snps <- read_delim("data/sheep_genome/aligned/final_probe_alignments/final_snp_sets/snp_set_ld_ram.txt", delim = " ")
table(ld_snps$chr, useNA = "always")

# some SNPs on the 50k and HD chip have different names (we only consider the X chromosome here)
# we want them to be the same for imputation
snp_names <- read_delim("data/sheep_genome/susie_ld_mapping_output/All_Chromosome_Common_SNPs_QCed.txt", 
                        delim = "\t") %>% 
                        filter(str_detect(SNP.Name.HD, "OARX")) %>% 
                        dplyr::select(SNP.Name.50K, SNP.Name.HD) %>% 
    write_delim("data/sheep_genome/aligned/final_probe_alignments/new_ids_for_x_snps_on_50K.txt", " ",
                col_names = FALSE)

# produce new HD SNP chip plink file
new_map_hd <- hd_snps %>% 
                dplyr::select(snp, pos) %>% 
                write_delim("data/sheep_genome/aligned/final_probe_alignments/new_hd_map_for_plink.txt", " ",
                            col_names = FALSE)

new_chr_hd <- hd_snps %>% 
                    dplyr::select(snp, chr) %>% 
                    write_delim("data/sheep_genome/aligned/final_probe_alignments/new_hd_chr_for_plink.txt", " ",
                                col_names = FALSE)
snps_hd <- hd_snps %>% 
            dplyr::select(snp)%>% 
            write_delim("data/sheep_genome/aligned/final_probe_alignments/new_hd_snps.txt", " ",
                col_names = FALSE)

system(paste0("~/programs/plink --bfile data/SNP_chip/Plates_1-2_HD_QC2 --sheep ",      # --out output/ROH/roh_nofilt
              "--make-bed --out data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram ",
              "--update-map data/sheep_genome/aligned/final_probe_alignments/new_hd_map_for_plink.txt ",
              "--update-chr data/sheep_genome/aligned/final_probe_alignments/new_hd_chr_for_plink.txt ",
              "--extract data/sheep_genome/aligned/final_probe_alignments/new_hd_snps.txt"))

# produce new LD SNP chip plink file
new_map_ld <- ld_snps %>% 
    dplyr::select(snp, pos) %>% 
    write_delim("data/sheep_genome/aligned/final_probe_alignments/new_ld_map_for_plink.txt", " ",
                col_names = FALSE)

new_chr_ld <- ld_snps %>% 
    dplyr::select(snp, chr) %>% 
    write_delim("data/sheep_genome/aligned/final_probe_alignments/new_ld_chr_for_plink.txt", " ",
                col_names = FALSE)

snps_ld <- ld_snps %>% 
    dplyr::select(snp)%>% 
    write_delim("data/sheep_genome/aligned/final_probe_alignments/new_ld_snps.txt", " ",
                col_names = FALSE)

system(paste0("~/programs/plink --bfile data/SNP_chip/Plates_1to87_QC3 --sheep ",      # --out output/ROH/roh_nofilt
              "--make-bed --out data/SNP_chip/ramb_mapping/Plates_1to87_QC4_ram ",
              "--update-map data/sheep_genome/aligned/final_probe_alignments/new_ld_map_for_plink.txt ",
              "--update-chr data/sheep_genome/aligned/final_probe_alignments/new_ld_chr_for_plink.txt ",
              "--extract data/sheep_genome/aligned/final_probe_alignments/new_ld_snps.txt "))

# change IDs of all
system(paste0("~/programs/plink --bfile data/SNP_chip/ramb_mapping/Plates_1to87_QC4_ram --sheep ",      # --out output/ROH/roh_nofilt
              "--make-bed --out data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram ",
              "--update-name data/sheep_genome/aligned/final_probe_alignments/new_ids_for_x_snps_on_50K.txt"))
