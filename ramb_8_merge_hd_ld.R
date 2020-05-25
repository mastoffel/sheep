# Merge plink HD and LD sheep datasets, keep individual order as in original datasets
# merge-mode is default so mismatches are set to missing
# if two variants have the same position, try to merge them (fails if allele names arent matching)


# merge chips based on rambouillet v1 mapping
system(paste0("/usr/local/bin/plink --bfile data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram --sheep ",      # --out output/ROH/roh_nofilt
              "--bmerge data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram.bed data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram.bim data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram.fam ",
              "--make-bed --out data/SNP_chip/ramb_mapping/merged_sheep_geno_ram ",
              "--indiv-sort 0 --merge-equal-pos "))

# merge chips based on rambouillet v1 mapping and discard wrongly mapped SNPs
system(paste0("/usr/local/bin/plink --bfile data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram --sheep ",      # --out output/ROH/roh_nofilt
              "--bmerge data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram.bed data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram.bim data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram.fam ",
              "--make-bed --out data/SNP_chip/ramb_mapping/merged_sheep_geno_ram_snpfilt ",
              "--indiv-sort 0 --exclude ../sheep_ID/data/wrongly_mapped_hdld_plink.txt --merge-equal-pos "))
# --merge-equal-pos

# merge chips based on oar 3.1 mapping
system(paste0("/usr/local/bin/plink --bfile data/SNP_chip/oar31_mapping/Plates_1-2_HD_QC2 --sheep ",      # --out output/ROH/roh_nofilt
              "--bmerge data/SNP_chip/oar31_mapping/Plates_1to87_QC3.bed data/SNP_chip/oar31_mapping/Plates_1to87_QC3.bim data/SNP_chip/oar31_mapping/Plates_1to87_QC3.fam ",
              "--make-bed --out data/SNP_chip/oar31_mapping/merged_sheep_geno_oar31 ",
              "--indiv-sort 0 --merge-equal-pos "))

