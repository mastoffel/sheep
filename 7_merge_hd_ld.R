# Merge plink HD and LD sheep datasets, keep individual order as in original datasets
# merge-mode is default so mismatches are set to missing
# if two variants have the same position, try to merge them (fails if allele names arent matching)

system(paste0("~/programs/plink --bfile data/SNP_chip/ramb_mapping/Plates_1-2_HD_QC3_ram --sheep ",      # --out output/ROH/roh_nofilt
              "--bmerge data/SNP_chip/ramb_mapping/Plates_1to87_QC5_ram ",
              "--make-bed --out data/SNP_chip/ramb_mapping/merged_sheep_geno_ram ",
              "--indiv-sort 0 --merge-equal-pos"))
# --merge-equal-pos
