#!/bin/bash

# Merge plink HD and LD sheep datasets, keep individual order as in original datasets
# merge-mode is default so mismatches are set to missing
# if two variants have the same position, try to merge them (fails if allele names arent matching)

plink --bfile 20140214_SheepHD_QC1_Polym --bmerge 20180209_SoayPlates1-83.bed 20180209_SoayPlates1-83.bim 20180209_SoayPlates1-83.fam --make-bed --indiv-sort 0 --sheep --merge-equal-pos --out merged_sheep_geno
