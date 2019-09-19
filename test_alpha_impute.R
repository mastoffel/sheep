library(snpStats)
library(tidyverse)
browseVignettes("snpStats")


# get low density SNP chip data and take only SNPs from chromosome 1
sheep_bed <- "data/SNP_chip/20180209_SoayPlates1-83.bed"
sheep_bim <- "data/SNP_chip/20180209_SoayPlates1-83.bim"
sheep_fam <- "data/SNP_chip/20180209_SoayPlates1-83.fam"
# 
sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)

# sex chromosome is 27
table(sample$map$chromosome)
head(sample$genotypes)

# coerce from raw to numeric
sheep_geno <- as(sample$genotypes, Class = "numeric")
sheep_geno[is.na(sheep_geno)] <- 9

# only take chromosome 1
sheep_geno_ch1 <- sheep_geno[, 1:4399]
sheep_geno_ch1[1:10, 1:10]
sheep_geno_ch1 <- as_tibble(sheep_geno_ch1, rownames = "ID")
rm("sheep_geno")

# get high density SNP chip data 
hd_sheep_bed <- "data/SNP_chip/20140214_SheepHD_QC1_Polym.bed"
hd_sheep_bim <- "data/SNP_chip/20140214_SheepHD_QC1_Polym.bim"
hd_sheep_fam <- "data/SNP_chip/20140214_SheepHD_QC1_Polym.fam"

hd_sample <- read.plink(hd_sheep_bed, hd_sheep_bim, hd_sheep_fam)
# coerce from raw to numeric
hd_sheep_geno <- as(hd_sample$genotypes, Class = "numeric")
hd_sheep_geno[is.na(hd_sheep_geno)] <- 9

# which are chromosome 1
table(hd_sample$map$chromosome)

# only take chromosome 1
hd_sheep_geno_chr1 <- hd_sheep_geno[, 1:44138]
hd_sheep_geno_chr1[1:10, 1:10]
hd_sheep_geno_chr1 <- as_tibble(hd_sheep_geno_chr1, rownames = "ID")

names(sheep_geno_ch1)[which(!(names(sheep_geno_ch1) %in% names(hd_sheep_geno_chr1)))]


# sample just a few individuals from low density
sheep_geno_ch1_few <- sample_n(sheep_geno_ch1, 5)

# create genotype data 
all_sheep_geno <- bind_rows(hd_sheep_geno_chr1, sheep_geno_ch1_few)

# create test dataset for imputation with true genotypes
loci_to_del <- sample(1:ncol(hd_sheep_geno_chr1), 40000)
sheep_to_miss <- sample(1:nrow(hd_sheep_geno_chr1), 10)

hd_sheep_missing_ch1 <- hd_sheep_geno_chr1
hd_sheep_missing_ch1[sheep_to_miss, loci_to_del] <- 9









