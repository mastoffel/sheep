library(snpStats)
plink_geno_path <- "data/SNP_chip/ramb_mapping/Plates_1to87_QC4_ram"
# read merged plink data
sheep_bed <- paste0(plink_geno_path, ".bed")
sheep_bim <- paste0(plink_geno_path, ".bim")
sheep_fam <- paste0(plink_geno_path, ".fam")
full_sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
summary(full_sample$genotypes)
snps_map <- full_sample$map$snp.name
snps_geno <- colnames(full_sample$genotypes@.Data)
sum(snps_geno == snps_map)

# are all positions increasing?
out <- full_sample$map %>% 
    group_by(chromosome) %>% 
    summarise(increasing = all(diff(position) >= 0))
out
