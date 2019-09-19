library(snpStats)
library(tidyverse)
browseVignettes("snpStats")


# get low density SNP chip data and take only SNPs from chromosome 1
sheep_bed <- "data/SNP_chip/20180209_SoayPlates1-83.bed"
sheep_bim <- "data/SNP_chip/20180209_SoayPlates1-83.bim"
sheep_fam <- "data/SNP_chip/20180209_SoayPlates1-83.fam"

sample <- read.plink(sheep_bed, sheep_bim, sheep_fam)
sample$genotypes
hist(col.summary(sample$genotypes)$MAF)
head(sample$fam)
head(sample$map)

# sex chromosome is 27
table(sample$map$chromosome)

# coerce from raw to numeric
sheep_geno <- as(sample$genotypes, Class = "numeric")
sheep_geno[is.na(sheep_geno)] <- 9

# only take chromosome 1
sheep_geno_ch1 <- sheep_geno[, 1:4399]
sheep_geno_ch1[1:10, 1:10]
sheep_geno_ch1 <- as_tibble(sheep_geno_ch1, rownames = "ID")


# get high density SNP chip data 








dim(sheep_geno)

sheep_geno[1:10, 1:10]


rownames(sample$genotypes@.Data)

dat <- as_tibble(sample$genotypes@.Data)


gtypes <- as.matrix(sample$genotypes@.Data)
is.raw(gtypes)
class(gtypes_sub) <- "SnpMatrix"
gtypes_sub
write.SnpMatrix(gtypes_sub, file = "test.txt")
snps <- read_delim("test.txt", delim = " ")

gtypes_sub <- gtypes[1:100, 1:100]
packBits(gtypes_sub[,1], type = "integer")

str(gtypes_sub[,1])

rawToChar(gtypes[1000:1100, 1000:1100])

rawToChar(gtypes[1,1])

str(gtypes)

str(matrix(c(2,2,2)))
sample

ld_sheep <- ld(sample$genotypes, stats = c("D.prime", "R.squared"), depth = 100)

image(ld_sheep$D.prime, lwd = 0)

lds <- quantile(ld_sheep$D.prime@x, na.rm = TRUE)
lds
options(scipen = 999)
lds


