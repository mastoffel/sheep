#
# Find Rambouillet SNP Positions
# SEJ, MAS
# Sep 2019
#

library(GenABEL)
library(plyr)
library(tidyverse)

load("../20190711 Soay Sheep 50K Data/Plates_1-2_HD_QC1.GenABEL.RData")

multimap <- read.table("multimapped_hd_chip.txt", header = T, stringsAsFactors = F)
fullmap <- read.table("all_mapped_hd_chip.txt", header = T, stringsAsFactors = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Format Tables                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Remove multi-mapped SNPs from single map SNPs

fullmap <- subset(fullmap, !snp %in% multimap$snp)

#~~ Examine the double aligns and keep those aligning to the assembled genome

table(fullmap$chr) %>% data.frame

multimap <- multimap[grep("CM", multimap$chr),]
multimap$chr <- gsub("SA:Z:", "", multimap$chr)

#~~ Make a table with the true chromosome ID and join to table

testchr <- data.frame(table(multimap$chr))
testchr$NewChr <- gsub(".1", "", testchr$Var1, fixed = T)
testchr$NewChr <- gsub("CM0084", "", testchr$NewChr) %>% as.numeric
testchr$NewChr <- testchr$NewChr - 71

names(testchr) <- c("chr", "freq", "newchr")

multimap <- join(multimap, testchr)
multimap <- subset(multimap, select = c(snp, pos, chr, newchr))

fullmap <- join(fullmap, testchr)
fullmap <- subset(fullmap, select = c(snp, pos, chr, newchr))

fullmap <- arrange(fullmap, newchr, pos)

multimap <- arrange(multimap, snp)

rm(testchr)

#~~ Create fields for the LD profiles

flanksize <- 10

multimap$MeanLD <- NA
multimap$VarLD  <- NA
multimap$LeftSNP <- NA
multimap$RightSNP <- NA

#~~ subset fullmap into maptab which includes all polymorphic SNPs passing QC.

maptab <- subset(fullmap, snp %in% snp.names(soayHD))

#~~ Get rid of multimap SNPs not in the soayHD

multimap <- subset(multimap, snp %in% snp.names(soayHD))

#~~ Get rid of multimap SNPs that only occur once in the table after getting rid
#   of those mapped to unmapped scaffolds

temp <- table(multimap$snp) %>% data.frame %>% subset(Freq > 1)

multimap <- subset(multimap, snp %in% temp$Var1)

rm(temp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Calculate LD                    #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

multigenos <- as.double.gwaa.data(soayHD[,unique(multimap$snp)])

LDhold <- list()

for(h in sort(unique(multimap$newchr))){
  
  print(paste("Running matches on Chromosome", h))
  
  snplines <- which(multimap$newchr == h)
  
  #~~ Pull out the SNPs on that chromosmoe
  
  chrmap <- maptab[which(maptab$newchr == h),]
  chrgenos <- as.double.gwaa.data(soayHD[,chrmap$snp])
  
  for(i in snplines){
    
    xpos <- which(chrmap$pos > multimap$pos[i])[1]
    
    #~~ Get the start and stop positions
    
    if(is.na(xpos)){
      xstart <- (nrow(chrmap) - flanksize)
      xstop <- nrow(chrmap)
      multimap$LeftSNP[i] <- chrmap$snp[xstop]
    } else {
      xstart <- ifelse(xpos - flanksize < 0, 1, xpos - flanksize)
      xstop  <- ifelse((xpos + flanksize - 1) > nrow(chrmap), nrow(chrmap), (xpos + flanksize - 1))
      
      if(xpos > 1) multimap$LeftSNP[i] <- chrmap$snp[xpos-1]
      multimap$RightSNP[i] <- chrmap$snp[xpos]
      
      
    }
    
    x <- chrmap[xstart:xstop,"snp"]
    
    #~~ Make the data frame of genotypes
    
    y <- data.frame(RefSNP = multimap$snp[i],
                    SNP = x,
                    LD = NA, stringsAsFactors = F)
    
    xref <- multigenos[,multimap$snp[i]]
    #~~ Calculate Rho
    
    for(j in 1:nrow(y)){
      
      y$LD[j] <- cor.test(xref, chrgenos[,y$SNP[j]])$estimate
      
    }
    
    y$LD <- ifelse(y$LD < 0, y$LD * -1, y$LD)
    
    LDhold[[(length(LDhold) + 1)]] <- y
    
    multimap$MeanLD[i] <- mean(y$LD)
    multimap$VarLD[i] <- var(y$LD)
    
    rm(xpos, xstart, xstop, x, y, xref, j)
  }
  
  rm(snplines, chrmap, chrgenos, i)  
}

rm(h, multigenos, flanksize)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 3. Look at the LD                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LDhold <- bind_rows(LDhold)

write.table(multimap, "10_Multimapped_SNPs_MeanLD.txt", row.names = F, sep = "\t", quote = F)

write.table(LDhold, "10_Multimapped_SNPs_FullLD.txt", row.names = F, sep = "\t", quote = F)

subset(multimap, snp == "OAR1_157363688.1")

library(seqinr)
chr1 <- seqinr::read.fasta("sequence (3).fasta")

seq1 <- chr1$CM008472.1[(159282631-200):(159282631+200)] %>% paste(collapse = "")

x <- subset(LDhold, RefSNP == "OAR1_157363688.1")


