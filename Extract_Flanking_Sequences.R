

x <- read.table("C:/Users/sjohns10/Google Drive/Soay Sheep Genomic Data/20190704 GenomeStudio Project Data/Plates 1-2 HD/Plates 1-2 HD/E10905_SheepHD_Plates1-2_310114_FullSNPTable.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
head(x)
x <- x[,2:4]

write.table(x, "20190912_HD_SNP_Positions.txt", row.names = F, sep = "\t", quote = F)

x$Seq1 <- NA
x$Seq2 <- NA

table(x$Chr)

library(seqinr)

library(dplyr)


for(i in c(1:26, "X")){
  
  print(paste("Loading Chromosome", i))
  
  chr <- read.fasta(paste0("Ovis_aries.Oar_v3.1.dna.chromosome.", i, ".fa"))
  chr <- chr[[1]]
  
  linevec <- which(x$Chr == i)
  
  for(j in 1:length(linevec)){
    
    if(j %in% seq(1, length(linevec), 1000)) print(paste("Extracting SNP", j, "of", length(linevec)))
      
    x$Seq1[linevec[j]] <- chr[(x$Position[linevec[j]]-100):(x$Position[linevec[j]]-1)] %>% paste(collapse = "") %>% toupper
    x$Seq2[linevec[j]] <- chr[(x$Position[linevec[j]]+1):(x$Position[linevec[j]]+100)] %>% paste(collapse = "") %>% toupper
    
  }
  
  rm(chr, linevec)
  gc()
  
 # save(x, file = "temp.RData")
  
}

