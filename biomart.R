library("biomaRt")

ensembl_oaries = useEnsembl(biomart="ensembl", dataset=c("oaries_gene_ensembl"))
x <- listAttributes(ensembl_oaries)

y <- getBM(attributes=c('external_gene_name','description'),
           filters = c("chromosome_name", "start", "end"),
           values = list(20, 277e5, 279e5), mart = ensembl_oaries)
