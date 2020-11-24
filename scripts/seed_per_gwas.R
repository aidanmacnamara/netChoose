library(tidyverse)


setwd("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/")

# import list of GWAS used
gwas.list <- read.delim("scripts/PASCAL/gwas_109traits.txt", header=FALSE, stringsAsFactors=FALSE, sep="\n")$V1
gwas.list.ed <- gsub(pattern = "gwas_results_", "", gwas.list)
# import MEGAtable
megatable <- read.delim(file="gwas/mega_tbl.txt", header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)


lapply(gwas.list.ed, function(gwasid){
  gwas.hits <- filter(megatable, gwas == gwasid, high_confidence_genetic == 1)
  #if(nrow(gwas.hits)!=0){
    write.table(gwas.hits$entrezgene, file = paste("./analysis/hit_genes/gwas_results_", gwasid, "_hits.txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  #}
})
