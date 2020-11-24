# objective - acquire V2G and clocH4 values for gene-SNP pairs that are considered by PASCAL and MAGMA
# PASCAL considers SNPs in a window +/-50kb from each gene (that falls into a pathway)
# MAGMA considers SNP in a window +/-0kb from each gene (dafault), or +/-10kb (most liberal settign that comutes in <200h)


setwd("Y://projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/various")


# set gene window coordinates - start with gene coordinates as used by Pascal
gc <- read.delim("Y://projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/scripts/PASCAL/PASCAL/resources/annotation/ucsc/ucsc_known_genes_2013-09-03.txt", header = TRUE, stringsAsFactors = FALSE)

gc.id <- unique(gc$hg19.knownToLocusLink.value)
gc.coord <- lapply(gc.id, function(x){
  gs <- subset(gc, hg19.knownToLocusLink.value == x)
  gs.start <- min(gs$hg19.knownGene.txStart)
  gs.end <- max(gs$hg19.knownGene.txEnd)
  c(gs$hg19.knownGene.chrom[[1]], gs.start, gs.end, x, gs$hg19.kgXref.geneSymbol[[1]])
})
gc.coord.df <- as.data.frame(do.call("rbind", gc.coord))
names(gc.coord.df) <- c("chr", "start", "end", "EntrezID", "GeneName")


# get coordinates for windows
gc.coord.df.50kb <- gc.coord.df
gc.coord.df.50kb$start <- sapply(gc.coord.df.50kb$start, function(x){ as.numeric(as.character(x))-50000})
gc.coord.df.50kb$start <- sapply(gc.coord.df.50kb$start, function(x){ if(as.numeric(as.character(x))<0){return(0)}else(return(x))})
gc.coord.df.50kb$end <- sapply(gc.coord.df.50kb$end, function(x){ as.numeric(as.character(x))+50000})
write.table(gc.coord.df.50kb, "gc.coord.df.50kb.bed", quote = FALSE, sep = "\t", row.names = FALSE)

gc.coord.df.10kb <- gc.coord.df
gc.coord.df.10kb$start <- sapply(gc.coord.df.10kb$start, function(x){ as.numeric(as.character(x))-10000})
gc.coord.df.10kb$start <- sapply(gc.coord.df.10kb$start, function(x){ if(as.numeric(as.character(x))<0){return(0)}else(return(x))})
gc.coord.df.10kb$end <- sapply(gc.coord.df.10kb$end, function(x){ as.numeric(as.character(x))+10000})
write.table(gc.coord.df.10kb, "gc.coord.df.10kb.bed", quote = FALSE, sep = "\t", row.names = FALSE)


# get coordinates for SNPs
#snps <- read.delim("Y://projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/scripts/MAGMA/data/g1000_eur.bim", header = TRUE, stringsAsFactors = FALSE)
# awk '{OFS="\t"; print "chr"$1, $4-1, $4, $2}' /GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/scripts/MAGMA/data/g1000_eur.bim > g1000_eur.bed
# did intersetcBed in shell (remove header from bed files)
# intersectBed -wa -wb -a g1000_eur.bed -b gc.coord.df.50kb.bed > gc.coord.df.50kb.g1000.SNPs.txt
# intersectBed -wa -wb -a g1000_eur.bed -b gc.coord.df.10kb.bed > gc.coord.df.10kb.g1000.SNPs.txt



