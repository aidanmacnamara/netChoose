library(tidyverse)


library('CBDD')
setwd("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/")


snps.bloodpressure <- read.table('gwas/May_gwas/gwas_results_GSK500KV3lin_Blood_pressure_medication.txt', sep='\t', stringsAsFactors = FALSE)
# need fake columns


##Toy gene regions
regions <- geneRegionsFromBED('extdata/toySamples/genes_pascal.bed')
##Toy gene sets
pathways <- ontologyFromGMT('extdata/toySamples/pathways_pascal.txt')

##Precomputed 'LD; matrices
load('extdata/toySamples/ld_pascal.rda'