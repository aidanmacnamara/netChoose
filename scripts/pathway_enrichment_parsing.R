library(tidyverse)


setwd("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/")

# import list of GWAS used
gwas.list <- read.delim("scripts/PASCAL/gwas_109traits.txt", header=FALSE, stringsAsFactors=FALSE, sep="\n")$V1
