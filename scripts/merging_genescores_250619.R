library(tidyverse)


setwd("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/")

# import list of GWAS used
gwas.list <- read.delim("scripts/PASCAL/gwas_109traits.txt", header=FALSE, stringsAsFactors=FALSE, sep="\n")$V1



# import genescores
# create and name file paths
# PASCAL
gwas.list.files.metabase <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/Metabase_GO_Maps_entrez.filt10_sum/", x, ".sum.genescores.txt", sep="")
})
gwas.list.files.ppi <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/DSDSpectral.consensus_top50.1_ppi_anonym_v2_Entrez_sum/", x, ".sum.genescores.txt", sep="")
})
gwas.list.files.coexp <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/DSDSpectral.consensus_top50.4_coexpr_anonym_v2_Entrez_sum/", x, ".sum.genescores.txt", sep="")
})
gwas.list.files.reactome <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/c2.cp.reactome.v6.2.entrez_sum/", x, ".sum.genescores.txt", sep="")
})
gwas.list.files.metabase_noghit <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/Metabase_GO_Maps_entrez.filt10_sum_noghit/", x, ".sum.genescores.txt", sep="")
})
# MAGMA
gwas.list.files.magma <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/magma/magma_109traits_50kb/", x, ".genes.out.tab", sep="")
})

# read in datasets
gwas.allFiles.metabase <- lapply(gwas.list.files.metabase, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE, colClasses = c("character", "integer", "integer", "character", "integer", "character", "numeric", "character"))})
gwas.allFiles.ppi <- lapply(gwas.list.files.ppi, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE, colClasses = c("character", "integer", "integer", "character", "integer", "character", "numeric", "character"))})
gwas.allFiles.coexp <- lapply(gwas.list.files.coexp, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE, colClasses = c("character", "integer", "integer", "character", "integer", "character", "numeric", "character"))})
gwas.allFiles.reactome <- lapply(gwas.list.files.reactome, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE, colClasses = c("character", "integer", "integer", "character", "integer", "character", "numeric", "character"))})
gwas.allFiles.metabase_noghit <- lapply(gwas.list.files.metabase_noghit, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE, colClasses = c("character", "integer", "integer", "character", "integer", "character", "numeric", "character"))})
gwas.allFiles.magma <- lapply(gwas.list.files.magma, function(x){read.delim(x, header=FALSE, row.names = NULL, sep="\t", stringsAsFactors = FALSE, colClasses = c("integer", "integer", "integer", "integer", "integer", "numeric", "numeric"))})

names(gwas.allFiles.metabase) <- gwas.list
names(gwas.allFiles.ppi) <- gwas.list
names(gwas.allFiles.coexp) <- gwas.list
names(gwas.allFiles.reactome) <- gwas.list
names(gwas.allFiles.metabase_noghit) <- gwas.list
names(gwas.allFiles.magma) <- gwas.list

# adjust colnames for magma imports
gwas.allFiles.magma <- lapply(gwas.allFiles.magma, setNames, nm = c("gene_id", "chromosome", "nSNPs", "size", "N", "stat", "pvalue"))

# create master list
gwas.allFiles <- list(gwas.allFiles.metabase, gwas.allFiles.ppi, gwas.allFiles.coexp, gwas.allFiles.reactome, gwas.allFiles.metabase_noghit, gwas.allFiles.magma)





#####################################################################################################################################
###   Comparing GENE SCORES                 
#####################################################################################################################################

## TYDING UP DATASETS
# 1st, for every runtype, extract only the pvalues column
gwas.allFiles.pvals <- lapply(gwas.allFiles, function(x){
  lapply(x, function(y){y[,c("gene_id", "pvalue")]})
})

#merge genescores across GWAS
gwas.allFiles.pvals.tidy <- lapply(gwas.allFiles.pvals, function(x){
  Reduce(function(...) merge(..., by="gene_id"), x)
})

# assign names
gwas.allFiles.pvals.tidy <- lapply(gwas.allFiles.pvals.tidy, function(x){setNames(x, nm = c("gene_id", gwas.list))})
gwas.allFiles.names <- list("Pascal_Metabase", "Pascal_PPI", "Pascal_coExp", "Pascal_Reactome", "Pascal_Metabase_noGenHit", "Magma_Metabase")
gwas.allFiles.pvals.tidyed <- list()
for(i in 1:length(gwas.allFiles.names)){
  gwas.allFiles.pvals.tidyed[[i]] <- gather(gwas.allFiles.pvals.tidy[[i]], key=gwas, value=tmp, gwas_results_GSK500kV3_BMD_Combined:gwas_results_GSK500kV3_Whole_body_fat_mass, na.rm = FALSE, convert = FALSE)
  names(gwas.allFiles.pvals.tidyed[[i]]) <- c("gene_id", "gwas", gwas.allFiles.names[[i]])
}

# merge across runtypes
gwas.allFiles.pvals.tidyed.all <- Reduce(function(...) full_join(..., by=c("gene_id"="gene_id","gwas"="gwas")), gwas.allFiles.pvals.tidyed)
gwas.allFiles.pvals.tidyed.all$gwas <- gsub("gwas_results_", "", gwas.allFiles.pvals.tidyed.all$gwas)

# MERGE with MEGATABLE - Pascal & Magma Metabase
megatable <- read.delim(file="gwas/mega_tbl.txt", header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)
megatable.updated <- left_join(megatable, gwas.allFiles.pvals.tidyed.all[,c("gene_id", "gwas", "Pascal_Metabase", "Magma_Metabase")], by=c("entrezgene"="gene_id","gwas"="gwas"))
#write.table(megatable.updated, "gwas/mega_tbl_250619.txt", quote = FALSE, sep = "\t", row.names = FALSE)






############################################
# parsing genescores according to gene categories
#load("r_data/gene_buckets.RData")
#load("r_data/gene_buckets_no_seed.RData")

# extract genescores per gwas per category
gene_buckets.mega <- lapply(gene_buckets, function(x){
  lapply(names(x), function(y){
    yg <- x[[y]]
    subset(megatable.updated, entrezgene %in% yg & gwas == y )
  })
})
gene_buckets_no_seed.mega <- lapply(gene_buckets_no_seed, function(x){
  lapply(names(x), function(y){
    yg <- x[[y]]
    subset(megatable.updated, entrezgene %in% yg & gwas == y )
  })
})

gene_buckets.mega.merged.h4 <- lapply(gene_buckets.mega, function(x){as.numeric(do.call("rbind",x)$h4)})
gene_buckets.mega.merged.Pascal <- lapply(gene_buckets.mega, function(x){as.numeric(do.call("rbind",x)$Pascal_Metabase)})
gene_buckets.mega.merged.Magma <- lapply(gene_buckets.mega, function(x){as.numeric(do.call("rbind",x)$Magma_Metabase)})
gene_buckets_no_seed.mega.merged.h4 <- lapply(gene_buckets_no_seed.mega, function(x){as.numeric(do.call("rbind",x)$h4)})
gene_buckets_no_seed.mega.merged.Pascal <- lapply(gene_buckets_no_seed.mega, function(x){as.numeric(do.call("rbind",x)$Pascal_Metabase)})
gene_buckets_no_seed.mega.merged.Magma <- lapply(gene_buckets_no_seed.mega, function(x){as.numeric(do.call("rbind",x)$Magma_Metabase)

#gene_buckets.mega.merged.h4 <- lapply(gene_buckets.mega, function(x){as.numeric(do.call("rbind",x)$h4)})
gene_buckets.mega.merged.Pascal.log10 <- lapply(gene_buckets.mega, function(x){-log10(as.numeric(do.call("rbind",x)$Pascal_Metabase))})
gene_buckets.mega.merged.Magma.log10 <- lapply(gene_buckets.mega, function(x){-log10(as.numeric(do.call("rbind",x)$Magma_Metabase))})
#gene_buckets_no_seed.mega.merged.h4 <- lapply(gene_buckets_no_seed.mega, function(x){as.numeric(do.call("rbind",x)$h4)})
gene_buckets_no_seed.mega.merged.Pascal.log10 <- lapply(gene_buckets_no_seed.mega, function(x){-log10(as.numeric(do.call("rbind",x)$Pascal_Metabase))})
gene_buckets_no_seed.mega.merged.Magma.log10 <- lapply(gene_buckets_no_seed.mega, function(x){-log10(as.numeric(do.call("rbind",x)$Magma_Metabase))})


# plot gene scores
par(mfrow=c(2,3))
par(mar=c(5.1,13.1,1.1,1.1))
boxplot(rev(gene_buckets.mega.merged.Pascal), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Pascal gene score")
boxplot(rev(gene_buckets.mega.merged.Magma), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Magma gene score")
boxplot(rev(gene_buckets.mega.merged.h4), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Colocalisation h4")
boxplot(rev(gene_buckets_no_seed.mega.merged.Pascal), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Pascal gene score")
boxplot(rev(gene_buckets_no_seed.mega.merged.Magma), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Magma gene score")
boxplot(rev(gene_buckets_no_seed.mega.merged.h4), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Colocalisation h4")
dev.off()

# plot gene scores log scale
par(mfrow=c(2,3))
par(mar=c(5.1,13.1,1.1,1.1))
boxplot(rev(gene_buckets.mega.merged.Pascal.log10), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="-log10(Pascal gene score)")
boxplot(rev(gene_buckets.mega.merged.Magma.log10), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="-log10(Magma gene score)")
boxplot(rev(gene_buckets.mega.merged.h4), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Colocalisation h4")
boxplot(rev(gene_buckets_no_seed.mega.merged.Pascal.log10), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="-log10(Pascal gene score)")
boxplot(rev(gene_buckets_no_seed.mega.merged.Magma.log10), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="-log10(Magma gene score)")
boxplot(rev(gene_buckets_no_seed.mega.merged.h4), notch = TRUE, horizontal = TRUE, col = "black", medcol="darkorange", border = "gray", las=2, outline=FALSE, xlab="Colocalisation h4")
dev.off()







####################
# PATHWAY selection
####################


#read in pathways
# PASCAL
gwas.list.pathwayfiles.metabase <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/Metabase_GO_Maps_entrez.filt10_sum/", x, ".PathwaySet--Metabase_GO_Maps_entrez.filt10--sum.txt", sep="")
})
gwas.list.pathwayfiles.ppi <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/DSDSpectral.consensus_top50.1_ppi_anonym_v2_Entrez_sum/", x, ".PathwaySet--Metabase_GO_Maps_entrez.filt10--sum.txt", sep="")
})
gwas.list.pathwayfiles.coexp <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/DSDSpectral.consensus_top50.4_coexpr_anonym_v2_Entrez_sum/", x, ".PathwaySet--Metabase_GO_Maps_entrez.filt10--sum.txt", sep="")
})
gwas.list.pathwayfiles.reactome <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/c2.cp.reactome.v6.2.entrez_sum/", x, ".PathwaySet--Metabase_GO_Maps_entrez.filt10--sum.txt", sep="")
})
gwas.list.pathwayfiles.metabase_noghit <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/pascal/selected_109traits/", x, "/Metabase_GO_Maps_entrez.filt10_sum_noghit/", x, ".PathwaySet--Metabase_GO_Maps_entrez.filt10--sum.txt", sep="")
})
# MAGMA
gwas.list.pathwayfiles.magma <- sapply(gwas.list, function(x){
  paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/analysis/magma/magma_109traits_50kb/genesets_out/", x, "_Metabase_GO_Maps_entrez.filt10.out.sets.out", sep="")
})

# read in datasets
gwas.allFiles.paths.metabase <- lapply(gwas.list.pathwayfiles.metabase, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)})
gwas.allFiles.paths.ppi <- lapply(gwas.list.pathwayfiles.ppi, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)})
gwas.allFiles.paths.coexp <- lapply(ggwas.list.pathwayfiles.coexp, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)})
gwas.allFiles.paths.reactome <- lapply(gwas.list.pathwayfiles.reactome, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)})
gwas.allFiles.paths.metabase_noghit <- lapply(gwas.list.pathwayfiles.metabase_noghit, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)})
gwas.allFiles.paths.magma <- lapply(gwas.list.pathwayfiles.magma, function(x){read.delim(x, header=TRUE, row.names = NULL, sep="\t", stringsAsFactors = FALSE)})
# read in Deepali/Sarah analysis
gwas.allFiles.paths.dps.tbl <- read.delim("PathEnrich_Deepali/pathway_enrich.deepali.109_traits.txt")

# adjust p-vals
gwas.allFiles.paths.metabase <- lapply(gwas.allFiles.paths.metabase, function(x){x$empPvalue.adj <- p.adjust(x$empPvalue, method = "BH");x})
gwas.allFiles.paths.ppi <- lapply(gwas.allFiles.paths.ppi, function(x){x$empPvalue.adj <- p.adjust(x$empPvalue, method = "BH");x})
gwas.allFiles.paths.coexp <- lapply(gwas.allFiles.paths.coexp, function(x){x$empPvalue.adj <- p.adjust(x$empPvalue, method = "BH");x})
gwas.allFiles.paths.reactome <- lapply(gwas.allFiles.paths.reactome, function(x){x$empPvalue.adj <- p.adjust(x$empPvalue, method = "BH");x})
gwas.allFiles.paths.magma <- lapply(gwas.allFiles.paths.magma, function(x){x$COMP_P.adj <- p.adjust(x$COMP_P, method = "BH");x})

# separate dps per gwas and adjust pval 
gwas.list.trim <- gsub("gwas_results_", "", gwas.list)
gwas.allFiles.paths.dps <- lapply(gwas.list.trim, function(gwas){
  subset(gwas.allFiles.paths.dps.tbl, Trait == gwas)
})
names(gwas.allFiles.paths.dps) <- gwas.list
gwas.allFiles.paths.dps <- lapply(gwas.allFiles.paths.dps, function(x){x$Pval.adj <- p.adjust(x$Pval, method = "BH");x})



# select pathways with adj p-val < 0.05
gwas.allFiles.paths.filt10.adj05 <- lapply(gwas.allFiles.paths.filt10.adj, function(x){subset(x, empPvalue.adj < 0.05)})
gwas.allFiles.paths.filt10_100up.adj05 <- lapply(gwas.allFiles.paths.filt10.adj, function(x){subset(x, empPvalue.adj < 0.05)})
gwas.allFiles.paths.magma_50kb.adj05 <- lapply(gwas.allFiles.paths.magma_50kb.adj, function(x){subset(x, COMP_P.adj < 0.05)})

# comparing numbers of sig. pathways per gwas
p <- sapply(gwas.allFiles.paths.filt10.adj05, nrow)
m <- sapply(gwas.allFiles.paths.magma_50kb.adj05, nrow)
plot(p, m, xlab="Pascal", ylab="Magma", main="Number of enriched pathways/GWAS\n(adj.pval < 0.05)", cex=0.7, pch=19, col="black")
# Pascal default and 100up same numbers of pathways - continue with just default

# adjust gwas names to match megatable
names(gwas.allFiles.paths.filt10.adj05) <- gsub(pattern = "gwas_results_", replacement = "", x = names(gwas.allFiles.paths.filt10.adj05))
names(gwas.allFiles.paths.magma_50kb.adj05) <- gsub(pattern = "gwas_results_", replacement = "", x = names(gwas.allFiles.paths.magma_50kb.adj05))




# get success & gen.hit targets
megatable.updated.success <- subset(megatable.updated, success==1)
megatable.updated.failure <- subset(megatable.updated, failure==1)
megatable.updated.genhit <- subset(megatable.updated, high_confidence_genetic==1)

# for every gwas, extract sig. pathways with success without genetic hit
# get Megabase pathway lists
library(qusage)
metabase <- read.gmt("genesets/Metabase_GO_Maps_entrez.filt10.gmt")
metabase.magma <- gsub(" ", "_", names(metabase))


# filter gwas for ones with success
gwas.list.trim <- gsub(pattern = "gwas_results_", replacement = "", x = gwas.list)

gwas.sigpath.withsuccess.pascal <- lapply(gwas.list.trim, function(x){
  #print(x)
   m <- metabase[names(metabase) %in% gwas.allFiles.paths.filt10.adj05[[x]]$Name]
   if(length(m)!=0){
   s <- subset(megatable.updated.success, gwas==x)$entrezgene
   m[sapply(m, function(z){length(intersect(s,z))!=0})]
   } 
})
gwas.sigpath.withsuccess.magma <- lapply(gwas.list.trim, function(x){
  m <- metabase[metabase.magma %in% gwas.allFiles.paths.magma_50kb.adj05[[x]]$SET]
  if(length(m)!=0){
    s <- subset(megatable.updated.success, gwas==x)$entrezgene
    m[sapply(m, function(z){length(intersect(s,z))!=0})]
  } 
})
names(gwas.sigpath.withsuccess.pascal) <- gwas.list.trim
names(gwas.sigpath.withsuccess.magma) <- gwas.list.trim

# get successes per enriched pathway
gwas.sigpath.withsuccess.pascal.s <- lapply(gwas.list.trim, function(x){
  #print(x)
  m <- metabase[names(metabase) %in% gwas.allFiles.paths.filt10.adj05[[x]]$Name]
  if(length(m)!=0){
    s <- subset(megatable.updated.success, gwas==x)$entrezgene
    ms <- m[sapply(m, function(z){length(intersect(s,z))!=0})]
    mss <- lapply(ms, function(z){intersect(s,z)})
  } 
})
gwas.sigpath.withsuccess.magma.s <- lapply(gwas.list.trim, function(x){
  m <- metabase[metabase.magma %in% gwas.allFiles.paths.magma_50kb.adj05[[x]]$SET]
  if(length(m)!=0){
    s <- subset(megatable.updated.success, gwas==x)$entrezgene
    m[sapply(m, function(z){length(intersect(s,z))!=0})]
    mss <- lapply(ms, function(z){intersect(s,z)})
  } 
})
names(gwas.sigpath.withsuccess.pascal.s) <- gwas.list.trim
names(gwas.sigpath.withsuccess.magma.s) <- gwas.list.trim



# calculate percentage of sig pathways with/out gen hit
gwas.sigpath.withsuccess.pascal.ratio <- lapply(gwas.list.trim, function(x){
  m <- gwas.sigpath.withsuccess.pascal[[x]]
  if(length(m)!=0){
    s <- subset(megatable.updated.genhit, gwas==x)$entrezgene
    s.nohit <- length(m[sapply(m, function(z){length(intersect(s,z))==0})])
    s.hit <- length(m[sapply(m, function(z){length(intersect(s,z))!=0})])
    c(s.hit, s.nohit)
  } 
})
gwas.sigpath.withsuccess.magma.ratio <- lapply(gwas.list.trim, function(x){
  m <- gwas.sigpath.withsuccess.magma[[x]]
  if(length(m)!=0){
    s <- subset(megatable.updated.genhit, gwas==x)$entrezgene
    s.nohit <- length(m[sapply(m, function(z){length(intersect(s,z))==0})])
    s.hit <- length(m[sapply(m, function(z){length(intersect(s,z))!=0})])
    c(s.hit, s.nohit)
  } 
})
names(gwas.sigpath.withsuccess.pascal.ratio) <- gwas.list.trim
names(gwas.sigpath.withsuccess.magma.ratio) <- gwas.list.trim

pr <- unlist(sapply(gwas.sigpath.withsuccess.pascal.ratio, function(x){x[1]/(x[1]+x[2])}))
mr <- unlist(sapply(gwas.sigpath.withsuccess.magma.ratio, function(x){x[1]/(x[1]+x[2])}))

# plot percentage
boxplot(list(pr, mr), names=c("Pascal", "Magma"), col = "black", medcol="darkorange", border = "black", outline=TRUE, notch=TRUE, pch=19, ylab="Percentage", main="Enriched pathways\nharbouring a genetic hit")


############################################################################
# identify gwas with sig path, succes target but no genetic hit
gwas.sigpath.withsuccess.pascal.nohit <- lapply(gwas.list.trim, function(x){
  m <- gwas.sigpath.withsuccess.pascal[[x]]
  if(length(m)!=0){
    s <- subset(megatable.updated.genhit, gwas==x)$entrezgene
    m[sapply(m, function(z){length(intersect(s,z))==0})]
  } 
})
gwas.sigpath.withsuccess.magma.nohit <- lapply(gwas.list.trim, function(x){
  m <- gwas.sigpath.withsuccess.magma[[x]]
  if(length(m)!=0){
    s <- subset(megatable.updated.genhit, gwas==x)$entrezgene
    m[sapply(m, function(z){length(intersect(s,z))==0})]
  } 
})
names(gwas.sigpath.withsuccess.pascal.nohit) <- gwas.list.trim
names(gwas.sigpath.withsuccess.magma.nohit) <- gwas.list.trim


# identify the most significantly enriched pathway per gwas
pas.mod <- lapply(gwas.allFiles.paths.filt10.adj05, function(x){
  x$Name <- gsub(" ", "_", x$Name)
  x
})

gwas.sigpath.withsuccess.pascal.nohit.mostsig <- lapply(gwas.list.trim, function(x){
  m <- gwas.sigpath.withsuccess.pascal.nohit[[x]]
  if(length(m)!=0){
    ms <- subset(gwas.allFiles.paths.filt10.adj05[[x]], Name %in% names(m))
    subset(ms, empPvalue==min(ms$empPvalue))
  } 
})
names(gwas.sigpath.withsuccess.pascal.nohit.mostsig) <- gwas.list.trim

# focusing on 
# "IAberrant production of IL-2 and IL-17 in SLE T cells" in "GSK500KV3lin_M2W_psoriasis"
# 30 Aberrant production of IL-2 and IL-17 in SLE T cells 5.088258e-07     4e-07  3.254118e-05
# get all genes
# forom Metabaser
#SLEttargets <- subset(get.map.drugtargets(7038), disease_name == "Psoriasis")
#SLEt <- get.map.genes(7038)
# get genescores
pascal.psoriasis <- subset(megatable.updated, gwas=="GSK500KV3lin_M2W_psoriasis")
pascal.psoriasis.slet <- subset(pascal.psoriasis, entrezgene %in% SLEt$gene)
pascal.psoriasis.slet <- pascal.psoriasis.slet[order(pascal.psoriasis.slet$entrezgene, decreasing=TRUE),]
pascal.psoriasis.slet <- pascal.psoriasis.slet[!duplicated(pascal.psoriasis.slet$entrezgene),]
pascal.psoriasis.slet <- pascal.psoriasis.slet[order(as.numeric(pascal.psoriasis.slet$Pascal_default)),]
##
# exclude HLS genes
hla <- read.table("scripts/PASCAL/PASCAL/resources/annotation/hla/hlaGenesEntrezIds.txt", header = TRUE)$gene_id
pascal.psoriasis.slet <- subset(pascal.psoriasis.slet, !(entrezgene %in% hla))
plot(rev(-log10(as.numeric(pascal.psoriasis.slet$Pascal_default))), 1:length(pascal.psoriasis.slet$Pascal_default), yaxt="n", xlab="-log10(Pascal gene score)", ylab="", col="black", pch=19)
axis(2,at=1:length(pascal.psoriasis.slet$Pascal_default),labels=pascal.psoriasis.slet$hgnc_symbol, las=2)


# focusing on 
# "Immune response_OX40L_OX40 signaling pathway" in "GSK500KV3lin_M2W_asthma"
# 30 Aberrant production of IL-2 and IL-17 in SLE T cells 5.088258e-07     4e-07  3.254118e-05
# get all genes
# forom Metabaser
#SLEttargets <- subset(get.map.drugtargets(7038), disease_name == "Psoriasis")
#SLEt <- get.map.genes(7038)
# get genescores
pascal.psoriasis <- subset(megatable.updated, gwas=="GSK500KV3lin_M2W_psoriasis")
pascal.psoriasis.slet <- subset(pascal.psoriasis, entrezgene %in% SLEt$gene)
pascal.psoriasis.slet <- pascal.psoriasis.slet[order(pascal.psoriasis.slet$entrezgene, decreasing=TRUE),]
pascal.psoriasis.slet <- pascal.psoriasis.slet[!duplicated(pascal.psoriasis.slet$entrezgene),]
pascal.psoriasis.slet <- pascal.psoriasis.slet[order(as.numeric(pascal.psoriasis.slet$Pascal_default)),]
##
# exclude HLS genes
hla <- read.table("scripts/PASCAL/PASCAL/resources/annotation/hla/hlaGenesEntrezIds.txt", header = TRUE)$gene_id
pascal.psoriasis.slet <- subset(pascal.psoriasis.slet, !(entrezgene %in% hla))
plot(rev(-log10(as.numeric(pascal.psoriasis.slet$Pascal_default))), 1:length(pascal.psoriasis.slet$Pascal_default), yaxt="n", xlab="-log10(Pascal gene score)", ylab="", col="black", pch=19)
axis(2,at=1:length(pascal.psoriasis.slet$Pascal_default),labels=pascal.psoriasis.slet$hgnc_symbol, las=2)



pn <- names(gwas.sigpath.withsuccess.pascal.nohit[sapply(gwas.sigpath.withsuccess.pascal.nohit, length) > 0])
mn <- names(gwas.sigpath.withsuccess.pascal.nohit[sapply(gwas.sigpath.withsuccess.magma.nohit, length) > 0])
intersect(pn, mn)




############################################
save.image(file="analysis/merging_genescores_250619.RData")


