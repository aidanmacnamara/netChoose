
#merge genescores across GWAS
gwas.allFiles.ppi.merged <- lapply(gwas.list, function(gwas){
  y <- lapply(gwas.allFiles.ppi, function(x){
    x[[gwas]]
  })
  Reduce(function(...) merge(..., by="Name"), y)
})
names(gwas.allFiles.ppi.merged) <- gwas.list.trim
# tidy names
gwas.allFiles.ppi.merged <- lapply(gwas.allFiles.ppi.merged, function(x){
  setNames(x, nm = c("Pathway", "Pascal.ppi.chiPvalue", "Pascal.ppi.Pvalue","Pascal.ppi.Pvalue.adj",
                     "Magma.ppi.NGENES", "Magma.ppi.SELF_P", "Magma.ppi.Pvalue","Magma.ppi.Pvalue.adj"))
})
gwas.allFiles.ppi.merged.tidyPval <- lapply(gwas.allFiles.ppi.merged, function(x){x[,c("Pathway", "Pascal.ppi.Pvalue", "Magma.ppi.Pvalue")]})
gwas.allFiles.ppi.merged.tidyPval.adj <- lapply(gwas.allFiles.ppi.merged, function(x){x[,c("Pathway", "Pascal.ppi.Pvalue.adj", "Magma.ppi.Pvalue.adj")]})


# pval.adj-to-pval across all GWAS
gwas.allFiles.ppi.merged.tidyPval.adj.sig <- lapply(gwas.allFiles.ppi.merged.tidyPval.adj, function(x){
  subset(x, Pascal.ppi.Pvalue.adj < 0.05 | Magma.ppi.Pvalue.adj < 0.05)
})
gwas.allFiles.ppi.merged.tidyPval.adj.sig.n <- lapply(gwas.allFiles.ppi.merged.tidyPval.adj, function(x){
  a <- nrow(subset(x, Pascal.ppi.Pvalue.adj < 0.05))
  b <- nrow(subset(x, Magma.ppi.Pvalue.adj < 0.05))
  c(a,b)
})

# first, pull all genes from pathways for all sig pathways
gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths <- lapply(gwas.allFiles.ppi.merged.tidyPval.adj.sig, function(x){
  p <- subset(x, Pascal.ppi.Pvalue.adj < 0.05)$Pathway
  m <- subset(x, Magma.ppi.Pvalue.adj < 0.05)$Pathway
  list(p,m)
})

gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes <- lapply(gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths, function(x){
  lapply(x, function(y){
    ppi.magma[y]
  })
})


# separate pathways/genes with genetic hit
gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes.wghit <- lapply(1:length(gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes), function(i){
  x <- gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes[[i]]
  xname <- names(gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes)[[i]]
  lapply(x, function(y){
    if(!is_empty(y)){y[sapply(y, function(z){any(z %in% gwas.list.trim.hits[[xname]])})]} else {return(NA)}
  })
})
gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes.woghit <- lapply(1:length(gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes), function(i){
  x <- gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes[[i]]
  xname <- names(gwas.allFiles.ppi.merged.tidyPval.adj.sig.paths.genes)[[i]]
  lapply(x, function(y){
    if(!is_empty(y)){y[sapply(y, function(z){!any(z %in% gwas.list.trim.hits[[xname]])})]} else {return(NA)}
  })
})
