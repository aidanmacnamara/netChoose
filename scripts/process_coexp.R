
#merge genescores across GWAS
gwas.allFiles.coexp.merged <- lapply(gwas.list, function(gwas){
  y <- lapply(gwas.allFiles.coexp, function(x){
    x[[gwas]]
  })
  Reduce(function(...) merge(..., by="Name"), y)
})
names(gwas.allFiles.coexp.merged) <- gwas.list.trim
# tidy names
gwas.allFiles.coexp.merged <- lapply(gwas.allFiles.coexp.merged, function(x){
  setNames(x, nm = c("Pathway", "Pascal.coexp.chiPvalue", "Pascal.coexp.Pvalue","Pascal.coexp.Pvalue.adj",
                     "Magma.coexp.NGENES", "Magma.coexp.SELF_P", "Magma.coexp.Pvalue","Magma.coexp.Pvalue.adj"))
})
gwas.allFiles.coexp.merged.tidyPval <- lapply(gwas.allFiles.coexp.merged, function(x){x[,c("Pathway", "Pascal.coexp.Pvalue", "Magma.coexp.Pvalue")]})
gwas.allFiles.coexp.merged.tidyPval.adj <- lapply(gwas.allFiles.coexp.merged, function(x){x[,c("Pathway", "Pascal.coexp.Pvalue.adj", "Magma.coexp.Pvalue.adj")]})


# pval.adj-to-pval across all GWAS
gwas.allFiles.coexp.merged.tidyPval.adj.sig <- lapply(gwas.allFiles.coexp.merged.tidyPval.adj, function(x){
  subset(x, Pascal.coexp.Pvalue.adj < 0.05 | Magma.coexp.Pvalue.adj < 0.05)
})
gwas.allFiles.coexp.merged.tidyPval.adj.sig.n <- lapply(gwas.allFiles.coexp.merged.tidyPval.adj, function(x){
  a <- nrow(subset(x, Pascal.coexp.Pvalue.adj < 0.05))
  b <- nrow(subset(x, Magma.coexp.Pvalue.adj < 0.05))
  c(a,b)
})

# first, pull all genes from pathways for all sig pathways
gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths <- lapply(gwas.allFiles.coexp.merged.tidyPval.adj.sig, function(x){
  p <- subset(x, Pascal.coexp.Pvalue.adj < 0.05)$Pathway
  m <- subset(x, Magma.coexp.Pvalue.adj < 0.05)$Pathway
  list(p,m)
})

gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes <- lapply(gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths, function(x){
  lapply(x, function(y){
    coexp.magma[y]
  })
})


# separate pathways/genes with genetic hit
gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes.wghit <- lapply(1:length(gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes), function(i){
  x <- gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes[[i]]
  xname <- names(gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes)[[i]]
  lapply(x, function(y){
    if(!is_empty(y)){y[sapply(y, function(z){any(z %in% gwas.list.trim.hits[[xname]])})]} else {return(NA)}
  })
})
gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes.woghit <- lapply(1:length(gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes), function(i){
  x <- gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes[[i]]
  xname <- names(gwas.allFiles.coexp.merged.tidyPval.adj.sig.paths.genes)[[i]]
  lapply(x, function(y){
    if(!is_empty(y)){y[sapply(y, function(z){!any(z %in% gwas.list.trim.hits[[xname]])})]} else {return(NA)}
  })
})
