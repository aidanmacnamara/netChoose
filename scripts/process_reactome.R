
#merge genescores across GWAS
gwas.allFiles.reactome.merged <- lapply(gwas.list, function(gwas){
  y <- lapply(gwas.allFiles.reactome, function(x){
    x[[gwas]]
  })
  Reduce(function(...) merge(..., by="Name"), y)
})
names(gwas.allFiles.reactome.merged) <- gwas.list.trim
# tidy names
gwas.allFiles.reactome.merged <- lapply(gwas.allFiles.reactome.merged, function(x){
  setNames(x, nm = c("Pathway", "Pascal.reactome.chiPvalue", "Pascal.reactome.Pvalue","Pascal.reactome.Pvalue.adj",
                     "Magma.reactome.NGENES", "Magma.reactome.SELF_P", "Magma.reactome.Pvalue","Magma.reactome.Pvalue.adj"))
})
gwas.allFiles.reactome.merged.tidyPval <- lapply(gwas.allFiles.reactome.merged, function(x){x[,c("Pathway", "Pascal.reactome.Pvalue", "Magma.reactome.Pvalue")]})
gwas.allFiles.reactome.merged.tidyPval.adj <- lapply(gwas.allFiles.reactome.merged, function(x){x[,c("Pathway", "Pascal.reactome.Pvalue.adj", "Magma.reactome.Pvalue.adj")]})


# pval.adj-to-pval across all GWAS
gwas.allFiles.reactome.merged.tidyPval.adj.sig <- lapply(gwas.allFiles.reactome.merged.tidyPval.adj, function(x){
  subset(x, Pascal.reactome.Pvalue.adj < 0.05 | Magma.reactome.Pvalue.adj < 0.05)
})
gwas.allFiles.reactome.merged.tidyPval.adj.sig.n <- lapply(gwas.allFiles.reactome.merged.tidyPval.adj, function(x){
  a <- nrow(subset(x, Pascal.reactome.Pvalue.adj < 0.05))
  b <- nrow(subset(x, Magma.reactome.Pvalue.adj < 0.05))
  c(a,b)
})

# first, pull all genes from pathways for all sig pathways
gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths <- lapply(gwas.allFiles.reactome.merged.tidyPval.adj.sig, function(x){
  p <- subset(x, Pascal.reactome.Pvalue.adj < 0.05)$Pathway
  m <- subset(x, Magma.reactome.Pvalue.adj < 0.05)$Pathway
  list(p,m)
})

gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes <- lapply(gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths, function(x){
  lapply(x, function(y){
    reactome.magma[y]
  })
})


# separate pathways/genes with genetic hit
gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes.wghit <- lapply(1:length(gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes), function(i){
  x <- gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes[[i]]
  xname <- names(gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes)[[i]]
  lapply(x, function(y){
    if(!is_empty(y)){y[sapply(y, function(z){any(z %in% gwas.list.trim.hits[[xname]])})]} else {return(NA)}
  })
})
gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes.woghit <- lapply(1:length(gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes), function(i){
  x <- gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes[[i]]
  xname <- names(gwas.allFiles.reactome.merged.tidyPval.adj.sig.paths.genes)[[i]]
  lapply(x, function(y){
    if(!is_empty(y)){y[sapply(y, function(z){!any(z %in% gwas.list.trim.hits[[xname]])})]} else {return(NA)}
  })
})
