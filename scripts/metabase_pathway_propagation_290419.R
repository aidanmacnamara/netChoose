library(dplyr)
library(tximport)
library(readr)
library(ggplot2)
library(RColorBrewer)

library(metabaser)
#metabase.connect(host = "xxxx", jdbc.url = "xxxxxx", uid = "xxxx", pwd = "XXXXX", driver = "jdbc")



## load interactions on all pathway maps
pwm_network <- distinct(get.map.interactions(get.maps()))
gids <- unique(c(pwm_network$id1, pwm_network$id2))

# get names of all genes in maps, human
gg <- get.genes(filter=list(has.maps = TRUE, species="Homo sapiens"))



# function for extracting 1st and 2nd interactors
get.interactors <- function(gene){
  ## Interactions use network objects, not Entrez Genes, so we'll convert this to a network object
  gene.nwobj <- get.gene.nwobjs(gene)
  
  ## Find first degree interactions downstream of gene
  ## To do this, I'll search for the network object gene as "id1" in the pathway map interaction network I created
  gene.out.inter1 <- subset(pwm_network, id1 == gene.nwobj$name)
  ## inverse if I want all incoming interactions
  gene.in.inter1 <- subset(pwm_network, id2 == gene.nwobj$name)
  ## get gene IDS for all in/out interactors
  in1 <- get.genes(filter=list(genesymbol=gene.in.inter1$id1))$gene
  out1 <- get.genes(filter=list(genesymbol=gene.out.inter1$id2))$gene
  genes.1st <- unique(c(in1, out1))
  
  ## take advantage of these lists to find the 2nd neighbor interactions downstream
  gene.out.inter2 <- subset(pwm_network, id1 %in% gene.out.inter1$id2)
  ## do similar for 2nd neighbor upstream interactions
  gene.in.inter2 <- subset(pwm_network, id2 %in% gene.in.inter1$id1)
  ## get gene IDS for all in/out interactors
  in2 <- get.genes(filter=list(genesymbol=gene.in.inter2$id1))$gene
  out2 <- get.genes(filter=list(genesymbol=gene.out.inter2$id2))$gene
  genes.2nd <- unique(c(in2, out2))
  
  return(list(genes.1st, genes.2nd))
}
  

# run function across all map genes
gg.interactors <- lapply(gg$gene, get.interactors)
gg.interactors <- lapply(gg.interactors, function(x){setNames(x, c("1st_interactor","2nd_interactor"))})
names(gg.interactors) <- gg$gene




# save object
save(gg.interactors, file = "Metabase_interactors.Rdata")




