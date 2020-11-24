###############################################################
### Selection of 100 Random Genes  (Ali Amin 15th Jan 2019) ###
###############################################################

library(biomaRt)
# to connect to the Ensembl live gene mart human dataset (GRCh38):
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# to get all human genes names
genes<-NULL
for (i in c(1:22,"X","Y","MT")){
  chr_genes <- getBM(attributes=c('ensembl_gene_id',
                                  'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters =
                       'chromosome_name', values =i, mart = ensembl)
  print (i)
  print(dim(chr_genes))
  genes<-rbind(genes,chr_genes)
}

genes$dup<-duplicated(genes$hgnc_symbol)
genes<-genes[genes$dup==FALSE,]
dim(genes)
#[1] 37550 7

## random selection of 100 Genes
random_genes<-sample(genes$hgnc_symbol, 100, replace = FALSE, prob = NULL)
length(random_genes)
#[1] 100
genes$include<-genes$hgnc_symbol%in%random_genes
random_genes<-genes[genes$include==TRUE,]
write.csv(random_genes,file="random_100genes.csv")
