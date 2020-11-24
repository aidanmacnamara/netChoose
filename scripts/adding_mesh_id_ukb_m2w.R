
#############################################################
## The MeSH.ID of 1470 phenotypes was provided by Nikolina ##
## UKB results was provided by Nikolina for 357 M2W after  ##
## I filtered UKB-M2W results by P< 5.0e-08. This script   ##
## adds the name of phenotypes and MeSH IDs to the results ##
## Ali Amin 21st Jan 2019                                  ##  
#############################################################

mesh<-read.table("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/mesh_to_analysis.csv",header=TRUE,sep=",")
mesh$dup<-duplicated(mesh$phenotype1)
mesh<-mesh[mesh$dup==FALSE,]
mesh<-mesh[,c(1:2)]

## I extract 382 disease names from UKB file names 
disease_name<-read.table("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/disease_names_m2w_filterd.csv",header=FALSE,sep=",")
dim(disease_name)
colnames(disease_name)<-c("a","phenotype1")
table(disease_name[,2]%in%mesh[,1])

disease_name<-merge(disease_name,mesh,by="phenotype1",all.x=TRUE)
disease_name<-disease_name[,c(1,3)]

## Adding Phenotype and MeSH_IDs to the UKB results
## The outcome files have an extention "...pvalfilter_phenotype.txt"
for (i in disease_name[,1]){
  print(i)
  Phenotype<-i
  MeSH_ID<-disease_name[disease_name$phenotype1==i,]$MeSH.ID
  file1<-paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/gwas_results_GSK500KV3lin_M2W_",i,"_pvalfilter.txt",sep="")
  file2<-paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/gwas_results_GSK500KV3lin_M2W_",i,"_pvalfilter_phenotype.txt",sep="")
  print(file1)
  print(file2)
  a<-read.table(file=file1,sep="\t",header=FALSE)
  a<-a[,c(1:2)]
  a[,3]<-Phenotype
  a[,4]<-MeSH_ID
  a<-a[,c(1:4)]
  colnames(a)<-c("SNP","pval","phenotype","MeSH_ID")
  write.table(a,file=file2,row.names = FALSE,col.names = TRUE,sep="\t",quote = FALSE)
}

q("no")
