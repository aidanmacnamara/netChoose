###############################################################
## Getting list of UKB SNPs at p< 5.0e10-8 level and finding ##
## thier corresponding gene with piccolo evidence of H4>=0.9 ##
## Ali Amin 22nd Jan 2019                                    ##
###############################################################

piccolo_90<-read.csv("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/GWAS_SNPS_V2G_90h4from_piccolo.csv")
disease_name<-read.table("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/disease_names_m2w_filterd.csv",header=FALSE,sep=",")
colnames(disease_name)<-c("a","phenotype1")
disease_name$common<-disease_name$phenotype1%in%piccolo_90$trait
table(disease_name$common)
## 34 phenotypes can be merged (name of disease is similar in both UKB and piccolo datasets, IBD is among them.).
## 323 phenotypes can't be merged (name of disease is not similar in both UKB and piccolo datasets.This needs more investigation.).
phenotype_passed<-disease_name[disease_name$common==TRUE,]
phenotype_notpased<-disease_name[disease_name$common==FALSE,]

##################################################################################
## Merging the results of UKB (p<5.0e-08) with Piccolo (H4>0.9) for 34          ##
## common phenotypes. Files generated with "..._v2g_piccolo_h90.txt" extention  ##
## include genes with high confidence for V2G using piccolo H4>=0.90            ##
##################################################################################

all_phenotypes1<-NULL

for (i in phenotype_passed[,2]){
  file1<-paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/gwas_results_GSK500KV3lin_M2W_",i,"_pvalfilter_phenotype.txt",sep="")
  file2<-paste("/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/gwas_results_GSK500KV3lin_M2W_",i,"_v2g_piccolo_h90.txt",sep="") 
  data1<-read.table(file=file1,header=TRUE,sep="\t")
  names(data1)<-c("rsID","pval_ukb","trait","MeSH_ID")
  data2<-merge(data1,piccolo_90,by=c("rsID","trait"),all.x=TRUE)
  data2$hgnc_idx<-as.character(data2$hgnc_idx)
  data2$rsID<-as.character(data2$rsID)
  data<-subset(data2,is.na(data2$hgnc_idx)==FALSE)
  write.table(data,file=file2,sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
  all_phenotypes1<-rbind(all_phenotypes1,data)
}

dim(all_phenotypes1)
write.table(all_phenotypes1,file="/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/gwas_results_GSK500KV3lin_M2W_all_phenotypes1_v2g_piccolo_h90.txt",sep="\t",row.names = FALSE,col.names = TRUE,quote = FALSE)

## to make a list for all disease and getting corresponding genes:

data<-all_phenotypes1
phenotype_UKBM2W_genes<-NULL

for (i in unique(data$trait)){
  print(i)
  a<-as.data.frame(table(as.character(unique(data[data$trait==i,]$hgnc_idx))))
  b<-list(as.character(a$Var1))
  phenotype_UKBM2W_genes<-append(phenotype_UKBM2W_genes,b)
}

name<-as.character(unique(data$trait))
names(phenotype_UKBM2W_genes)<-name
save(phenotype_UKBM2W_genes,file="/GWD/appbase/projects/rd-scratchdata-gdc/TSCI/CB/network_analysis/gwas/ukb500kv3_m2w_filter/phenotype_UKBM2W_genes_v2g_piccolo_90h4.Rdata")

q("no")

