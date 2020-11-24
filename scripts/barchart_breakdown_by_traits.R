library(tidyverse)

info1<-read.table("mega_tbl.txt",header=TRUE,sep="\t")
gwas_info<-read.table("gwas_list.txt",header=TRUE,sep="\t")
info1$gwas_flag<-info1$gwas%in%gwas_info$gwas

df<-info1 %>% filter(gwas_flag==TRUE) %>% group_by(disease) %>% summarise(No.GWAS=n_distinct(gwas))
df<-df[order(df$No.GWAS,decreasing = TRUE),]
df$No.Disease<-c(1:nrow(df))
df$Disease<-paste(df$No.GWAS,df$disease,sep = "-")

# Barplot (PDF and PNG format)
pdf(file="Barplot_GWAS_per_disease.pdf",paper = "a4r")
ggplot(df, aes(x=No.Disease, y=No.GWAS, color=Disease)) +  geom_bar(stat="identity", fill="white") + ggtitle("Number Of GWAS Per Disease")
dev.off()

png(file="Barplot_GWAS_per_disease.png",width = 640, height = 320)
ggplot(df, aes(x=No.Disease, y=No.GWAS, color=Disease)) +  geom_bar(stat="identity", fill="white") + ggtitle("Number Of GWAS Per Disease")
dev.off()

q("no")


