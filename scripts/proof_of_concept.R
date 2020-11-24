
# LOAD --------------------------------------------------------------------

require(tidyverse)
require(readxl)
require(biomaRt)
require(eulerr)
require(jsonlite)
require(fgsea)
require(UpSetR)


# GWAS --------------------------------------------------------------------

# gwas catalog for proof of concept
# start off with ibd

ibd_gwas = read_csv("gwas/download_IBD.csv")

# apply filters
p_value = as.numeric(str_replace(ibd_gwas$`P-value`, " x 10", "e"))
ibd_gwas$`P-value` = p_value
ibd_gwas_filt = filter(ibd_gwas, `Reported trait`=="Inflammatory bowel disease", `P-value`<=5e-8)
tail(sort(table(ibd_gwas_filt$`Mapped gene`)))

# define hits as above for the moment, mapping to entrez
hits = sapply(ibd_gwas_filt$`Mapped gene`, function(x) (str_extract_all(x, "[[:alnum:]-]+")))
hits_df = data.frame(cbind(p_value=as.numeric(ibd_gwas_filt$`P-value`), hits, location=ibd_gwas_filt$Location))
hits_df = hits_df %>% separate_rows(hits, sep=",")
hits_df$hits = str_replace_all(hits_df$hits, "(\\))|(c\\()|\"", "")
hits_df$p_value = as.numeric(unlist(hits_df$p_value))
hits_df$location = as.character(unlist(hits_df$location))

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mapping_1 = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id"), mart=mart)

hits_df$entrez = mapping_1$entrezgene[match(hits_df$hits, mapping_1$hgnc_symbol)]

hits_df_entrez = hits_df[!is.na(hits_df$entrez),]
hits_df_entrez$score = -log(hits_df_entrez$p_value, base=10)
hits_df_entrez = hits_df_entrez %>% group_by(hits) %>% slice(which.max(score)) %>% arrange(desc(score))


# HOTNET2 PREP ------------------------------------------------------------

gene_ix = read_tsv("/Volumes/am673712/links/network_analysis/interactions/metabase_plus_ppi/metabase_plus_ppi_gene_index.txt", col_names=FALSE)
heat_scores = data.frame(gene=gene_ix$X2, score=NA)
heat_scores$score = hits_df_entrez$score[match(heat_scores$gene, hits_df_entrez$hits)]
heat_scores$score[is.na(heat_scores$score)] = 0
heat_scores = arrange(heat_scores, desc(score))
write_tsv(heat_scores, "/Volumes/am673712/links/network_analysis/hotnet/test_scores.tsv", col_names=FALSE)

# run on slurm
# sbatch -n 1 -c 12 --mem=50G -t 360 -o slurm_log.out run_test.sh


# TIPS DATA ---------------------------------------------------------------

tips_clin = read_excel("tips/Pprojects_1Nov2018_TIPSinterpretation.xlsx")
tips_clin %>% filter(`Clinical Label_PP`=="Succeeded", DiseaseType=="Non-Neoplasm") %>% group_by(`MeSH ID`) %>% summarise(N=n()) %>% arrange(desc(N))

# data for ibd
filter(tips_clin, `MeSH ID`=="D015212") %>% group_by(`Clinical Label_PP`) %>% summarise(N=n())
ibd_ts = filter(tips_clin, `MeSH ID` %in% c("D015212","D003424","D003093")) # get children as well - need to code this
ibd_ts$gene = mapping_1$hgnc_symbol[match(ibd_ts$`EntrezGene ID`, mapping_1$entrezgene)]
# ibd_ts = ibd_ts %>% filter(`Clinical Label_PP`=="Succeeded") %>% dplyr::select("EntrezGene ID") 
ibd_ts = ibd_ts %>% dplyr::select("EntrezGene ID", gene) 
names(ibd_ts)[1] = "entrez"

ibd_ts = distinct(ibd_ts)


# OVERLAPS ----------------------------------------------------------------

# gwas results (closest gene)

plot(euler(
  list(
    gwas=hits_df_entrez$entrez,
    drug_targets=ibd_ts$entrez
  )
), quantities=TRUE)

# hotnet results

hotnet_res = read_tsv("hotnet/hotnet/results/consensus/subnetworks.tsv", col_names=FALSE)
hotnet_res = unlist(str_split(hotnet_res$X1[3:length(hotnet_res$X1)], " "))
hotnet_res = data.frame(gene=hotnet_res, entrez=mapping_1$entrezgene[match(hotnet_res, mapping_1$hgnc_symbol)])

# all hotnet results - all

dirs = list.dirs("hotnet/hotnet/results/full_network-test_scores/")[-1]

res_all = vector("list", length(dirs))
for(i in 1:length(dirs)) {
  res = read_json(paste0(dirs[i], "/results.json"))
  res = unique(unlist(res$components))
  res_all[[i]] = data.frame(gene=res, entrez=mapping_1$entrezgene[match(res, mapping_1$hgnc_symbol)])
}

# hierarchical hotnet

hier_res = read_csv("hotnet/hierarchical_hotnet/results/clusters_metabase_plus_ppi_scores_1.tsv", skip=7, col_names=FALSE)
hier_res = unlist(str_split(hier_res$X1, "\t"))
hier_res = data.frame(gene=hier_res, entrez=mapping_1$entrezgene[match(hier_res, mapping_1$hgnc_symbol)])

# complex data

net_complex = read_tsv("interactions/complex_portal/net_complex.tsv")

c_ix = apply(net_complex[,1:2], 1, function(x) any(x %in% hits_df_entrez$entrez))
complex_res = data.frame(entrez = unique(unlist(net_complex[c_ix,1:2])))

# ligand receptor data

lr_dat = read.table("interactions/ligand_receptor_db/LigandReceptor.tsv")
lr_ix = apply(dplyr::select(lr_dat, Ligand.Entrez, Receptor.Entrez), 1, function(x) any(x %in% hits_df_entrez$entrez))
table(lr_ix)
lr_res = data.frame(entrez = unique(unlist(lr_dat[lr_ix,which(names(lr_dat) %in% c("Ligand.Entrez","Receptor.Entrez"))])))

# network neighbours

source("/GWD/appbase/projects/RD-TSci-Software/CB/packages/R/MetaBase.R")

trust = show.interaction.trusts()$trust_name[1:4]
mechanism = show.interaction.mechanisms()
mechanism = mechanism$mechanism_name[mechanism$direct]

first_n = rbind(
  find_neighbours(hits_df_entrez$entrez, direction="downstream", trust=trust, mechanism=mechanism),
  find_neighbours(hits_df_entrez$entrez, direction="upstream", trust=trust, mechanism=mechanism)
)

first_n$gene_in = as.numeric(first_n$gene_in)
first_n$gene_out = as.numeric(first_n$gene_out)

x = first_n %>% dplyr::filter(direction=="downstream") %>% group_by(gene_in) %>% summarise(N=n(), gene=genesymbol_in[1]) %>% arrange(desc(N))
table(x$gene_in %in% hits_df_entrez$entrez)
x = first_n %>% dplyr::filter(direction=="upstream") %>% group_by(gene_out) %>% summarise(N=n(), gene=genesymbol_out[1]) %>% arrange(desc(N))
table(x$gene_out %in% hits_df_entrez$entrez)

first_n_res = data.frame(entrez=as.numeric(unique(unlist(first_n[,3:4]))))
table(hits_df_entrez$entrez %in% first_n_res$entrez)

# pathways

metabase = gmtPathways("genesets/Metabase_GO_Maps_entrez.filt10.gmt")
p_ix = unlist(lapply(metabase, function(x) any(hits_df_entrez$entrez %in% x)))
table(p_ix)

plot(euler(
  list(
    gwas=hits_df_entrez$entrez,
    drug_targets=ibd_ts$entrez,
    hotnet=hotnet_res$entrez,
    #    hotnet_delta=res_all[[4]]$entrez
    complex=complex_res$entrez,
    ligand_receptor=lr_res$entrez
    # neighbours=first_n_res$entrez
  )
), quantities=TRUE)

all_data = list(gwas=hits_df_entrez$entrez,
                drug_targets=ibd_ts$entrez,
                hotnet=hotnet_res$entrez,
                hierarchical=hier_res$entrez,
                # hotnet_delta=res_all[[4]]$entrez,
                complex=complex_res$entrez,
                ligand_receptor=lr_res$entrez,
                neighbours=first_n_res$entrez
)

all_data_tbl = as.data.frame.matrix((table(stack(all_data))))
all_data_tbl = cbind(rownames(all_data_tbl), all_data_tbl)
rownames(all_data_tbl) = NULL
upset(all_data_tbl, order.by="freq", nsets=20)

