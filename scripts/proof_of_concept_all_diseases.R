
# LOAD --------------------------------------------------------------------

require(tidyverse)
require(biomaRt)
require(eulerr)
require(jsonlite)
require(readxl)


# MAPPING -----------------------------------------------------------------

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mapping = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id"), mart=mart)


# GWAS --------------------------------------------------------------------

# traits: significant hits from piccolo/co-localisation
traits = read_excel("gwas/Selected_77traits_info.xlsx", sheet=3)

# tips
tips_clin = read_excel("tips/Pprojects_1Nov2018_TIPSinterpretation.xlsx")
table(traits$MESH_ID %in% tips_clin$`MeSH ID`) # all co-localisation traits have clinical data in tips

mesh = unique(traits$MESH_ID)
mesh_df = data.frame(mesh=mesh, gwas=NA, clinical=NA)
for(i in 1:length(mesh)) {
  mesh_df[i,] = c(mesh[i], length(traits$MESH_ID[traits$MESH_ID==mesh[i]]), filter(tips_clin, `MeSH ID`==mesh[i], `Clinical Label_PP`=="Succeeded") %>% summarise(n()))
}
ggplot(mesh_df, aes(gwas, clinical)) + geom_point() + theme_thesis() + xlab("GWAS") + ylab("Success")

# add entrez id
traits$entrez = mapping$entrezgene[match(traits$ensembl_id, mapping$ensembl_gene_id)] # 432 missing values
traits = traits[!is.na(traits$entrez),]
traits_summ = traits %>% group_by(MESH_ID) %>% summarise(N=n(), trait=trait[1]) %>% arrange(desc(N))


# HOTNET2 PREP ------------------------------------------------------------

gene_ix = read_tsv("/Volumes/am673712/links/network_analysis/interactions/metabase_plus_ppi/metabase_plus_ppi_gene_index.txt", col_names=FALSE)

for(i in 1:dim(traits_summ)[1]) {
  
  heat_scores = data.frame(gene=gene_ix$X2, score=NA)
  traits_filt = filter(traits, MESH_ID==traits_summ$MESH_ID[i])
  heat_scores$score = traits_filt$log_gene_score_pph[match(heat_scores$gene, traits_filt$hgncid)]
  heat_scores$score[is.na(heat_scores$score)] = 0
  heat_scores = arrange(heat_scores, desc(score))
  write_tsv(heat_scores, paste0("/Volumes/am673712/links/network_analysis/hotnet/hotnet/test_all_traits/input_all_traits/input_", i, ".tsv"), col_names=FALSE)
}

# run on slurm
# sbatch -n 1 -c 12 --mem=50G -t 360 -o slurm_log.out run_test.sh


# RESULTS DATA FRAME ------------------------------------------------------

res = tbl_df(data.frame(gene=rep(gene_ix$X2, dim(traits_summ)[1]), disease=rep(traits_summ$MESH_ID, each=dim(gene_ix)[1]), gwas=0, success=0))
res$entrez = mapping$entrezgene[match(res$gene, mapping$hgnc_symbol)]
res$index = 1:dim(res)[1]

for(i in 1:length(mesh)) {
  print(i)
  successes = filter(tips_clin, `Clinical Label_PP`=="Succeeded", `MeSH ID`==mesh[i]) %>% dplyr::select(`EntrezGene ID`) %>% unlist()
  gwas = filter(traits, MESH_ID==mesh[i]) %>% dplyr::select(entrez) %>% unlist()
  res$success[(res$disease==mesh[i] & res$entrez %in% successes)] = 1
  res$gwas[(res$disease==mesh[i] & res$entrez %in% gwas)] = 1
}


# OVERLAPS ----------------------------------------------------------------

# gwas results (closest gene)

plot(euler(
  list(
    gwas=res$index[res$gwas==1],
    successes=res$index[res$success==1]
  )
), quantities=TRUE)

# hotnet results

res$hotnet = 0
for(i in 1:dim(traits_summ)) {
  print(i)
  hotnet_res = read_tsv(paste0("/Volumes/am673712/links/network_analysis/hotnet/hotnet/test_all_traits/results/output_", i, "/consensus/subnetworks.tsv"), col_names=FALSE)
  hotnet_res = unlist(str_split(hotnet_res$X1[3:length(hotnet_res$X1)], " "))
  res$hotnet[(res$disease==traits_summ$MESH_ID[i] & res$gene %in% hotnet_res)] = 1
}

plot(euler(
  list(
    gwas=res$index[res$gwas==1],
    successes=res$index[res$success==1],
    hotnet=res$index[res$hotnet==1]
  )
), quantities=TRUE)

# all hotnet results - all

res$hotnet_full = 0

for(i in 1:dim(traits_summ)) {
  
  print(i)
  dirs = list.dirs(paste0("/Volumes/am673712/links/network_analysis/hotnet/hotnet/test_all_traits/results/output_", i, "/full_network-input_", i))[-1]
  
  hotnet_res_all = c()
  for(j in 1:length(dirs)) {
    res_delta = read_json(paste0(dirs[j], "/results.json"))
    res_delta = unique(unlist(res_delta$components))
    hotnet_res_all = c(hotnet_res_all, res_delta)
  }
  
  hotnet_res_all = unique(hotnet_res_all)
  res$hotnet_full[(res$disease==traits_summ$MESH_ID[i] & res$gene %in% hotnet_res_all)] = 1
  
  Sys.sleep(5)
  
}

plot(euler(
  list(
    gwas=res$index[res$gwas==1],
    successes=res$index[res$success==1],
    # hotnet=res$index[res$hotnet==1],
    hotnet_full=res$index[res$hotnet_full==1]
  )
), quantities=TRUE)

# hierarchical hotnet

res$hotnet_hierarchical = 0

for(i in 1:dim(traits_summ)[1]) {
  
  print(i)
  
  if(file.exists(paste0("/Volumes/am673712/links/network_analysis/hotnet/hierarchical_hotnet/test_all_traits/results/output_", i, "/results/clusters_metabase_plus_ppi_score_", i, ".tsv"))) {
    hier_res = read_csv(paste0("/Volumes/am673712/links/network_analysis/hotnet/hierarchical_hotnet/test_all_traits/results/output_", i, "/results/clusters_metabase_plus_ppi_score_", i, ".tsv"), skip=7, col_names=FALSE)
    hier_res = unlist(str_split(hier_res$X1, "\t"))
    res$hotnet_hierarchical[(res$disease==traits_summ$MESH_ID[i] & res$gene %in% hier_res)] = 1
  }
  
}

# complex data

res$complex = 0
net_complex = read_tsv("interactions/complex_portal/net_complex.tsv")

for(i in 1:dim(traits_summ)[1]) {
  
  hits = as.numeric(unlist(res[res$disease==traits_summ$MESH_ID[i] & res$gwas==1, 'entrez']))
  c_ix = apply(net_complex[,1:2], 1, function(x) any(x %in% hits))
  table(c_ix)
  if(all(!c_ix)) {
    next
  } else {
    res$complex[(res$disease==traits_summ$MESH_ID[i] & res$entrez %in% unique(unlist(net_complex[c_ix,1:2])))] = 1
  }
  
}

plot(euler(
  list(
    gwas=res$index[res$gwas==1],
    successes=res$index[res$success==1],
    hotnet=res$index[res$hotnet==1],
    # hotnet_full=res$index[res$hotnet_full==1],
    complex=res$index[res$complex==1]
  )
), quantities=TRUE)

# ligand receptor data

res$ligand_receptor = 0
lr_dat = read.table("interactions/ligand_receptor_db/LigandReceptor.tsv")

for(i in 1:dim(traits_summ)[1]) {
  
  hits = as.numeric(unlist(res[res$disease==traits_summ$MESH_ID[i] & res$gwas==1, 'entrez']))
  lr_ix = apply(dplyr::select(lr_dat, Ligand.Entrez, Receptor.Entrez), 1, function(x) any(x %in% hits))
  table(lr_ix)
  if(all(!lr_ix)) {
    next
  } else {
    res$ligand_receptor[(res$disease==traits_summ$MESH_ID[i] & res$entrez %in% unique(unlist(lr_dat[lr_ix,which(names(lr_dat) %in% c("Ligand.Entrez","Receptor.Entrez"))])))] = 1
  }
  
}

plot(euler(
  list(
    gwas=res$index[res$gwas==1],
    successes=res$index[res$success==1],
    hotnet=res$index[res$hotnet==1],
    # hotnet_full=res$index[res$hotnet_full==1],
    # complex=res$index[res$complex==1],
    lr=res$index[res$ligand_receptor==1]
  )
), quantities=TRUE)

# network neighbours

res$first_neighbors = 0

source("/GWD/appbase/projects/RD-TSci-Software/CB/packages/R/MetaBase.R")

trust = show.interaction.trusts()$trust_name[1:4]
mechanism = show.interaction.mechanisms()
mechanism = mechanism$mechanism_name[mechanism$direct]

for(i in 1:dim(traits_summ)[1]) {
  
  print(i)  
  hits = as.numeric(unlist(res[res$disease==traits_summ$MESH_ID[i] & res$gwas==1, 'entrez']))
  
  first_n = rbind(
    find_neighbours(hits, direction="downstream", trust=trust, mechanism=mechanism),
    find_neighbours(hits, direction="upstream", trust=trust, mechanism=mechanism)
  )
  
  first_n$gene_in = as.numeric(first_n$gene_in)
  first_n$gene_out = as.numeric(first_n$gene_out)
  
  x = first_n %>% dplyr::filter(direction=="downstream") %>% group_by(gene_in) %>% summarise(N=n(), gene=genesymbol_in[1]) %>% arrange(desc(N))
  table(x$gene_in %in% hits)
  x = first_n %>% dplyr::filter(direction=="upstream") %>% group_by(gene_out) %>% summarise(N=n(), gene=genesymbol_out[1]) %>% arrange(desc(N))
  table(x$gene_out %in% hits)
  
  first_n_res = data.frame(entrez=as.numeric(unique(unlist(first_n[,3:4]))))
  table(hits %in% first_n_res$entrez)
  
  res$first_neighbors[(res$disease==traits_summ$MESH_ID[i] & res$entrez %in% first_n_res$entrez)] = 1
  
}

# pathways

res$pathways = 0
metabase = gmtPathways("genesets/Metabase_GO_Maps_entrez.filt10.gmt")

for(i in 1:dim(traits_summ)[1]) {
  
  hits = as.numeric(unlist(res[res$disease==traits_summ$MESH_ID[i] & res$gwas==1, 'entrez']))
  
  p_ix = unlist(lapply(metabase, function(x) any(hits %in% x)))
  table(p_ix)
  pathway_hits = as.numeric(unique(unlist(metabase[p_ix])))
  res$pathways[(res$disease==traits_summ$MESH_ID[i] & res$entrez %in% pathway_hits)] = 1
  
}

venn.diagram(
  list(
    gwas=res$index[res$gwas==1],
    successes=res$index[res$success==1],
    # hotnet=res$index[res$hotnet==1]
    # hotnet_full=res$index[res$hotnet_full==1],
    # complex=res$index[res$complex==1]
    # lr=as.character(res$index[res$ligand_receptor==1]),
    hotnet_hier=res$index[res$hotnet_hierarchical==1]
    # pathway=res$index[res$pathways==1]
  ), filename="out.png", imagetype="png"
)


# SCORING -----------------------------------------------------------------

# add clinical failure

res$failure = 0

for(i in 1:dim(traits_summ)[1]) {
  
  print(i)
  failures = filter(tips_clin, `Clinical Label_PP`=="Clinical Failure", `MeSH ID`==traits_summ$MESH_ID[i]) %>% dplyr::select(`EntrezGene ID`) %>% unlist()
  res$failure[(res$disease==traits_summ$MESH_ID[i] & res$entrez %in% failures)] = 1
  
}

res_prop_or = res %>% group_by(disease) %>%
  summarise(
    gwas_success = (sum(success==1 & gwas==1)/sum(success==1 & gwas==0)) / (sum(failure==1 & gwas==1)/sum(failure==1 & gwas==0)),
    hotnet_success = (sum(success==1 & hotnet==1)/sum(success==1 & hotnet==0)) / (sum(failure==1 & hotnet==1)/sum(failure==1 & hotnet==0)),
    complex_success = (sum(success==1 & complex==1)/sum(success==1 & complex==0)) / (sum(failure==1 & complex==1)/sum(failure==1 & complex==0)),
    lr_success = (sum(success==1 & ligand_receptor==1)/sum(success==1 & ligand_receptor==0)) / (sum(failure==1 & ligand_receptor==1)/sum(failure==1 & ligand_receptor==0)),
    pathway_success = (sum(success==1 & pathways==1)/sum(success==1 & pathways==0)) / (sum(failure==1 & pathways==1)/sum(failure==1 & pathways==0))
  )
res_prop_or = res_prop_or[res_prop_or$gwas_success>0,]
res_prop_or = res_prop_or[!is.infinite(res_prop_or$gwas_success),]

res_prop = res %>% group_by(disease) %>%
  summarise(
    base_success = (sum(success) / (sum(success)+sum(failure))) * 100,
    gwas_success = sum(success==1 & gwas==1) / (sum(success==1 & gwas==1) + sum(failure==1 & gwas==1)) * 100,
    hotnet_success = sum(success==1 & hotnet==1) / (sum(success==1 & hotnet==1) + sum(failure==1 & hotnet==1)) * 100,
    complex_success = sum(success==1 & complex==1) / (sum(success==1 & complex==1) + sum(failure==1 & complex==1)) * 100,
    lr_success = sum(success==1 & ligand_receptor==1) / (sum(success==1 & ligand_receptor==1) + sum(failure==1 & ligand_receptor==1)) * 100,
    pathway_success = sum(success==1 & pathways==1) / (sum(success==1 & pathways==1) + sum(failure==1 & pathways==1)) * 100
  )
res_prop = res_prop[!is.nan(res_prop$gwas_success),]

res_plot_1 = as.data.frame(t(res_prop[,2:7]))
res_plot_2 = data.frame(label=names(res_prop)[2:7], mean=NA, ci_low=NA, ci_high=NA) 

res_plot_2$mean = apply(res_plot_1, 1, mean, na.rm=TRUE)
res_plot_2$ci_low = apply(res_plot_1, 1, function(x) mean(x, na.rm=TRUE) - (1.96*(sd(x, na.rm=TRUE)/sqrt(sum(!is.nan(x))))))
res_plot_2$ci_high = apply(res_plot_1, 1, function(x) mean(x, na.rm=TRUE) + (1.96*(sd(x, na.rm=TRUE)/sqrt(sum(!is.nan(x))))))

ggplot(res_plot_2, aes(y=str_replace(label, "_success", ""), x=mean)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis() + xlab("") + ylab("")

res_plot_1 = apply(res_plot_1, 2, function(x) x/x[1])
res_plot_1 = res_plot_1[-1,]
res_plot_3 = data.frame(label=names(res_prop)[3:7], odds=NA, ci_low=NA, ci_high=NA) 
res_plot_3$odds = apply(res_plot_1, 1, mean, na.rm=TRUE)
res_plot_3$ci_low = apply(res_plot_1, 1, function(x) mean(x, na.rm=TRUE) - (1.96*(sd(x, na.rm=TRUE)/sqrt(sum(!is.nan(x))))))
res_plot_3$ci_high = apply(res_plot_1, 1, function(x) mean(x, na.rm=TRUE) + (1.96*(sd(x, na.rm=TRUE)/sqrt(sum(!is.nan(x))))))

# *** add raw numbers to the y-axis
# *** control for there method (e.g. hotnet) is not returning results

ggplot(res_plot_3, aes(y=str_replace(label, "_success", ""), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis() + xlab("") + ylab("")

