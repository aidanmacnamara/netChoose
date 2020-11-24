
# LOAD --------------------------------------------------------------------

require(tidyverse)
require(biomaRt)
require(eulerr)
require(jsonlite)
require(readxl)
require(igraph)
require(sparklyr)


# MAPPING -----------------------------------------------------------------

mapping = read_excel("gwas/HumanGeneList_17Sep2018_workup_betterensembl_list.xlsx") # use mark's mapping file
names(mapping)[1:3] = c("ensembl_gene_id","entrezgene","hgnc_symbol")

# use karsten's mapping
mesh_to_analysis = read_excel("gwas/ukbb_ttam_finngen_phenotype_view_current_20190326.xlsx", sheet="UKBB") 

# clean up names
names(mesh_to_analysis) <- str_replace_all(names(mesh_to_analysis), " ", "_") %>% str_replace_all(., regex('[:punct:]'), "_") %>% str_replace_all(., regex('_{2,}'), "_") %>% str_to_lower()

# tidy and select relevant columns
mesh_to_analysis <- mesh_to_analysis %>% filter(!str_detect(analysis, "n/a")) %>% distinct(analysis, mesh_id)

# successful drug information
data_clin = read_excel("gwas/Pprojects_1Nov2018.xlsx") # fully mapped to entrez id
data_clin_summ = data_clin %>% filter(`Clinical Label_PP`=="Succeeded") %>% group_by(`MeSH ID`) %>% summarise(N=n()) %>% arrange(desc(N))

# fuzzy mesh mapping
fuzz = read_excel("gwas/PossibleMeSHrelationships_summary_rev.xlsx", sheet="PossibleMeSHrelationships")

# clean column names
names(fuzz) <- names(fuzz) %>% str_replace_all(., " ", "_") %>% str_replace_all(., regex('[:punct:]'), "_") %>% str_replace_all(., regex('_{2,}'), "_") %>% str_to_lower()

# filter for robust relationships (mark hurle suggested using 6 as a cutoff)
fuzz_filt <- fuzz %>% filter(relationship_class <= 6)

# export to cytoscape
fuzz_nodes = fuzz_filt[,1:2]; names(fuzz_nodes) = c("mesh","phenotype")
fuzz_nodes_2 = fuzz_filt[,3:4]; names(fuzz_nodes_2) = c("mesh","phenotype")
fuzz_nodes = rbind(fuzz_nodes, fuzz_nodes_2); fuzz_nodes = distinct(fuzz_nodes)
rm(fuzz_nodes_2)
write_tsv(fuzz_nodes, "cytoscape/fuzz_nodes.txt")

fuzz_edges = fuzz_filt %>% dplyr::select(-one_of(c("mesh_name_2","mesh_name_4")))
write_tsv(fuzz_edges, "cytoscape/fuzz_edges.txt")

# connect to gwas and clean
fuzz_to_gwas = data_clin %>% filter(`Clinical Label_PP`=="Succeeded" | `Clinical Label_PP`=="Clinical Failure") %>% inner_join(fuzz_filt, by=c("MeSH ID"="msh1")) %>% dplyr::select("EntrezGene ID", "Clinical Label_PP", "MeSH ID", msh2, relationship_class) %>% inner_join(mesh_to_analysis, by=c("msh2"="mesh_id")) %>% dplyr::select(`MeSH ID`,"msh2","relationship_class","analysis") %>% distinct()
names(fuzz_to_gwas)[1] = "msh1"


# GWAS --------------------------------------------------------------------

# full coloc results
coloc = read_tsv("~/Downloads/2019_03_12_ukb_all_colocs_overlap_GWAS.tsv") # all (karsten's table)

# karsten's computed coloc_hc
coloc_hc = read_tsv("~/Downloads/2019_03_12_ukb_positive_target_indications.tsv") # high-confidence

# filtered coloc_hc
# coloc_hc = filter(coloc, minpval2 <= 5e-8, minpval1 <= 1e-4, p12 >= 0.8)
# coloc_hc = coloc_hc %>% group_by(entity) %>% slice(which.max(p12))

# coloc edits
coloc_hc$entrez = mapping$entrezgene[match(coloc_hc$entity, mapping$ensembl_gene_id)]

# TMEM133  
coloc_hc$entity[which(coloc_hc$hgncid=="TMEM133")] = "ENSG00000165895"
coloc_hc$entrez[which(coloc_hc$hgncid=="TMEM133")] = 143872
coloc_hc$hgncid[which(coloc_hc$hgncid=="TMEM133")] = "ARHGAP42"

filter(coloc_hc, is.na(entrez)) %>% dplyr::select(hgncid) %>% distinct() # all entrez/ensembl mapped

coloc$entity[which(coloc$hgncid=="TMEM133")] = "ENSG00000165895"
coloc$hgncid[which(coloc$hgncid=="TMEM133")] = "ARHGAP42"

# filter
my_gwas = coloc_hc %>% group_by(analysis2) %>% summarise(N=n()) %>% arrange(desc(N)) %>% dplyr::select(analysis2)
names(my_gwas) = "gwas"
table(my_gwas$gwas %in% fuzz_to_gwas$analysis) # which gwas have a mesh id
my_gwas$gwas[!my_gwas$gwas %in% fuzz_to_gwas$analysis]
my_gwas$mesh = fuzz_to_gwas$msh1[match(my_gwas$gwas, fuzz_to_gwas$analysis)] # remove any without a mapping to mesh
my_gwas = my_gwas[!is.na(my_gwas$mesh),] # only keep those with a mesh mapping
my_gwas$disease = fuzz$mesh_name_2[match(my_gwas$mesh, fuzz$msh1)]


# RESULTS DATA FRAME ------------------------------------------------------

mapping_filt = mapping %>% filter(type_of_gene=="protein-coding") %>% dplyr::select(1:3)
mega_tbl = do.call("rbind", replicate(dim(my_gwas)[1], mapping_filt, simplify=FALSE))
mega_tbl$gwas = rep(my_gwas$gwas, each=dim(mapping_filt)[1])
mega_tbl$mesh = rep(my_gwas$mesh, each=dim(mapping_filt)[1])
mega_tbl$disease = rep(my_gwas$disease, each=dim(mapping_filt)[1])
mega_tbl$high_confidence_genetic=0; mega_tbl$success=0; mega_tbl$failure=0; mega_tbl$index=1:dim(mega_tbl)[1]

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  successes = filter(data_clin, `Clinical Label_PP`=="Succeeded", `MeSH ID` %in% unique(fuzz_to_gwas$msh2[which(fuzz_to_gwas$msh1 == my_gwas$mesh[i])])) %>% dplyr::select(`EntrezGene ID`) %>% unlist()
  failures =  filter(data_clin, `Clinical Label_PP`=="Clinical Failure", `MeSH ID` %in% unique(fuzz_to_gwas$msh2[which(fuzz_to_gwas$msh1 == my_gwas$mesh[i])])) %>% dplyr::select(`EntrezGene ID`) %>% unlist()
  
  gwas = filter(coloc_hc, analysis2==my_gwas$gwas[i]) %>% dplyr::select(entity) %>% unlist()
  
  # take directly from full coloc table
  # gwas = filter(coloc, analysis2==my_gwas$gwas[i], p12 >= 0.8) %>% dplyr::select(entity) %>% unlist() %>% as.character()
  
  gwas = gwas[!is.na(gwas)]; gwas = unique(gwas)
  
  # add randomizer
  # gwas = sample(x=filter(mapping, type_of_gene=="protein-coding") %>% dplyr::select(ensembl_gene_id) %>% unlist() %>% as.character(), size=length(gwas), replace=FALSE)
  
  mega_tbl$success[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% successes)] = 1
  mega_tbl$failure[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% failures)] = 1
  mega_tbl$high_confidence_genetic[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$ensembl_gene_id %in% gwas)] = 1
  
}

# check numbers
coloc_hc_summ = coloc_hc %>% group_by(analysis2) %>% summarise(N=n()) %>% arrange(desc(N))
mega_tbl_n = mega_tbl %>% filter(high_confidence_genetic==1) %>% group_by(gwas) %>% summarise(N=n()) %>% arrange(desc(N))

# total number of hcghs
sum(mega_tbl_n$N) # 14374

# total number of drug targets
data_clin %>% filter((`Clinical Label_PP`=="Succeeded"|`Clinical Label_PP`=="Clinical Failure"), `MeSH ID` %in% unique(fuzz_to_gwas$msh2[fuzz_to_gwas$msh1 %in% my_gwas_nsp$mesh])) %>% distinct(`Target|Indication`) # 1268

data.frame(coloc=coloc_hc_summ$N, mega_tbl=mega_tbl_n$N[match(coloc_hc_summ$analysis2, mega_tbl_n$gwas)]) %>% ggplot(aes(x=coloc, y=mega_tbl)) + geom_point() + theme_thesis()


# ADD GENE SCORE ----------------------------------------------------------

mega_tbl$h4 = 0
mega_tbl$h4_log = 0

# require(odbc)
# conn = dbConnect(odbc::odbc(), dsn="impaladsn")
# dbListFields(conn, "2019_03_12_ukb_all_target_indications_no_filter_max_p12")
# test_sql = dbSendQuery(conn, "SELECT * FROM `2019_03_12_ukb_all_target_indications_no_filter_max_p12` WHERE analysis2 = 'GSK500kV3_FEV1_maximumValue'")

# conf <- sparklyr::spark_config()
# conf$sparklyr.log.console <- FALSE
# conf$spark.executor.memory <- "4g"
# conf$spark.yarn.am.memory <- "4g" 
# conf$spark.driver.maxResultSize <- "4g" 

# connect to full coloc results
# sc <- spark_connect(
#   master = "yarn-client",
#   spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
#   version = "2.3",
#   config = conf
# ) 

# sc %>% tbl_change_db("am673712")
# src_tbls(sc)
# coloc <- tbl(sc, "2019_03_12_ukb_all_target_indications_no_filter_max_p12")

# system.time({
#   coloc_head = coloc %>% head() %>% collect()
# })

for(i in 1:dim(my_gwas)[1]) {
  
  print(i)
  
  # get the minimum eqtl score per tissue
  # multiple gwas results here
  coloc_filt = filter(coloc, analysis2==my_gwas$gwas[i]) %>% group_by(entity) %>% dplyr::slice(which.max(p12))
  coloc_filt_hits = paste(coloc_filt$entity, coloc_filt$analysis2, sep="_")
  hit_ix = match(coloc_filt_hits, paste(mega_tbl$ensembl_gene_id, mega_tbl$gwas, sep="_"))
  
  # some not matching because id not present in mega_tbl
  na_ix = which(is.na(hit_ix))
  
  if(!is_empty(na_ix)) {
    hit_ix = hit_ix[!is.na(hit_ix)] # remove na indexes
    mega_tbl$h4[hit_ix] = coloc_filt$p12[-na_ix]
    mega_tbl$h4_log[hit_ix] = -log(1-coloc_filt$p12[-na_ix], base=2)
  } else {
    mega_tbl$h4[hit_ix] = coloc_filt$p12
    mega_tbl$h4_log[hit_ix] = -log(1-coloc_filt$p12, base=2)
  }
  
}


# HOTNET SETUP ------------------------------------------------------------

# gene_ix = read_tsv("/Volumes/am673712/links/network_analysis/interactions/metabase_plus_ppi/metabase_plus_ppi_gene_index.txt", col_names=FALSE)
load("r_data/gene_ix.RData")

for(i in 1:dim(my_gwas)[1]) {
  
  print(i)
  
  heat_scores = data.frame(gene=gene_ix$X2, score=NA)
  
  coloc_filt = filter(coloc, analysis2==my_gwas$gwas[i]) %>% group_by(entity) %>% dplyr::slice(which.max(p12))
  
  # add heat scores
  heat_scores$score = -log(1-coloc_filt$p12[match(heat_scores$gene, coloc_filt$hgncid)], base=2)
  heat_scores$score[is.na(heat_scores$score)] = 0
  heat_scores = arrange(heat_scores, desc(score))
  write_tsv(heat_scores, paste0("/Volumes/am673712/links/network_analysis/hotnet/hotnet/test_all_gwas/input/metabase_plus_ppi/input_", i, ".tsv"), col_names=FALSE)
  
}

# run on slurm
# sbatch -n 1 -c 12 --mem=50G -t 360 -o slurm_log.out run_test.sh


# OTHER NETWORK SOURCES ---------------------------------------------------

networks = c("omnipath","huri","string","intomics")

for(i in 1:dim(my_gwas)[1]) {
  
  print(i)
  coloc_filt = filter(coloc, analysis2==my_gwas$gwas[i]) %>% group_by(entity) %>% dplyr::slice(which.max(p12))
  coloc_filt$entrez = mapping$entrezgene[match(coloc_filt$entity, mapping$ensembl_gene_id)]
  
  for(j in 4:length(networks)) {
    
    # load gene score
    gene_ix = read_tsv(paste0("interactions/",networks[j],"/",networks[j],"_nodes.txt"), col_names=FALSE)
    
    # add heat scores
    heat_scores = data.frame(gene=gene_ix$X2, score=NA)
    heat_scores$score = -log(1-coloc_filt$p12[match(heat_scores$gene, coloc_filt$entrez)], base=2)
    heat_scores$score[is.na(heat_scores$score)] = 0
    heat_scores = arrange(heat_scores, desc(score))
    write_tsv(heat_scores, paste0("/Volumes/am673712/links/network_analysis/hotnet/hotnet/test_all_gwas/input/", networks[j], "/input_", i, ".tsv"), col_names=FALSE)
    
  }
  
}


# OVERLAPS ----------------------------------------------------------------

plot(euler(
  list(
    `high confidence genetics` = mega_tbl$index[mega_tbl$high_confidence_genetic==1],
    successes = mega_tbl$index[mega_tbl$success==1],
    failures = mega_tbl$index[mega_tbl$failure==1]
  )
), quantities=TRUE)

# remove target/gwas pairs where success==1 and failure==1 - this results from fuzzy mesh matching
mega_tbl = mega_tbl %>% filter(!(success==1 & failure==1))

plot(euler(
  list(
    `high confidence genetics` = mega_tbl$index[mega_tbl$high_confidence_genetic==1],
    successes = mega_tbl$index[mega_tbl$success==1],
    failures = mega_tbl$index[mega_tbl$failure==1]
  )
), quantities=TRUE)


# COMPLEX DATA ------------------------------------------------------------

mega_tbl$complex = 0
net_complex = read_tsv("interactions/complex_portal/net_complex.tsv")

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  hits = as.numeric(unlist(mega_tbl[mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$high_confidence_genetic==1, 'entrezgene']))
  c_ix = apply(net_complex[,1:2], 1, function(x) any(x %in% hits))
  table(c_ix)
  if(all(!c_ix)) {
    next
  } else {
    mega_tbl$complex[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% unique(unlist(net_complex[c_ix,1:2])))] = 1
  }
  
}

plot(euler(
  list(
    `high confidence genetics` = mega_tbl$index[mega_tbl$high_confidence_genetic==1],
    successes = mega_tbl$index[mega_tbl$success==1],
    failures = mega_tbl$index[mega_tbl$failure==1],
    complex = mega_tbl$index[mega_tbl$complex==1]
  )
), quantities=TRUE)


# LIGAND RECEPTOR ---------------------------------------------------------

mega_tbl$ligand_receptor = 0
mark_proxy_filt = read_tsv("interactions/ligand_receptor_db/mark_proxy_filt.txt")

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  hits = as.numeric(unlist(mega_tbl[mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$high_confidence_genetic==1, 'entrezgene']))
  lr_ix = apply(mark_proxy_filt, 1, function(x) any(x %in% hits))
  table(lr_ix)
  if(all(!lr_ix)) {
    next
  } else {
    mega_tbl$ligand_receptor[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% unique(unlist(mark_proxy_filt[lr_ix,])))] = 1
  }
  
}

plot(euler(
  list(
    `high confidence genetics` = mega_tbl$index[mega_tbl$high_confidence_genetic==1],
    successes = mega_tbl$index[mega_tbl$success==1],
    failures = mega_tbl$index[mega_tbl$failure==1],
    complex = mega_tbl$index[mega_tbl$complex==1],
    ligand_receptor = mega_tbl$index[mega_tbl$ligand_receptor==1]
  )
), quantities=TRUE)

method_overlaps =   list(
  `high confidence genetics` = mega_tbl$index[mega_tbl$high_confidence_genetic==1],
  successes = mega_tbl$index[mega_tbl$success==1],
  failures = mega_tbl$index[mega_tbl$failure==1],
  complex = mega_tbl$index[mega_tbl$complex==1],
  ligand_receptor = mega_tbl$index[mega_tbl$ligand_receptor==1]
)

my_dat_genes_all_tbl = as.data.frame.matrix((table(stack(method_overlaps))))
my_dat_genes_all_tbl = cbind(rownames(my_dat_genes_all_tbl), my_dat_genes_all_tbl)
rownames(my_dat_genes_all_tbl) = NULL
upset(my_dat_genes_all_tbl, order.by="freq")


# NETWORK NEIGHBOURS ------------------------------------------------------

find_network_second_neighbor = TRUE # boolean: true for second neighbors, false for first neighbors

mega_tbl$network_first_neighbor = 0
mega_tbl$network_second_neighbor = 0

mb = read_tsv("interactions/metabase_plus_ppi/metabase/MetaBaseFullNetworkHML.txt")
mb_filt = mb %>% filter(TrustLevel=="high", !Mechanism %in% c("Competition","Influence on expression","Pharmacological effect", "Toxic effect","Unspecified"))
mb_graph = graph_from_data_frame(mb_filt[,c(1,3,5:11)])

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  hits = as.numeric(unlist(mega_tbl[mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$high_confidence_genetic==1, 'entrezgene']))
  
  if(is_empty(hits)) next
  
  my_neighbors = c()
  
  for(j in 1:length(hits)) {
    possible_error <- tryCatch(
      {
        if(find_network_second_neighbor) {
          neighborhood(mb_graph, order=2, nodes=hits[j])
        } else {
          neighbors(mb_graph, hits[j], mode="all")    
        }
      },
      error=function(e) e
    )
    if(inherits(possible_error, "error")) next
    if(find_network_second_neighbor) {
      my_neighbors = c(my_neighbors, as.numeric(possible_error[[1]]$name))
    } else {
      my_neighbors = c(my_neighbors, as.numeric(possible_error$name))
    }
  }
  
  my_neighbors = unique(my_neighbors)
  
  if(find_network_second_neighbor) {
    mega_tbl$network_second_neighbor[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% my_neighbors)] = 1
  } else {
    mega_tbl$network_first_neighbor[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% my_neighbors)] = 1
  }
  
}


# PATHWAYS ----------------------------------------------------------------

mega_tbl$pathways = 0
metabase = gmtPathways("genesets/Metabase_GO_Maps_entrez.filt10.gmt")

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  hits = as.numeric(unlist(mega_tbl[mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$high_confidence_genetic==1, 'entrezgene']))
  
  p_ix = unlist(lapply(metabase, function(x) any(hits %in% x)))
  table(p_ix)
  pathway_hits = as.numeric(unique(unlist(metabase[p_ix])))
  mega_tbl$pathways[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% pathway_hits)] = 1
  
}


# PATHWAY NEIGHBOURS ------------------------------------------------------

load("r_data/Metabase_interactors.Rdata")

mega_tbl$pathway_first_neighbor = 0
mega_tbl$pathway_second_neighbor = 0

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  hits = as.numeric(unlist(mega_tbl[mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$high_confidence_genetic==1, 'entrezgene']))
  
  p_ix = which(as.numeric(names(gg.interactors)) %in% hits)
  first_hits = as.numeric(unlist(lapply(gg.interactors[p_ix], function(x) x[[1]])))
  second_hits = as.numeric(unlist(lapply(gg.interactors[p_ix], function(x) x[[2]])))
  
  mega_tbl$pathway_first_neighbor[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% first_hits)] = 1
  mega_tbl$pathway_second_neighbor[(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% second_hits)] = 1
  
}


# NETWORK PROPAGATION -----------------------------------------------------

my_networks = c(
  "intomics_pascal_lfdr_0.01",
  "intomics_pascal_lfdr_0.05",
  "intomics_pascal_lfdr_0.1"
)

missing_dat = data.frame(matrix(0, nrow=length(my_gwas$mesh), ncol=length(my_networks)))
names(missing_dat) = my_networks

for(j in 1:length(my_networks)) {
  
  mega_tbl = mega_tbl %>% mutate(!!my_networks[j] := 0)
  
  for(i in 1:length(my_gwas$mesh)) {
    
    print(i)
    
    if(file.exists(paste0("hotnet/hotnet/test_all_gwas/results/", my_networks[j], "/output_", i, "/consensus/subnetworks.tsv"))) {
      hotnet_res = read_tsv(paste0("hotnet/hotnet/test_all_gwas/results/", my_networks[j], "/output_", i, "/consensus/subnetworks.tsv"), col_names=FALSE)
      hotnet_res = as.numeric(unlist(str_split(hotnet_res$X1[3:length(hotnet_res$X1)], " ")))
      mega_tbl[[my_networks[j]]][(mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$entrezgene %in% hotnet_res)] = 1
    } else {
      missing_dat[[my_networks[j]]][i] = 1
    }
    
  }
}

# RANDOM ------------------------------------------------------------------

mega_tbl = mega_tbl %>% group_by(gwas) %>% mutate(random = ifelse(1:n() %in% sample(1:n(),10000), 1, 0))

# check pathways or < 1
metabase = gmtPathways("genesets/Metabase_GO_Maps_entrez.filt10.gmt")
metabase = sort(as.numeric(unlist(metabase)))
mega_tbl = mega_tbl %>% group_by(gwas) %>% mutate(random = ifelse(entrezgene %in% sample(metabase,10000,replace=FALSE), 1, 0))


# OMNIPATH HEAT -----------------------------------------------------------

mega_tbl$hotnet_omnipath_heat = NA

for(i in 1:length(my_gwas$mesh)) {
  
  print(i)
  
  if(file.exists(paste0("hotnet/hotnet/test_all_gwas/results/omnipath/output_", i, "/heat.txt"))) {
    heat = read_tsv(paste0("hotnet/hotnet/test_all_gwas/results/omnipath/output_", i, "/heat.txt"), col_names=FALSE)
  } else {
    print(paste("File:", i, "doesn't exist ..."))
    next
  }
  
  # pull out the index and the entrez id
  
  match_dat = mega_tbl %>% filter(gwas==my_gwas$gwas[i]) %>% dplyr::select(entrezgene, index)
  heat_ix = match(heat$X1, match_dat$entrezgene)
  missing_ix = which(is.na(heat_ix))
  heat_ix = heat_ix[-missing_ix]
  
  # update using index
  
  mega_tbl$hotnet_omnipath_heat[match_dat$index[heat_ix]] = heat$X2[-missing_ix]
  
}


# USE PASCAL GENE SCORES FOR HOTNET ---------------------------------------

require(twilight)

# define a cutoff for the pascal scores

networks = c("intomics")
lfdr_range = c(0.1, 0.05, 0.01)
for(j in 1:length(networks)) {
  for(k in 1:length(lfdr_range)) {
    if(!dir.exists(paste0("hotnet/hotnet/test_all_gwas/input/", networks[j], "_pascal_lfdr_", lfdr_range[k]))) dir.create(paste0("hotnet/hotnet/test_all_gwas/input/", networks[j], "_pascal_lfdr_", lfdr_range[k]))
  }
}

for(i in 1:dim(my_gwas)[1]) {
  
  hits_filt = mega_tbl %>% filter(gwas==my_gwas$gwas[i])
  if(all(is.na(hits_filt$pascal_default))) { # no pascal data for this gwas
    print(paste("No Pascal data for GWAS:", my_gwas$gwas[i]))
    next
  }
  
  for(j in 1:length(networks)) {
    
    # load gene score
    gene_ix = read_tsv(paste0("interactions/",networks[j],"/",networks[j],"_nodes.txt"), col_names=FALSE)
    
    # add heat scores
    heat_scores = data.frame(gene=gene_ix$X2, score=NA)
    
    # get the gene p values
    heat_scores$p = hits_filt$pascal_default[match(heat_scores$gene, hits_filt$entrezgene)]
    
    # find the lfdr for the p values
    heat_scores$lfdr = NA
    t_res = tbl_df(twilight(heat_scores$p[!is.na(heat_scores$p)])$result)
    t_res = t_res %>% arrange(index)
    
    for(k in 1:length(lfdr_range)) {
      
      heat_scores_copy = heat_scores
      heat_scores_copy$lfdr[!is.na(heat_scores_copy$p)] = t_res$fdr
      
      # heat_scores_copy %>% ggplot(aes(x=p, y=1-lfdr)) + geom_line() + theme_thesis() + xlab("Pascal p-value") + ylab("1-FDR")
      
      heat_scores_copy$score = -log(hits_filt$pascal_default[match(heat_scores_copy$gene, hits_filt$entrezgene)], base=10)
      heat_scores_copy$score[heat_scores_copy$lfdr > lfdr_range[k]] = 0
      heat_scores_copy$score[is.na(heat_scores_copy$score)] = 0
      heat_scores_copy = arrange(heat_scores_copy, desc(score))
      write_tsv(heat_scores_copy[,1:2], paste0("hotnet/hotnet/test_all_gwas/input/", networks[j], "_pascal_lfdr_", lfdr_range[k], "/input_", i, ".tsv"), col_names=FALSE)
      
    }
  }
  
}


