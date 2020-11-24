require(jsonlite)
require(readxl)
require(igraph)
require(yardstick)
require(plotROC)

mega_tbl = read_tsv("mega_tbl.txt")


# SCORING -----------------------------------------------------------------

# filter out the gwas with no pascal results

my_gwas_has_pascal = mega_tbl %>% group_by(gwas) %>% summarise(pascal=any(!is.na(pascal_default)))
table(my_gwas_has_pascal$pascal)

# take out the 'surgical procedures' gwas and missing pascal (oa) gwas

gwas_filt_ix = !grepl("surgical", my_gwas$disease, ignore.case=TRUE) & my_gwas_has_pascal$pascal[match(my_gwas$gwas, my_gwas_has_pascal$gwas)]
table(gwas_filt_ix)

my_gwas_filt = my_gwas[gwas_filt_ix,]

# this leaves 648 gwas

mega_tbl_filt = filter(mega_tbl, gwas %in% my_gwas_filt$gwas)

# check result with only pascal sig genes
pascal_true_false = mega_tbl$gwas %in% my_gwas_filt | mega_tbl$pascal_default < 0.05
# mega_tbl_filt_pascal = mega_tbl %>% filter(pascal_default < 0.05)
# saved data will be 'res_plot_without_hcg_with_0s_pascal_filtered' and 'method_avg_without_hcg_pascal_filtered'

# group by mesh
# mega_tbl_filt = mega_tbl %>% filter(gwas %in% my_gwas_filt$gwas) %>% dplyr::select(-gwas) %>% group_by(ensembl_gene_id,entrezgene,hgnc_symbol,mesh,disease) %>% summarise_all(funs(max)) %>% arrange(disease, hgnc_symbol)
# mega_tbl_filt$gwas = mega_tbl_filt$mesh

methods_tested = c("high_confidence_genetic","complex","ligand_receptor","network_first_neighbor","pathways","pathway_first_neighbor","pathway_second_neighbor","random","hotnet_metabase_plus_ppi","hotnet_omnipath","network_second_neighbor","hotnet_huri","hotnet_string","hotnet_intomics","intomics_pascal_lfdr_0.01","intomics_pascal_lfdr_0.05","intomics_pascal_lfdr_0.1")

res_main_without_seed = mega_tbl_out(mega_tbl_filt, col_names=methods_tested, include_seed=FALSE)
res_main_with_seed = mega_tbl_out(mega_tbl_filt, col_names=methods_tested, include_seed=TRUE) # run with 'include_seed=true' to get the hcgh or

res_plot_cmh = res_main_without_seed$res_plot_cmh # maintain name for markdown
res_plot_cmh = rbind(
  res_plot_cmh,
  res_main_with_seed$res_plot_cmh[res_main_with_seed$res_plot_cmh$method=="high_confidence_genetic",]
)
res_plot_cmh = arrange(res_plot_cmh, odds)

level_order = res_plot_cmh$method
plot_a = ggplot(res_plot_cmh, aes(y=factor(method, levels=level_order), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(25) + xlab("") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a


# METHOD COLLATION --------------------------------------------------------

method_avg = method_delta(mega_tbl_filt, col_names=methods_tested, my_gwas=my_gwas_filt)
# method_avg$Hits[method_avg$method=="high_confidence_genetic"] = 0

plot_b = method_avg %>% ggplot(aes(x=factor(method,levels=level_order), y=log(Hits, base=10))) + geom_bar(stat="identity") + theme_thesis(10) + coord_flip() + xlab("") + ylab("Average Number of Additional Hits (Log10)")
plot_b

plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# FIRST NEIGHBORS ALL NETWORKS --------------------------------------------

mega_tbl_first_neighbors = mega_tbl %>% dplyr::select(1:10)
find_network_second_neighbor = FALSE # boolean: true for second neighbors, false for first neighbors

col_names = c("metabase_plus_ppi","string","omnipath","huri","intomics")
mega_tbl_first_neighbors = tbl_df(cbind(mega_tbl_first_neighbors, matrix(0, nrow=dim(mega_tbl_first_neighbors)[1], ncol=length(col_names))))
names(mega_tbl_first_neighbors)[11:15] = col_names

for(i in 1:length(col_names)) {
  
  print(paste("Network:", col_names[i]))
  mega_tbl_ix = which(names(mega_tbl_first_neighbors)==col_names[i])
  net = read_tsv(paste0("interactions/", col_names[i], "/", col_names[i], "_entrez.txt")) 
  net = graph_from_data_frame(net)
  
  for(j in 1:length(my_gwas$mesh)) {
    
    print(paste0("GWAS ", my_gwas$gwas[j], ", index: ", j))
    
    hits = as.numeric(unlist(mega_tbl_first_neighbors[mega_tbl_first_neighbors$gwas==my_gwas$gwas[j] & mega_tbl_first_neighbors$high_confidence_genetic==1, 'entrezgene']))
    
    if(is_empty(hits)) next
    
    my_neighbors = c()
    
    for(k in 1:length(hits)) {
      possible_error <- tryCatch(
        {
          if(find_network_second_neighbor) {
            neighborhood(net, order=2, nodes=hits[k])
          } else {
            neighbors(net, hits[k], mode="all")    
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
    mega_tbl_first_neighbors[(mega_tbl_first_neighbors$gwas==my_gwas$gwas[j] & mega_tbl_first_neighbors$entrezgene %in% my_neighbors),mega_tbl_ix] = 1
    
  }
}

for(which_neighbors in 1:2) { # get data from first or second neighbors
  
  if(which_neighbors==1) mega_tbl_neighbors = read_tsv("mega_tbl_first_neighbors.txt")
  if(which_neighbors==2) mega_tbl_neighbors = read_tsv("mega_tbl_second_neighbors.txt")
  
  # mega_tbl_neighbors = filter(mega_tbl_neighbors, gwas %in% my_gwas_filt$gwas)
  mega_tbl_neighbors = filter(mega_tbl_neighbors, pascal_true_false)
  
  methods_tested = col_names
  res_neighbors = mega_tbl_out(mega_tbl_neighbors, col_names=methods_tested, include_seed=FALSE)
  res_neighbors = res_neighbors$res_plot_cmh
  if(which_neighbors==1) res_neighbors$method = paste0(res_neighbors$method, "_first_neighbors")
  if(which_neighbors==2) res_neighbors$method = paste0(res_neighbors$method, "_second_neighbors")
  res_plot_cmh = rbind(res_plot_cmh, res_neighbors) # combine with the main results
  
  method_avg_neighbors = method_delta(mega_tbl_neighbors, col_names=methods_tested, my_gwas=my_gwas_filt)
  if(which_neighbors==1) method_avg_neighbors$method = paste0(method_avg_neighbors$method, "_first_neighbors")
  if(which_neighbors==2) method_avg_neighbors$method = paste0(method_avg_neighbors$method, "_second_neighbors")
  method_avg = rbind(method_avg, method_avg_neighbors) # combine with the main results
  
}

save(res_plot_cmh, file="final_report/res_plot_without_hcg_with_0s_pascal_filtered.RData") # savepoint
save(method_avg, file="final_report/method_avg_without_hcg_pascal_filtered.RData") # savepoint

level_order = res_plot_cmh$method
plot_a = ggplot(res_plot_cmh, aes(y=factor(method, levels=level_order), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(25) + xlab("") + ylab("") # + coord_cartesian(xlim=c(0,5))
plot_a

plot_b = method_avg %>% ggplot(aes(x=factor(method,levels=level_order), y=log(Hits, base=10))) + geom_bar(stat="identity") + theme_thesis(10) + coord_flip() + xlab("") + ylab("Average Number of Hits per GWAS (Log10)")
plot_b

plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# ASSESS PASCAL SIGNIFICANT GENE SETS -------------------------------------

load("gwas_698_sigGenesets_sigGenes_immunonavigator.RData")
load("gwas_650_ii_pathways.RData")
load("gwas_650_sigGenesets_sigGenes.RData")
load("gwas_650_sigGenesets_sigGenes_magma.RData")
load("gwas_650_sigGenesets_sigGenes_ghits.RData")

pascal_cell_type_combined = vector("list", length(pascal.cd4.sig.genes.clean))
names(pascal_cell_type_combined) = names(pascal.cd4.sig.genes.clean)
all(names(pascal_cell_type_combined)==names(pascal.mph.sig.genes.clean))
for(i in 1:length(pascal_cell_type_combined)) {
  if(is_empty(pascal.cd4.sig.genes.clean[[i]]) & is_empty(pascal.mph.sig.genes.clean[[i]])) next
  pascal_cell_type_combined[[i]] = c(pascal.cd4.sig.genes.clean[[i]],pascal.mph.sig.genes.clean[[i]])
}
all(names(pascal_cell_type_combined)==names(pascal.mph.sig.genes.clean))

gb_list = list(
  pascal.cd4.sig.genes.clean = pascal.cd4.sig.genes.clean,
  pascal.mph.sig.genes.clean = pascal.mph.sig.genes.clean,
  pascal.metabase.sig.genes.clean = pascal.metabase.sig.genes.clean,
  pascal.ppi.sig.genes.clean = pascal.ppi.sig.genes.clean,
  pascal.ppi.sig.genes.ghits = pascal.ppi.sig.genes.ghits,
  pascal.coexp.sig.genes.clean = pascal.coexp.sig.genes.clean,
  pascal.coexp.sig.genes.ghits = pascal.coexp.sig.genes.ghits,
  pascal.reactome.sig.genes.clean = pascal.reactome.sig.genes.clean,
  pascal.reactome.sig.genes.ghits = pascal.reactome.sig.genes.ghits,
  pascal_cell_type_combined = pascal_cell_type_combined,
  magma.coexp.sig.genes.clean = magma.coexp.sig.genes.clean,
  magma.metabase.sig.genes.clean = magma.metabase.sig.genes.clean,
  magma.ppi.sig.genes.clean = magma.ppi.sig.genes.clean,
  magma.reactome.sig.genes.clean = magma.reactome.sig.genes.clean
)

# initialise new mega_tbl

mega_tbl_gb = mega_tbl_filt %>% dplyr::select(1:10)
mega_tbl_gb = tbl_df(cbind(mega_tbl_gb, matrix(0, nrow=dim(mega_tbl_gb)[1], ncol=length(gb_list))))
names(mega_tbl_gb)[11:24] = names(gb_list)

for(i in 1:length(gb_list)) {
  
  mega_tbl_ix = which(names(mega_tbl_gb)==names(gb_list)[i])
  print(paste("Table ix:", mega_tbl_ix))
  total_hits = 0
  
  for(j in 1:length(my_gwas_filt$gwas)) {
    
    if(!my_gwas_filt$gwas[j] %in% names(gb_list[[i]])) next
    
    hits = as.numeric(unlist(gb_list[[i]][[which(names(gb_list[[i]])==my_gwas_filt$gwas[j])]]))
    
    if(is_empty(hits)) {
      # print(paste("Skipping index:", j, "Trait:", traits_109$gwas[j]))
      next
    } else {
      # print(paste0("Index: ", j, ", Trait: ", traits_109$gwas[j], ", N: ", length(hits)))
      mega_tbl_gb[(mega_tbl_gb$gwas==my_gwas_filt$gwas[j] & mega_tbl_gb$entrezgene %in% hits),mega_tbl_ix] = 1
      total_hits = total_hits + length(hits)
    }
    
  }
}

methods_tested = names(gb_list)

res_pascal = mega_tbl_out(mega_tbl_gb, col_names=methods_tested, include_seed=FALSE)
res_plot_cmh_pascal = res_pascal$res_plot_cmh
save(res_plot_cmh_pascal, file="final_report/res_plot_pascal_without_hcg_with_0s.RData") # savepoint

level_order = res_pascal$res_plot_cmh$method
plot_a = ggplot(res_pascal$res_plot_cmh, aes(y=factor(method, levels=level_order), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("") + ylab("") + coord_cartesian(xlim=c(0,14))
plot_a

method_avg_pascal = method_delta(mega_tbl=mega_tbl_gb, col_names=methods_tested, my_gwas=my_gwas_filt)
method_avg_pascal = method_avg_pascal[method_avg_pascal$method %in% level_order,]
save(method_avg_pascal, file="final_report/method_avg_pascal_without_hcg.RData") # savepoint

plot_b = method_avg_pascal %>% ggplot(aes(x=factor(method,levels=level_order), y=Hits)) + geom_bar(stat="identity") + theme_thesis(10) + coord_flip() + xlab("") + ylab("Average Number of Hits per GWAS")
plot_b

plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))

# immune signature of each gwas

to_plot = unlist(pascal.metabase.genes.unionii.ordered.list)
to_plot = tbl_df(data.frame(gwas=names(to_plot), immune_fraction=to_plot, row.names=NULL))
to_plot$gwas = str_remove(to_plot$gwas, "[0-9]+$")
ggplot(to_plot, aes(x=reorder(gwas, immune_fraction, FUN=median), y=immune_fraction)) + geom_boxplot()

to_plot_summ = to_plot %>% group_by(gwas) %>% summarise(median=median(immune_fraction)) %>% arrange(desc(median))

immune_gwas = to_plot_summ %>% filter(median > 0.8) %>% select(gwas) %>% unlist() %>% as.character()

res_pascal_immune = mega_tbl_out(filter(mega_tbl_gb, gwas %in% immune_gwas), col_names=methods_tested, include_seed=FALSE)
level_order = res_pascal_immune$res_plot_cmh$method
plot_a = ggplot(res_pascal_immune$res_plot_cmh, aes(y=factor(method, levels=level_order), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("") + ylab("") + coord_cartesian(xlim=c(0,14))
plot_a


# NEW ENRICHMENT HITS -----------------------------------------------------

# this section has some worked examples of where proxy enrichment occurs

hits_complex = mega_tbl %>% filter(complex==1, high_confidence_genetic==0, success==1)

# filter the mega_tbl on the complex data and the example

complex_portal = read_tsv("interactions/complex_portal/net_complex.tsv")
complex_nodes = unique(unlist(complex_portal[,1:2]))
mega_tbl_complex = mega_tbl %>% filter(gwas=="GSK500KV3lin_M2W_psoriasis", entrezgene %in% complex_nodes)
write_tsv(mega_tbl_complex, "tmp/mega_tbl_complex.txt")

# hotnet omnipath

hits_hotnet = mega_tbl %>% filter(hotnet_omnipath==1, high_confidence_genetic==0, success==1)

# pick out the example of TGFB1 vs. osteoarthritis

hits_hotnet %>% filter(hgnc_symbol=="TGFB1", disease=="Osteoarthritis") %>% View()

# get the network

net = read_tsv("interactions/omnipath/omnipath_entrez.tsv")
nodes = read_tsv("interactions/omnipath/omnipath_nodes.txt", col_names=FALSE)
net = graph_from_data_frame(net[,1:2])

# get hotnet results

res_ix = which(my_gwas$gwas=="GSK500KV3lin_M2W_osteoarthritis")

hotnet_res = read_tsv(paste0("/Volumes/am673712/links/network_analysis/hotnet/hotnet/test_all_gwas/results/omnipath/output_", res_ix, "/consensus/subnetworks.tsv"), col_names=FALSE)

hotnet_res = as.numeric(unlist(str_split(hotnet_res$X1[3:length(hotnet_res$X1)], " ")))

table(hotnet_res %in% V(net)$name)

net_filt = subgraph(net, v=as.character(hotnet_res))
write_tsv(as.data.frame(get.edgelist(net_filt)), "tmp/omnipath_edges.txt")
mega_tbl_omnipath = mega_tbl %>% filter(gwas=="GSK500KV3lin_M2W_osteoarthritis", entrezgene %in% V(net_filt)$name)
write_tsv(mega_tbl_omnipath , "tmp/mega_tbl_omnipath.txt")

mega_tbl %>% filter(gwas=="GSK500KV3lin_M2W_osteoarthritis", high_confidence_genetic==1) %>% View()

# colocalisation data for tgfb1 and osteoarthritis

coloc_tgfb1 = coloc %>% filter(grepl("osteo", analysis2, ignore.case=TRUE), entity=="ENSG00000105329") 
table(coloc_tgfb1$analysis2)

coloc_tgfb1 %>% filter(analysis2=="GSK500KV3lin_M2W_osteoarthritis") %>% dplyr::slice(which.max(p12)) %>% View()
coloc_tgfb1 %>% filter(minpval2 <= 5e-8, minpval1 <= 1e-4, p12 >= 0.8)


# HOTNET HEAT -------------------------------------------------------------

mega_example = mega_tbl %>% filter(gwas==my_gwas$gwas[1]) %>% dplyr::select(h4_log, hotnet_omnipath_heat, hgnc_symbol, entrezgene)
mega_example = mega_example[apply(mega_example, 1, function(x) !any(is.na(x))),]
save(mega_example, file="final_report/mega_example.RData") # savepoint

mega_example %>% ggplot(aes(h4_log, hotnet_omnipath_heat)) + geom_point(alpha=0.5, size=0.5) + theme_thesis(20) + xlab("H4") + ylab("Heat-diffused")

# look at the gene with the largest heat

my_genes = mega_example$entrezgene[mega_example$hotnet_omnipath_heat>4]
my_genes = "440689"

net_filt = graph.union(make_ego_graph(graph=net, order=2, nodes=as.character(my_genes)))
write_tsv(as.data.frame(get.edgelist(net_filt)), "tmp/omnipath_edges.txt")
mega_tbl %>% filter(gwas==my_gwas$gwas[1], entrezgene %in% V(net_filt)$name) %>% dplyr::select(h4_log, hotnet_omnipath_heat, hgnc_symbol, entrezgene) %>% write_tsv("tmp/mega_tbl_omnipath.txt")

# get all the hits

hits = mega_tbl %>% dplyr::select(h4_log, hotnet_omnipath_heat, success)
hits = arrange(hits, desc(hotnet_omnipath_heat))
hits = hits[apply(hits, 1, function(x) !any(is.na(x))),]
table(hits$success)

hits_long = hits %>% gather("method","score",1:2)

ggplot(hits_long, aes(m=score, d=success, color=method)) + geom_roc() + geom_rocci()


# ALEX'S CODE -------------------------------------------------------------

hits = mega_tbl %>% dplyr::select(h4_log, hotnet_omnipath_heat, success, pascal_default)
# hits = hits[apply(hits, 1, function(x) !any(is.na(x))),]
hits$pascal_default = -log(hits$pascal_default, base=10)
hits$random_score = rgamma(dim(hits)[1], 10)

hits_long = hits %>% gather("method","score",c(1:2,4:5))
hits_long$success = factor(hits_long$success, levels=c(1,0))
dim(hits_long) # 64,806,620
hits_long = hits_long[!is.na(hits_long$score),]
dim(hits_long) # 39,520,258
hits_long = hits_long[!is.infinite(hits_long$score),]
dim(hits_long) # 39,520,258
range(hits_long$score)

# does the success distribution make sense?
# i.e. are the more successes with higher genetic scores?

hits_dist = hits_long %>% group_by(method) %>% arrange(desc(score)) %>% group_modify(~ data.frame(hit_rank=which(.x$success==1)))
ggplot(hits_dist, aes(x=hit_rank)) + geom_histogram(bins=100) + facet_wrap(~method, scales="free") + theme_thesis(15)

hits_long_eval = hits_long %>% group_by(method) %>% arrange(score) %>% group_modify(~ pr_curve(.x, truth=success, score))
hits_long_eval = hits_long_eval %>% filter(!is.na(precision))

hits_long_eval_filt = hits_long_eval[seq(from=1, to=dim(hits_long_eval)[1], length=3e5),] # sample for faster plotting
save(hits_long_eval_filt, file="final_report/hits_long_eval_filt.RData") # savepoint

hits_long_eval_filt %>% filter(recall != 0, precision != 0) %>% ggplot(aes(x=recall, y=precision, color=method)) + geom_path(size=1) + xlab("Recall (% of known targets recovered)") + ylab("Precision (% of predictions that are known targets)") + scale_x_log10() + scale_y_log10() + scale_color_discrete(name="Gene Score", labels=c("GWAS","Propagated GWAS","Pascal","Random")) + theme_thesis(10)

hits_long %>% group_by(method) %>% arrange(score) %>% group_modify(~ pr_auc(.x, truth=success, score))


# DO GENETIC HITS CLUSTER PER METHOD / GWAS? ------------------------------

# pull out the network for each method

net = read_tsv("interactions/omnipath/omnipath_entrez.tsv")
net = graph_from_data_frame(net)
mean_distance(net)

all_dists = data.frame()
iters = 1e3

for(i in 1:length(my_gwas$gwas)) {
  
  print(i)
  
  # get the hits
  hits = as.numeric(unlist(mega_tbl[mega_tbl$gwas==my_gwas$gwas[i] & mega_tbl$high_confidence_genetic==1, 'entrezgene']))
  
  # store the results
  gwas_dist = data.frame(distance=rep(NA, iters), method="random")
  
  for(j in 1:length(gwas_dist$distance)) {
    sample_hits = sample(V(net), length(hits), replace=FALSE)
    sample_dists = distances(net, v=V(net)[V(net) %in% sample_hits], to=V(net)[V(net) %in% sample_hits])
    gwas_dist$distance[j] = mean(sample_dists[!is.infinite(sample_dists)])
  }
  
  hits_dist = distances(net, v=V(net)[V(net) %in% hits], to=V(net)[V(net) %in% hits])
  hits_dist = mean(hits_dist[!is.infinite(hits_dist)])
  
  gwas_dist = rbind(gwas_dist, data.frame(distance=hits_dist, method="hits"))
  gwas_dist$gwas = my_gwas$gwas[i]
  
  all_dists = rbind(all_dists, gwas_dist)
  
}

all_dists = tbl_df(all_dists)
all_dists_filt = all_dists %>% filter(gwas %in% my_gwas$gwas[2:20]) 

all_dists_filt %>% ggplot(aes(x=distance)) + stat_density(geom="line") + theme_thesis(10) + xlab("") + ylab("") + facet_wrap(~gwas) + geom_point(data=all_dists_filt[all_dists_filt$method=="hits",], aes(y=0), shape=17, color="red", size=2)


# CROSS PLOT - GENE SCORES VS OR ------------------------------------------

traits_109 = read_tsv("~/Downloads/gwas_list.txt")
mega_tbl_filt = mega_tbl %>% filter(gwas %in% traits_109$gwas)

methods_tested = c("high_confidence_genetic","complex","ligand_receptor","network_first_neighbor","pathways","pathway_first_neighbor","pathway_second_neighbor","random","hotnet_metabase_plus_ppi","hotnet_omnipath","network_second_neighbor","hotnet_huri","hotnet_string","hotnet_intomics")

res_or = mega_tbl_out(mega_tbl=mega_tbl_filt, col_names=methods_tested, include_seed=TRUE)[[3]]
gene_methods = c("pascal_default","pascal_100kb_up","magma_50kb")
my_measures = c("median","lower_quartile","upper_quartile")
col_names = paste(rep(gene_methods, each=length(my_measures)), my_measures, sep="_")
res_or = tbl_df(cbind(res_or, matrix(0, nrow=dim(res_or)[1], ncol=length(col_names))))
names(res_or)[6:14] = col_names

for(i in 1:length(gene_methods)) {
  for(j in 1:length(res_or$method)) {
    
    col_ix = which(str_replace(names(res_or),"_median","")==gene_methods[i])
    method_scores = mega_tbl_filt[mega_tbl_filt[,which(names(mega_tbl_filt)==res_or$method[j])]==1,which(names(mega_tbl_filt)==gene_methods[i])] %>% unlist() %>% as.numeric()
    q_scores = quantile(method_scores[!is.na(method_scores)])
    res_or[j,col_ix:(col_ix+2)] = q_scores[c(3,2,4)]
    
  }
}

i=1
for(i in 1:length(gene_methods)) {
  my_plot = res_or %>% ggplot(aes_string(x="odds", y=paste0(gene_methods[i],"_median"))) + geom_point(size=3.5, aes(color=method)) + theme_thesis(10) + xlab("Odds Ratio") + ylab("Median Gene Score") + ggtitle(gene_methods[i]) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.02, color="gray50") + geom_errorbar(aes_string(ymax=paste0(gene_methods[i],"_upper_quartile"), ymin=paste0(gene_methods[i],"_lower_quartile")), size=0.5, width=0.02, color="gray50")
  print(my_plot)
}

save(res_or, file="final_report/res_or.RData") # savepoint


# LIST FOR NIKOLINA -------------------------------------------------------

gene_buckets_no_seed = vector("list", length(methods_tested))
names(gene_buckets_no_seed) = methods_tested

for(i in 1:length(gene_buckets_no_seed)) {
  
  print(i)
  per_gwas = vector("list", length(traits_109$gwas))
  names(per_gwas) = traits_109$gwas
  
  for(j in 1:length(traits_109$gwas)) {
    per_gwas[[j]] = mega_tbl %>% filter(gwas==traits_109$gwas[j], !!sym(methods_tested[i])==1, high_confidence_genetic==0) %>% dplyr::select(entrezgene) %>% unlist %>% as.numeric() %>% unique()
  }
  
  gene_buckets_no_seed[[i]] = per_gwas
  
}


# PLOTS FOR PRESENTATION --------------------------------------------------

# 1. genetic hits, lr pairs, complex, ppi 1/2 neighbors, pathway, pathway 1 neighbors, random

load("final_report/res_plot_with_hcg_with_0s.RData")
load("final_report/method_avg.RData")

res_plot_1 = res_plot_cmh[c(1,3,4,5,7,11,13,14),]
method_avg_1 = method_avg[c(1,2,3,4,5,6,8,11),]
level_order = res_plot_1$method

plot_a = ggplot(res_plot_1, aes(y=factor(method, levels=level_order), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, color="red", alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

plot_b = method_avg_1 %>% ggplot(aes(x=factor(name,levels=level_order), y=log(Hits, base=10))) + geom_bar(stat="identity") + theme_thesis(15) + coord_flip() + xlab("") + ylab("Average Number of Hits (Log10) per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size=rel(0.8)))
plot_b

pdf(file="~/Downloads/plot_1.pdf", width=12, height=6)
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))
dev.off()

# plot 2

load("final_report/res_plot_without_hcg_with_0s.RData")
load("final_report/method_avg_without_hcg.RData")
res_plot_2 = res_plot_cmh[c(1,3,4,5,7,11,13),]
method_avg_2 = method_avg[c(2,3,4,5,6,8,11),]
level_order = res_plot_2$method

plot_a = ggplot(res_plot_2, aes(y=factor(method, levels=level_order), x=odds)) + geom_point(size=3.5, color="orange") + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, color="red", alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

plot_b = method_avg_2 %>% ggplot(aes(x=factor(name,levels=level_order), y=log(Hits, base=10))) + geom_bar(stat="identity") + theme_thesis(15) + coord_flip() + xlab("") + ylab("Average Number of Hits (Log10) per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size=rel(0.8)))
plot_b

pdf(file="~/Downloads/plot_2.pdf", width=12, height=6)
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))
dev.off()


# COMPARE GENE SCORES -----------------------------------------------------

y = mega_tbl %>% filter(gwas %in% my_gwas$gwas[1:10])
ggplot(y, aes(x=h4_log)) + geom_density() + facet_wrap(~gwas)
ggplot(y, aes(x=-log(pascal_default, base=10))) + geom_density() + facet_wrap(~gwas)


# FIRST AND SECOND NEIGHBOR PASCAL SCORES ---------------------------------

pascal_neighbors = mega_tbl_filt %>% select(entrezgene,gwas,mesh,disease,network_first_neighbor,network_second_neighbor,pascal_default,high_confidence_genetic)

to_plot = rbind(
  data.frame(
    score = pascal_neighbors %>% filter(network_first_neighbor==1) %>% select(pascal_default),
    group = "first_neighbors"
  ),
  data.frame(
    score = pascal_neighbors %>% filter(network_first_neighbor==0, network_second_neighbor==0) %>% select(pascal_default),
    group = "background"
  ),
  data.frame(
    score = pascal_neighbors %>% filter(network_first_neighbor==0, network_second_neighbor==1) %>% select(pascal_default),
    group = "second_neighbors"
  )
)

to_plot = to_plot[sample(1:dim(to_plot)[1], 1e4),]
to_plot %>% ggplot(aes(x=group, y=pascal_default)) + geom_boxplot()


# LAYERED PLOT WITH/WITHOUT PASCAL CUT-OFF --------------------------------

load("final_report/res_plot_without_hcg_with_0s.RData")
load("final_report/method_avg_without_hcg.RData")
res_plot_1 = res_plot_cmh
method_avg_1 = method_avg
rm(res_plot_cmh)
rm(method_avg)

res_plot_1$analysis = "total"
method_avg_1$analysis = "total"

load("final_report/res_plot_without_hcg_with_0s_pascal_filtered.RData")
load("final_report/method_avg_without_hcg_pascal_filtered.RData")
res_plot_2 = res_plot_cmh
method_avg_2 = method_avg
rm(res_plot_cmh)
rm(method_avg)

res_plot_2$analysis = "pascal cut-off"
method_avg_2$analysis = "pascal cut-off"

res_plot = rbind(res_plot_1, res_plot_2)
method_avg = rbind(method_avg_1, method_avg_2)

rm(res_plot_1); rm(res_plot_2); rm(method_avg_1); rm(method_avg_2); 

# subset of methods to plot for this plot

methods_to_use = c(
  "high_confidence_genetic",
  "complex",
  "ligand_receptor",
  "string_first_neighbors",
  "string_second_neighbors",
  "pathways",
  "pathway_first_neighbor",
  "pathway_second_neighbor",
  "random"
)

# res_plot_basic_methods = res_plot %>% filter(method %in% methods_to_use)
res_plot_basic_methods = res_plot
res_plot_basic_methods$method = as.character(res_plot_basic_methods$method)

# method_avg_basic_methods = method_avg %>% filter(method %in% methods_to_use)
method_avg_basic_methods = method_avg

res_plot_basic_methods = res_plot_basic_methods %>% arrange(method, analysis)
method_avg_basic_methods = method_avg_basic_methods %>% arrange(method, analysis)

all(as.character(res_plot_basic_methods$method) == method_avg_basic_methods$method)
all_basic_methods = full_join(res_plot_basic_methods, method_avg_basic_methods)

# clean up labels for plotting

plot_a = ggplot(all_basic_methods, aes(y=factor(method), x=odds, group=analysis, color=analysis)) + geom_point(size=2,  position=ggstance::position_dodgev(height=0.5)) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50", position=ggstance::position_dodgev(height=0.5)) + theme_thesis(10) + xlab("") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

level_order = levels(factor(all_basic_methods$method))

plot_b = ggplot(all_basic_methods, aes(x=factor(method, levels=level_order), y=Hits, fill=analysis)) + geom_errorbar(aes(ymax=Hits+SD, ymin=2), size=0.75, width=0.2, position=position_dodge(width=1)) + geom_bar(stat="identity", position="dodge") + theme_thesis(10) + scale_y_log10() + coord_flip() + xlab("") + ylab("Average Number of Hits per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_text(size=rel(0.8)))
plot_b

plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# OVERLAP ACROSS NETWORKS -------------------------------------------------

# what is the intersection across the different networks

networks = c("omnipath","huri","string","intomics","metabase_plus_ppi")
net_list = vector("list", length(networks))
names(net_list) = networks

for(i in 1:length(networks)) {
  
  print(networks[i])
  if(networks[i]=="string") next
  
  net = read_tsv(paste0("interactions/", networks[i], "/", networks[i], "_entrez.txt"))
  if(i==4) {
    net_edges = paste(unlist(net[,7]), unlist(net[,8]), sep="_")
  } else {
    net_edges = paste(unlist(net[,1]), unlist(net[,2]), sep="_")
  }
  
  net_edges = unique(net_edges)
  net_list[[i]] = net_edges
  
}

require(eulerr)
plot(euler(net_list))
plot(euler(net_list[-3]))


# TOTAL OR PLOT -----------------------------------------------------------

load("final_report/res_plot_without_hcg_with_0s.RData")
load("final_report/res_plot_pascal_without_hcg_with_0s.RData")
load("final_report/method_avg_without_hcg.RData")
load("final_report/method_avg_pascal_without_hcg.RData")

method_avg = rbind(method_avg, method_avg_pascal)
res_plot_cmh = rbind(res_plot_cmh, res_plot_cmh_pascal)

res_plot_cmh = res_plot_cmh %>% arrange(odds)
method_avg = method_avg[match(as.character(res_plot_cmh$method), method_avg$method),]

# remove 'network n neighbors' and cell-type specific pascal
res_plot_cmh = res_plot_cmh %>% filter(!method %in% c("network_first_neighbor","network_second_neighbor","pascal.mph.sig.genes.clean","pascal.cd4.sig.genes.clean"))
method_avg = method_avg %>% filter(!method %in% c("network_first_neighbor","network_second_neighbor","pascal.mph.sig.genes.clean","pascal.cd4.sig.genes.clean"))

map_names = data.frame(
  orig = c("string_second_neighbors","pathways","huri_second_neighbors","pathway_second_neighbor","metabase_plus_ppi_first_neighbors","network_first_neighbor","metabase_plus_ppi_second_neighbors","omnipath_second_neighbors","omnipath_first_neighbors","string_first_neighbors","random","hotnet_huri","network_second_neighbor","intomics_second_neighbors","intomics_first_neighbors","huri_first_neighbors","hotnet_string","pathway_first_neighbor","pascal.reactome.sig.genes.clean","hotnet_omnipath","hotnet_metabase_plus_ppi","pascal.metabase.sig.genes.clean","hotnet_intomics","pascal.ppi.sig.genes.clean","ligand_receptor","pascal.coexp.sig.genes.clean","pascal.cd4.sig.genes.clean","complex","high_confidence_genetic","pascal.mph.sig.genes.clean","pascal_cell_type_combined"),
  map = c("STRING 1st/2nd Neighbor","All Pathway","HuRI 1st/2nd Neighbor","Pathway 1st/2nd Neighbor","Metabase+ 1st Neighbor","Network 1st Neighbor","Metabase+ 1st/2nd Neighbor","OmniPath 1st/2nd Neighbor","Omnipath 1st Neighbor","STRING 1st Neighbor","Random","HotNet2 HuRI","Network 1st/2nd Neighbor","Intomics 1st/2nd Neighbor","Intomics 1st Neighbor","HuRI 1st Neighbor","HotNet2 STRING","Pathway 1st Neighbor","Pascal Reactome","HotNet2 OmniPath","HotNet2 Metabase+","Pascal Metabase","HotNet2 Intomics","Pascal PPI","Ligand Receptor","Pascal Co-Expression","Pascal CD4+","Complex","HCGH","Pascal Macrophage","Pascal Combined Cell-Type"),
  type = c("Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","HotNet2","Basic","Basic","Basic","Basic","HotNet2","Basic","Pascal","HotNet2","HotNet2","Pascal","HotNet2","Pascal","Basic","Pascal","Pascal","Basic","Basic","Pascal","Pascal")
)

res_plot_cmh$label = map_names$map[match(res_plot_cmh$method, map_names$orig)]
res_plot_cmh$`Method Type` = factor(map_names$type[match(res_plot_cmh$method, map_names$orig)])
method_avg$label = map_names$map[match(method_avg$method, map_names$orig)]

level_order = res_plot_cmh$label
plot_a = ggplot(res_plot_cmh, aes(y=factor(label, levels=level_order), x=odds)) + geom_point(aes(color=`Method Type`), size=3.5) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(25) + xlab("") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

plot_b = method_avg %>% ggplot(aes(x=factor(label,levels=level_order), y=log(Hits, base=10))) + geom_bar(stat="identity") + theme_thesis(20) + coord_flip() + xlab("") + ylab("Average Number of Additional Hits (Log10)") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))

