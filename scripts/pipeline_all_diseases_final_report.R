
require(tidyverse)
require(epiChoose)
require(cowplot)


# PREP DATA ---------------------------------------------------------------

load(".RData")

# mega_tbl = read_tsv("mega_tbl.txt")

load("final_report/res_plot_without_hcg_with_0s_095.RData")
load("final_report/method_avg_without_hcg_095.RData")

load("final_report/res_plot_pascal_without_hcg_with_0s.RData")
load("final_report/method_avg_pascal_without_hcg.RData")

method_avg = rbind(method_avg, method_avg_pascal)
res_plot_cmh = rbind(res_plot_cmh, res_plot_cmh_pascal)

res_plot_cmh = res_plot_cmh %>% arrange(odds)
method_avg = method_avg[match(as.character(res_plot_cmh$method), method_avg$method),]

# remove 'network n neighbors' and cell-type specific pascal
res_plot_cmh = res_plot_cmh %>% filter(!method %in% c("network_first_neighbor","network_second_neighbor","pascal.mph.sig.genes.clean","pascal.cd4.sig.genes.clean"))
method_avg = method_avg %>% filter(!method %in% c("network_first_neighbor","network_second_neighbor","pascal.mph.sig.genes.clean","pascal.cd4.sig.genes.clean"))

map_names = data.frame(
  orig = c("string_second_neighbors","pathways","huri_second_neighbors","pathway_second_neighbor","metabase_plus_ppi_first_neighbors","network_first_neighbor","metabase_plus_ppi_second_neighbors","omnipath_second_neighbors","omnipath_first_neighbors","string_first_neighbors","random","hotnet_huri","network_second_neighbor","intomics_second_neighbors","intomics_first_neighbors","huri_first_neighbors","hotnet_string","pathway_first_neighbor","pascal.reactome.sig.genes.clean","hotnet_omnipath","hotnet_metabase_plus_ppi","pascal.metabase.sig.genes.clean","hotnet_intomics","pascal.ppi.sig.genes.clean","ligand_receptor","pascal.coexp.sig.genes.clean","pascal.cd4.sig.genes.clean","complex","high_confidence_genetic","pascal.mph.sig.genes.clean","pascal_cell_type_combined","magma.metabase.sig.genes.clean","magma.ppi.sig.genes.clean","magma.reactome.sig.genes.clean","magma.coexp.sig.genes.clean"),
  map = c("STRING 1st/2nd Neighbor","All Pathway","HuRI 1st/2nd Neighbor","Pathway 1st/2nd Neighbor","Metabase+ 1st Neighbor","Network 1st Neighbor","Metabase+ 1st/2nd Neighbor","OmniPath 1st/2nd Neighbor","Omnipath 1st Neighbor","STRING 1st Neighbor","Random","HotNet2 HuRI","Network 1st/2nd Neighbor","InBio Map 1st/2nd Neighbor","InBio Map 1st Neighbor","HuRI 1st Neighbor","HotNet2 STRING","Pathway 1st Neighbor","Pascal Reactome","HotNet2 OmniPath","HotNet2 Metabase+","Pascal Metabase","HotNet2 InBio Map","Pascal PPI","Ligand Receptor","Pascal Co-Expression","Pascal CD4+","Complex","HCGH","Pascal Macrophage","Pascal Combined Cell-Type","Magma Metabase","Magma PPI","Magma Reactome","Magma Co-Expression"),
  type = c("Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","Basic","HotNet2","Basic","Basic","Basic","Basic","HotNet2","Basic","Pascal","HotNet2","HotNet2","Pascal","HotNet2","Pascal","Basic","Pascal","Pascal","Basic","Basic","Pascal","Pascal","Magma","Magma","Magma","Magma")
)

rm(res_plot_cmh_pascal)
rm(method_avg_pascal)


# SUPPLEMENTARY FIGURE 3 --------------------------------------------------

# index the relevant rows
ix_full = c(5,6,7,15,19,22,23,24,25,30,32,35)
res_plot_cmh_full = res_plot_cmh[-ix_full,]
method_avg_full = method_avg[-ix_full,]

res_plot_cmh_full$label = map_names$map[match(res_plot_cmh_full$method, map_names$orig)]
res_plot_cmh_full$`Method Type` = factor(map_names$type[match(res_plot_cmh_full$method, map_names$orig)])
method_avg_full$label = map_names$map[match(method_avg_full$method, map_names$orig)]

level_order = res_plot_cmh_full$label
plot_a = ggplot(res_plot_cmh_full, aes(y=factor(label, levels=level_order), x=odds)) + geom_point(aes(color=`Method Type`), size=3.5) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

plot_b = method_avg_full %>% ggplot(aes(x=factor(label,levels=level_order), y=Hits, base=10)) + geom_bar(stat="identity") + theme_thesis(15) + scale_y_log10() + coord_flip() + xlab("") + ylab("Average Number of Additional Hits per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# FIGURE 2 ----------------------------------------------------------------

# index the relevant rows
ix_1 = c(1,2,4,10,11,18,29,33,34)
res_plot_cmh_1 = res_plot_cmh[ix_1,]
method_avg_1 = method_avg[ix_1,]

res_plot_cmh_1$label = map_names$map[match(res_plot_cmh_1$method, map_names$orig)]
res_plot_cmh_1$`Method Type` = factor(map_names$type[match(res_plot_cmh_1$method, map_names$orig)])
method_avg_1$label = map_names$map[match(method_avg_1$method, map_names$orig)]
res_plot_cmh_1$`Method Type` = c("Network","Pathway","Pathway","Network","Random","Pathway","HC Interactions","HC Interactions","HCGH")

level_order = res_plot_cmh_1$label
plot_a = ggplot(res_plot_cmh_1, aes(y=factor(label, levels=level_order), x=odds)) + geom_point(aes(color=`Method Type`), size=3.5) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) + theme(legend.position="none") # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

method_avg_1$`Method Type` = res_plot_cmh_1$`Method Type`[match(method_avg_1$method, res_plot_cmh_1$method)]
plot_b = method_avg_1 %>% ggplot(aes(x=factor(label,levels=level_order), y=Hits, base=10)) + geom_bar(stat="identity", aes(fill=`Method Type`)) + theme_thesis(15) + scale_y_log10() + coord_flip() + xlab("") + ylab("Average Number of Additional Hits per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# SUPPLEMENTARY FIGURE 1 --------------------------------------------------

# naive all other networks (first/second neighbors)

# index the relevant rows
ix_2 = c(3,8,9,13,14,16)
res_plot_cmh_2 = res_plot_cmh[ix_2,]
method_avg_2 = method_avg[ix_2,]

res_plot_cmh_2$label = map_names$map[match(res_plot_cmh_2$method, map_names$orig)]
res_plot_cmh_2$`Method Type` = factor(map_names$type[match(res_plot_cmh_2$method, map_names$orig)])
method_avg_2$label = map_names$map[match(method_avg_2$method, map_names$orig)]

level_order = res_plot_cmh_2$label
plot_a = ggplot(res_plot_cmh_2, aes(y=factor(label, levels=level_order), x=odds)) + geom_point(color="gray50", size=3.5) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

plot_b = method_avg_2 %>% ggplot(aes(x=factor(label,levels=level_order), y=Hits, base=10)) + geom_bar(stat="identity") + theme_thesis(15) + scale_y_log10() + coord_flip() + xlab("") + ylab("Average Number of Additional Hits per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# FIGURE 3 ----------------------------------------------------------------

# advanced methods - hotnet and pascal

# index the relevant rows
ix_3 = c(9,10,12,14,16,17,20,21,26,27,28,31,34)
res_plot_cmh_3 = res_plot_cmh[ix_3,]
method_avg_3 = method_avg[ix_3,]

res_plot_cmh_3$label = map_names$map[match(res_plot_cmh_3$method, map_names$orig)]
res_plot_cmh_3$`Method Type` = factor(map_names$type[match(res_plot_cmh_3$method, map_names$orig)])
method_avg_3$label = map_names$map[match(method_avg_3$method, map_names$orig)]
res_plot_cmh_3$`Method Type` = as.character(res_plot_cmh_3$`Method Type`)
res_plot_cmh_3$`Method Type`[res_plot_cmh_3$`Method Type`=="Basic"] = "First Neighbor"
res_plot_cmh_3$`Method Type`[res_plot_cmh_3$method=="high_confidence_genetic"] = "HCGH"
res_plot_cmh_3$Network = factor(c("OmniPath","STRING","HuRI","InBio Map","HuRI","STRING","Pascal Reactome","OmniPath","Pascal Metabase","InBio Map","Pascal PPI","Pascal Co-Expression","HCGH"))
method_avg_3$Network = res_plot_cmh_3$Network
method_avg_3$`Method Type` = res_plot_cmh_3$`Method Type`

level_order = rev(c("HCGH","Pascal Co-Expression","Pascal PPI","Pascal Metabase","Pascal Reactome","InBio Map","OmniPath","HuRI","STRING"))
plot_a = ggplot(res_plot_cmh_3, aes(y=factor(Network, levels=level_order), x=odds, group=`Method Type`, color=`Method Type`)) + geom_point(size=2, position=ggstance::position_dodgev(height=0.5)) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.75, height=0.2, position=ggstance::position_dodgev(height=0.5)) + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2)+ theme(legend.position="none") # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

plot_b = ggplot(method_avg_3, aes(x=factor(Network, levels=level_order), y=Hits, fill=`Method Type`)) + geom_errorbar(aes(ymax=Hits+SD, ymin=2), size=0.75, width=0.2, position=position_dodge(width=1)) + geom_bar(stat="identity", position= "dodge") + theme_thesis(15) + scale_y_log10() + coord_flip() + xlab("") + ylab("Average Number of Hits per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size=rel(0.8)))
plot_b

plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# FIGURE 4 ----------------------------------------------------------------

mega_tbl = read_tsv("mega_tbl.txt")
include_ix = which(mega_tbl$gwas %in% my_gwas_filt$gwas)
pascal_vec = mega_tbl$pascal_default[include_ix]
rm(mega_tbl)

# need to pull the pascal scores from the different tables

to_plot_all = data.frame()
tables = c("mega_tbl.txt","mega_tbl_first_neighbors.txt","mega_tbl_second_neighbors.txt","mega_tbl_gb.txt")

for(j in 1:length(tables)) {

  print(j)
    
  pascal_dat = data.frame()
  mt = read_tsv(tables[j])
  if(j!=4) mt = mt[include_ix,] # the pascal mega table is already filtered
  
  if(j==2) names(mt)[(dim(mt)[2]-4):dim(mt)[2]] = paste(names(mt)[(dim(mt)[2]-4):dim(mt)[2]], "first_neighbors", sep="_")
  if(j==3) names(mt)[(dim(mt)[2]-4):dim(mt)[2]] = paste(names(mt)[(dim(mt)[2]-4):dim(mt)[2]], "second_neighbors", sep="_")
  
  for(i in 1:length(map_names$orig)) {
    
    if(map_names$orig[i] %in% names(mt)) {
      print(map_names$orig[i])
      to_add = data.frame(
        method = map_names$orig[i],
        pascal = pascal_vec[mt[,which(names(mt)==map_names$orig[i])]==1]
      )
      pascal_dat = rbind(pascal_dat, to_add)
    }
  }
 
  pascal_dat$pascal_log = -log(pascal_dat$pascal, base=10)
  to_plot = pascal_dat %>% group_by(method) %>% summarise(median=quantile(pascal_log, probs=0.50, na.rm=TRUE), lower_q=quantile(pascal_log, probs=0.25, na.rm=TRUE), upper_q=quantile(pascal_log, probs=0.75, na.rm=TRUE), min_q=0)
  to_plot$max_q = to_plot$median + 1.5*(to_plot$upper_q - to_plot$lower_q)
  to_plot$label = map_names$map[match(to_plot$method, map_names$orig)]
  to_plot_all = rbind(to_plot_all, to_plot)
  
  rm(to_plot)
  rm(mt)
  rm(pascal_dat)
   
}

level_order_a = res_plot_cmh_1$label
to_plot = to_plot_all[to_plot_all$label %in% level_order_a,]
to_plot = distinct(to_plot)
plot_a = ggplot(to_plot, aes(x=factor(label, levels=level_order_a), middle=median, lower=lower_q, upper=upper_q, ymin=min_q, ymax=max_q)) + geom_boxplot(stat="identity", width=0.5) + geom_hline(yintercept=to_plot_all$median[to_plot_all$method=="random"], alpha=0.5, linetype=2) + theme_thesis(15) + ylab("-log10(GeneScore)") + xlab("") + coord_flip() + scale_y_sqrt()

level_order_b = c("HotNet2 HuRI","HotNet2 STRING","HotNet2 OmniPath","HotNet2 Intomics","Pascal Reactome","Pascal Metabase","Pascal PPI","Pascal Co-Expression","HCGH")
to_plot = to_plot_all[to_plot_all$label %in% level_order_b,]
to_plot = distinct(to_plot)
plot_b = ggplot(to_plot, aes(x=factor(label, levels=level_order_b), middle=median, lower=lower_q, upper=upper_q, ymin=min_q, ymax=max_q)) + geom_boxplot(stat="identity", width=0.5) + geom_hline(yintercept=to_plot_all$median[to_plot_all$method=="random"], alpha=0.5, linetype=2) + theme_thesis(15) + ylab("-log10(GeneScore)") + xlab("") + coord_flip() + scale_y_sqrt()

plot_grid(plot_a, plot_b, align="h", rel_widths=c(1,1))


# SUPPLEMENTARY FIGURE 5 --------------------------------------------------

# pascal scores of first/second neighbors across networks

mega_tbl = read_tsv("mega_tbl.txt")
include_ix = which(mega_tbl$gwas %in% my_gwas_filt$gwas)
mega_tbl = mega_tbl[include_ix,]
mega_tbl = mega_tbl %>% select(entrezgene,gwas,mesh,disease,pascal_default,high_confidence_genetic)

mega_tbl_first_neighbors = read_tsv("mega_tbl_first_neighbors.txt")
names(mega_tbl_first_neighbors)[12:15] = paste(names(mega_tbl_first_neighbors)[12:15], "first_neighbors", sep="_")
mega_tbl_first_neighbors = mega_tbl_first_neighbors[include_ix,]
mega_tbl = cbind(mega_tbl, mega_tbl_first_neighbors[,12:15])
rm(mega_tbl_first_neighbors)

mega_tbl_second_neighbors = read_tsv("mega_tbl_second_neighbors.txt")
names(mega_tbl_second_neighbors)[12:15] = paste(names(mega_tbl_second_neighbors)[12:15], "second_neighbors", sep="_")
mega_tbl_second_neighbors = mega_tbl_second_neighbors[include_ix,]
mega_tbl = cbind(mega_tbl, mega_tbl_second_neighbors[,12:15])
rm(mega_tbl_second_neighbors)

mega_tbl = tbl_df(mega_tbl)

networks = c("omnipath","huri","string","intomics")
to_plot_all = data.frame()

for(i in 1:length(networks)) {
  
  print(networks[i])
  
  to_plot = rbind(
    data.frame(
      score = mega_tbl %>% filter(!!sym(paste0(networks[i],"_first_neighbors"))==1) %>% select(pascal_default),
      group = "first_neighbors",
      network = networks[i]
    ),
    data.frame(
      score = mega_tbl %>% filter(!!sym(paste0(networks[i],"_first_neighbors"))==0, !!sym(paste0(networks[i],"_second_neighbors"))==0) %>% select(pascal_default),
      group = "background",
      network = networks[i]
    ),
    data.frame(
      score = mega_tbl %>% filter(!!sym(paste0(networks[i],"_first_neighbors"))==0, !!sym(paste0(networks[i],"_second_neighbors"))==1) %>% select(pascal_default),
      group = "second_neighbors",
      network = networks[i]
    )
  )
  
  to_plot_all = rbind(to_plot_all, to_plot)
  
}

to_plot_pascal_n  = to_plot_all %>% group_by(network) %>% do(sample_n(.,1e4))

to_plot_pascal_n$network = as.character(to_plot_pascal_n$network)
to_plot_pascal_n$network[to_plot_pascal_n$network=="omnipath"] = "OmniPath"
to_plot_pascal_n$network[to_plot_pascal_n$network=="string"] = "STRING"
to_plot_pascal_n$network[to_plot_pascal_n$network=="intomics"] = "InBio Map"
to_plot_pascal_n$network[to_plot_pascal_n$network=="huri"] = "HuRI"

to_plot_pascal_n %>% ggplot(aes(x=group, y=-log(pascal_default, base=10))) + geom_boxplot() + facet_wrap(~network) + xlab("") + ylab("-log10(GeneScore)") + coord_flip() + theme_thesis(20)
save(to_plot_pascal_n, file="r_data/to_plot_pascal_n.RData")


# SUPPLEMENTARY FIGURE 4 --------------------------------------------------

# what is the intersection across the different networks

networks = c("omnipath","huri","string","intomics")
net_list = vector("list", length(networks))
names(net_list) = networks

for(i in 1:length(networks)) {
  
  print(networks[i])
  # if(networks[i]=="string") next
  
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
names(net_list) = c("OmniPath","HuRI","STRING","InBio Map")
plot(euler(net_list[-3]), quantities=TRUE)
plot(euler(net_list))


# SUPPLEMENTARY FIGURE 6 --------------------------------------------------

load("gwas_650_sigGenesets_sigGenes.RData")
load("gwas_650_sigGenesets_sigGenes_magma.RData")
mega_tbl_gb = read_tsv("mega_tbl_gb.txt")

gb_list = list(
  magma.coexp.sig.genes.clean = magma.coexp.sig.genes.clean,
  magma.metabase.sig.genes.clean = magma.metabase.sig.genes.clean,
  magma.ppi.sig.genes.clean = magma.ppi.sig.genes.clean,
  magma.reactome.sig.genes.clean = magma.reactome.sig.genes.clean
)

res_magma = mega_tbl_out(mega_tbl_gb, col_names=names(gb_list), include_seed=FALSE)
res_plot_cmh_magma = res_magma$res_plot_cmh

method_avg_magma = method_delta(mega_tbl=mega_tbl_gb, col_names=names(gb_list), my_gwas=my_gwas_filt)
method_avg_magma = method_avg_magma[match(res_plot_cmh_magma$method, method_avg_magma$method),]

res_plot_cmh = rbind(res_plot_cmh, res_plot_cmh_magma)
method_avg = rbind(method_avg, method_avg_magma)

# index the relevant rows
ix_4 = c(24,22,25,30,20,28,26,31,34)
res_plot_cmh_4 = res_plot_cmh[ix_4,]
method_avg_4 = method_avg[ix_4,]

res_plot_cmh_4$label = map_names$map[match(res_plot_cmh_4$method, map_names$orig)]
res_plot_cmh_4$`Method Type` = factor(map_names$type[match(res_plot_cmh_4$method, map_names$orig)])
method_avg_4$label = map_names$map[match(method_avg_4$method, map_names$orig)]

level_order = res_plot_cmh_4$label
plot_a = ggplot(res_plot_cmh_4, aes(y=factor(label, levels=level_order), x=odds)) + geom_point(aes(color=`Method Type`), size=3.5) + geom_errorbarh(aes(xmax=ci_high, xmin=ci_low), size=0.5, height=0.2, color="gray50") + theme_thesis(15) + xlab("Odds Ratio") + ylab("") + geom_vline(xintercept=1, alpha=0.5, linetype=2) + theme(legend.position="none") # + coord_cartesian(xlim=c(0, ci_high+1))
plot_a

method_avg_4$`Method Type` = res_plot_cmh_4$`Method Type`[match(method_avg_4$method, res_plot_cmh_4$method)]
plot_b = method_avg_4 %>% ggplot(aes(x=factor(label,levels=level_order), y=Hits, base=10)) + geom_bar(stat="identity", aes(fill=`Method Type`)) + theme_thesis(15) + scale_y_log10() + coord_flip() + xlab("") + ylab("Average Number of Additional Hits per GWAS") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
plot_grid(plot_a, plot_b, align="h", rel_widths=c(1.5,1))


# MODULES PREP ------------------------------------------------------------

load(".RData")
mega_tbl = read_tsv("mega_tbl.txt")

mega_tbl_modules = mega_tbl %>% dplyr::select(1:10)

networks = c("metabase_plus_ppi","string","omnipath","huri","intomics")
col_names = paste0(networks, "_modules")
mega_tbl_modules = tbl_df(cbind(mega_tbl_modules, matrix(0, nrow=dim(mega_tbl_modules)[1], ncol=length(col_names))))
names(mega_tbl_modules)[11:15] = col_names

for(i in 1:length(col_names)) {
  
  print(paste("Network:", col_names[i]))
  mega_tbl_ix = which(names(mega_tbl_modules)==col_names[i])
  
  for(j in 1:length(my_gwas$mesh)) {
    
    if(file.exists(paste0("hotnet/hotnet/test_all_gwas/results/", networks[i], "/output_", j, "/consensus/subnetworks.tsv"))) {
      res = read_tsv(paste0("hotnet/hotnet/test_all_gwas/results/", networks[i], "/output_", j, "/consensus/subnetworks.tsv"), col_names=FALSE)
      res = str_split(res$X1[3:length(res$X1)], " ")
      for(k in 1:length(res)) {
        if(networks[i]=="metabase_plus_ppi") { # use gene symbols
          mega_tbl_modules[(mega_tbl_modules$gwas==my_gwas$gwas[j] & mega_tbl_modules$hgnc_symbol %in% res[[k]]), mega_tbl_ix] = k
        } else {
          mega_tbl_modules[(mega_tbl_modules$gwas==my_gwas$gwas[j] & mega_tbl_modules$entrezgene %in% as.numeric(res[[k]])), mega_tbl_ix] = k
        }
      }
    } else {
      mega_tbl_modules[mega_tbl_modules$gwas==my_gwas$gwas[j], mega_tbl_ix] = NA
    }
    
  }
}


# MODULES ANALYSIS --------------------------------------------------------

mega_tbl_modules = read_tsv("mega_tbl_modules.txt", col_types=cols(intomics_modules=col_double()))

module_dat = mega_tbl_modules %>% select(hgnc_symbol, gwas, high_confidence_genetic, matches("_modules")) %>% gather("method","module", 4:8) %>% filter(module!=0 & !is.na(module)) %>% group_by(gwas, module, hgnc_symbol, method) %>% summarise(high_confidence_genetic=max(high_confidence_genetic))
module_dat$module = factor(module_dat$module)

i=2
p = module_dat %>% filter(gwas==my_gwas$gwas[i]) %>% ggplot(aes(x=forcats::fct_infreq(module))) + geom_bar() + facet_wrap(~method, scales="free") + theme_thesis(10)
p + geom_bar(aes(fill=factor(high_confidence_genetic))) + xlab("Module") + ylab("Size") # + coord_cartesian(xlim=c(0,20)) 

# pick out examples of modules where hgch >= 2 and drug successes >=1
network = "omnipath"
moi = mega_tbl_modules %>%
  select(high_confidence_genetic,gwas,success,!!sym(paste0(network,"_modules"))) %>%
  filter(success!=high_confidence_genetic) %>%
  group_by(!!sym(paste0(network,"_modules")), gwas) %>%
  summarise(hcghs=sum(high_confidence_genetic), drugs=sum(success), module_size=n()) %>%
  filter(!!sym(paste0(network,"_modules"))!=0, hcghs>=2, drugs>=1) %>%
  arrange(desc(drugs),desc(hcghs))

# work through example

net = read_tsv(paste0("interactions/", network, "/", network, "_entrez.txt"))
net = graph_from_data_frame(net[,1:2]) # get the network into igraph

which_gwas = "GSK500KV3lin_HES_Disorders_of_lipoprotein_metabolism_lipidaemias" # pick the gwas of interest
which_module = 1 # pick hotnet2 module

# pick example
module = mega_tbl_modules %>% filter(gwas==which_gwas, !!sym(paste0(network,"_modules"))==which_module)
net_filt = subgraph(net, v=as.character(module$entrezgene))
write_tsv(as.data.frame(get.edgelist(net_filt)), "tmp/module_edges.txt")
write_tsv(module, "tmp/module_nodes.txt")

# what are the enriched pascal pathways for which_gwas

load("gwas_650_sigGenesets_sigGenes.RData")
load("gwas_650_sigGenesets.RData")

# significant genes
pascal_entrez = pascal.metabase.sig.genes.clean[[which(names(pascal.metabase.sig.genes.clean)==which_gwas)]]

# significant pathways
pascal_pathways = pascal.metabase.sig[[which(names(pascal.metabase.sig)==which_gwas)]]

# all metabase pathways
metabase_pathways = gmtPathways("genesets/Metabase_GO_Maps_entrez.filt10.gmt")

# top pascal pathway nodes
my_p = as.numeric(metabase_pathways[[which(names(metabase_pathways)==pascal_pathways$Name[1])]])
plot(euler(list(pascal=my_p, hotnet=unique(module$entrezgene))), quantities=TRUE)

module_pascal = mega_tbl_modules %>% filter(gwas==which_gwas, entrezgene %in% my_p)
mb = read_tsv("interactions/metabase_plus_ppi/metabase/MetaBaseFullNetworkHML.txt")
mb_filt = mb %>% filter(TrustLevel=="high", !Mechanism %in% c("Competition","Influence on expression","Pharmacological effect", "Toxic effect","Unspecified"))
mb_graph = graph_from_data_frame(mb_filt[,c(1,3,5:11)])
mb_filt = subgraph(mb_graph, v=as.character(my_p))
write_tsv(as.data.frame(get.edgelist(mb_filt)), "tmp/pascal_edges.txt")
write_tsv(module_pascal, "tmp/pascal_nodes.txt")
intersect(unique(module_pascal$hgnc_symbol), unique(module$hgnc_symbol))


# SUPPLEMENTARY FIGURE 2 --------------------------------------------------

# what are the 10 most "unsuccessful" targets?

d = 500
target_drug = mega_tbl_modules %>%
  select(success,failure,hgnc_symbol,entrezgene) %>%
  group_by(hgnc_symbol) %>%
  summarise(all_success=sum(success), all_failure=sum(failure), entrezgene=entrezgene[1]) %>%
  mutate(diff=all_failure-all_success) %>% arrange(desc(diff))

wins = target_drug %>% tail(d)
fails = target_drug %>% head(d)

network = "string"
net = read_tsv(paste0("interactions/", network, "/", network, "_entrez.txt"))
net = graph_from_data_frame(net[,1:2])

n_degree = degree(net)
fails_score = as.numeric(n_degree[V(net)$name %in% fails$entrezgene])
background_score = as.numeric(n_degree[!V(net)$name %in% fails$entrezgene])
to_plot = data.frame(
  group = factor(c(rep("Fails",length(fails_score)), rep("Background",length(background_score)))),
  degree = c(fails_score, background_score)
)

to_plot %>% ggplot(aes(group, log(degree))) + geom_boxplot() + theme_thesis(20) + ylab("Target Degree (Log2)") + xlab("Target Type")
t.test(log(fails_score), log(background_score))


# OVERLAP OF HOTNET2 MODULES AND NEIGHBOURS -------------------------------

# does hotnet2 draw its modules from first/second neighbors?

method_overlaps =   list(
  `high confidence genetics` = mega_tbl$index[mega_tbl$high_confidence_genetic==1],
  network_first_neighbor = mega_tbl$index[mega_tbl$network_first_neighbor==1],
  network_second_neighbor = mega_tbl$index[mega_tbl$network_second_neighbor==1],
  hotnet_omnipath = mega_tbl$index[mega_tbl$hotnet_omnipath==1]
)

my_dat_genes_all_tbl = as.data.frame.matrix((table(stack(method_overlaps))))
my_dat_genes_all_tbl = cbind(rownames(my_dat_genes_all_tbl), my_dat_genes_all_tbl)
rownames(my_dat_genes_all_tbl) = NULL
upset(my_dat_genes_all_tbl, order.by="freq")

