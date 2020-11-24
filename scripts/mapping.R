
require(biomaRt)
require(AnnotationDbi)
require(igraph)


# CONVERT ENTREZ IDS TO ENSEMBL -------------------------------------------

# import network/nodes
net = read_tsv("/Volumes/am673712/links/network_analysis/networks/metabase_plus_ppi/metabase_plus_ppi_edge_list.txt", col_names=FALSE)
genes = read_tsv("/Volumes/am673712/links/network_analysis/networks/metabase_plus_ppi/metabase_plus_ppi_gene_index.txt", col_names=FALSE)

require(biomaRt)
mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mapping_1 = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id"), filters="hgnc_symbol", values=genes$X2, mart=mart)
table(genes$X2 %in% mapping_1$hgnc_symbol)
net_mapped = data.frame(from=rep(NA,dim(net)[1]), to=rep(NA,dim(net)[1]))
net_mapped$from = mapping_1$hgnc_symbol[match(net$X1, mapping_1$entrezgene)]
net_mapped$to = mapping_1$hgnc_symbol[match(net$X2, mapping_1$entrezgene)]
net_mapped = net_mapped[apply(net_mapped, 1, function(x) !any(is.na(x))),]


# COMPLEX NETWORK ---------------------------------------------------------

complex = read_tsv("interactions/complex_portal/homo_sapiens.tsv")
x = complex$`Identifiers (and stoichiometry) of molecules in complex`[1]
c_ids = sapply(complex$`Identifiers (and stoichiometry) of molecules in complex`, function(x) str_extract_all(x, "[A-Z][[:alnum:]]{5}"))
uniprot = Rkeys(org.Hs.egUNIPROT)
mapping_2 = select(org.Hs.eg.db, uniprot, "ENTREZID", "UNIPROT")

net_complex = data.frame(from=NULL, to=NULL, complex=NULL)
for(i in 1:length(c_ids)) {
  c_entrez = mapping_2$ENTREZID[match(c_ids[[i]], mapping_2$UNIPROT)]
  if(length(c_entrez)>1) {
    to_add = data.frame(t(combn(c_entrez, 2)), complex$`Recommended name`[i])
    names(to_add) = c("from","to","complex")
    net_complex = rbind(net_complex, to_add)
  }
}

write_tsv(net_complex, "interactions/complex_portal/net_complex.tsv")


# STRING ------------------------------------------------------------------

string = tbl_df(read.table("interactions/string/9606.protein.links.v11.0.txt", sep=" ", skip=1))
string$V1 =  str_replace_all(string$V1, "9606\\.", "")
string$V2 =  str_replace_all(string$V2, "9606\\.", "")

mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
mapping_1 = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","ensembl_peptide_id"), mart=mart)

string_entrez = string
string_entrez$V1 = mapping_1$entrezgene[match(string$V1, mapping_1$ensembl_peptide_id)]
string_entrez$V2 = mapping_1$entrezgene[match(string$V2, mapping_1$ensembl_peptide_id)]
string_entrez = string_entrez[!(is.na(string_entrez$V1) | is.na(string_entrez$V2)),]
names(string_entrez) = c("gene_1", "gene_2", "score")
string_entrez = string_entrez %>% distinct(gene_1, gene_2)
write_tsv(string_entrez, "interactions/string/string_entrez.txt")

string_nodes_entrez = unique(unlist(string_entrez[,1:2]))
string_nodes = data.frame(entrez=string_nodes_entrez) %>% arrange(entrez)
string_nodes = filter(string_nodes, !is.na(entrez))
string_nodes = string_nodes %>% mutate(index=1:length(entrez)) %>% dplyr::select(index, entrez)
write_tsv(string_nodes, "interactions/string/string_nodes.txt", col_names=FALSE)

# convert edge nodes to indices
string_edges = string_entrez[,1:2]
string_edges$gene_1 = string_nodes$index[match(string_entrez$gene_1, string_nodes$entrez)]
string_edges$gene_2 = string_nodes$index[match(string_entrez$gene_2, string_nodes$entrez)]
write_tsv(string_edges, "interactions/string/string_edges.txt", col_names=FALSE)


# OMNIPATH ----------------------------------------------------------------

omnipath = read_tsv("interactions/omnipath/omnipath.txt")
uniprot = Rkeys(org.Hs.egUNIPROT)
mapping = select(org.Hs.eg.db, uniprot, "ENTREZID", "UNIPROT")

omnipath_entrez = data.frame(source=NULL, target=NULL, is_directed=NULL, is_stimulation=NULL, is_inhibition=NULL, dip_url=NULL)

for(i in 1:dim(omnipath)[1]) {
  
  entrez_in = mapping$ENTREZID[which(mapping$UNIPROT==omnipath$source[i])]
  entrez_out = mapping$ENTREZID[which(mapping$UNIPROT==omnipath$target[i])]
  if(length(entrez_in)>1 | length(entrez_out)>1) {
    print(i)
  }
  if(!(is_empty(entrez_in) | is_empty(entrez_out))) {
    for(j in 1:length(entrez_in)) {
      for(k in 1:length(entrez_out)) {
        omnipath_entrez = rbind(omnipath_entrez, data.frame(source=entrez_in[j], target=entrez_out[k], omnipath[i,3:6]))
      }
    }
  }
}

write_tsv(omnipath_entrez, "interactions/omnipath/omnipath_entrez.txt")

omnipath_nodes_entrez = sort(as.numeric(as.character(unique(unlist(omnipath_entrez[,1:2])))))
omnipath_nodes = data.frame(entrez=omnipath_nodes_entrez) %>% arrange(entrez)
omnipath_nodes = filter(omnipath_nodes, !is.na(entrez))
omnipath_nodes = omnipath_nodes %>% mutate(index=1:length(entrez)) %>% dplyr::select(index, entrez)
write_tsv(omnipath_nodes, "interactions/omnipath/omnipath_nodes.txt", col_names=FALSE)

# convert edge nodes to indices
omnipath_edges = omnipath_entrez[,1:2]
omnipath_edges$source = omnipath_nodes$index[match(omnipath_entrez$source, omnipath_nodes$entrez)]
omnipath_edges$target = omnipath_nodes$index[match(omnipath_entrez$target, omnipath_nodes$entrez)]
write_tsv(omnipath_edges, "interactions/omnipath/omnipath_edges.txt", col_names=FALSE)


# HURI --------------------------------------------------------------------

mapping = read_excel("gwas/HumanGeneList_17Sep2018_workup_betterensembl_list.xlsx")

huri = read_tsv("interactions/huri/HuRI.tsv")
huri_nodes = tbl_df(data.frame(ensembl=unique(as.character(unlist(huri[,1:2])))))
huri_nodes$entrez = mapping$EntrezGeneID[match(huri_nodes$ensembl, mapping$ENSEMBL_ID)]
huri_nodes = huri_nodes[!is.na(huri_nodes$entrez),]
huri_nodes$index = 1:dim(huri_nodes)[1]
huri_nodes %>% dplyr::select(index, entrez) %>% write_tsv("interactions/huri/huri_nodes.txt", col_names=FALSE)

huri_edges = data.frame(
  a = mapping$EntrezGeneID[match(huri$Ensembl_gene_id_a, mapping$ENSEMBL_ID)],
  b = mapping$EntrezGeneID[match(huri$Ensembl_gene_id_b, mapping$ENSEMBL_ID)]
)
huri_edges = huri_edges[apply(huri_edges, 1, function(x) !any(is.na(x))),]

write_tsv(huri_edges, "interactions/huri/huri_entrez.txt")

# convert edge nodes to indices
huri_edges$a = huri_nodes$index[match(huri_edges$a, huri_nodes$entrez)]
huri_edges$b = huri_nodes$index[match(huri_edges$b, huri_nodes$entrez)]
write_tsv(huri_edges, "interactions/huri/huri_edges.txt", col_names=FALSE)


# METABASE ----------------------------------------------------------------

metabase = read_tsv("/Volumes/am673712/links/network_analysis/interactions/metabase_plus_ppi/live/metabase_plus_ppi.tsv", col_names=FALSE)

metabase_nodes = tbl_df(data.frame(hgnc=unique(as.character(unlist(metabase)))))
metabase_nodes$entrez = mapping$entrezgene[match(metabase_nodes$hgnc, mapping$hgnc_symbol)]
metabase_nodes = metabase_nodes[!is.na(metabase_nodes$entrez),]
metabase_nodes$index = 1:dim(metabase_nodes)[1]
metabase_nodes %>% dplyr::select(index, entrez) %>% write_tsv("interactions/metabase_plus_ppi/metabase_nodes.txt", col_names=FALSE)

metabase_edges = data.frame(
  a = mapping$EntrezGeneID[match(metabase$X1, mapping$Symbol)],
  b = mapping$EntrezGeneID[match(metabase$X2, mapping$Symbol)]
)
metabase_edges = metabase_edges[apply(metabase_edges, 1, function(x) !any(is.na(x))),]

write_tsv(metabase_edges, "interactions/metabase_plus_ppi/metabase_plus_ppi_entrez.txt")

# convert edge nodes to indices
metabase_edges$a = metabase_nodes$index[match(metabase_edges$a, metabase_nodes$entrez)]
metabase_edges$b = metabase_nodes$index[match(metabase_edges$b, metabase_nodes$entrez)]
write_tsv(metabase_edges, "interactions/metabase_plus_ppi/metabase_edges.txt", col_names=FALSE)


# INTOMICS ----------------------------------------------------------------

intomics = read_graph("interactions/intomics/2019_04_01.high_confidence.graphml", format="graphml")
intomics_nodes = tbl_df(as.data.frame(get.vertex.attribute(intomics)))
intomics_nodes$index = 1:dim(intomics_nodes)[1]
intomics_nodes = intomics_nodes %>% separate_rows(`Ensembl.gene`, sep=",") 

mapping_gene = read_excel("gwas/HumanGeneList_17Sep2018_workup_betterensembl_list.xlsx")

intomics_nodes$entrez = mapping_gene$EntrezGeneID[match(intomics_nodes$Ensembl.gene, mapping_gene$ENSEMBL_ID)]
table(!is.na(intomics_nodes$entrez))

intomics_edges = tbl_df(as_data_frame(intomics, what="edges"))
intomics_edges$entrez_in = intomics_nodes$entrez[match(intomics_edges$from, intomics_nodes$index)]
intomics_edges$entrez_out = intomics_nodes$entrez[match(intomics_edges$to, intomics_nodes$index)]
intomics_edges = intomics_edges[apply(intomics_edges, 1, function(x) !any(is.na(x[7:8]))),]

intomics_nodes_filt = tbl_df(data.frame(entrez=intomics_nodes$entrez) %>% arrange(entrez))
intomics_nodes_filt = filter(intomics_nodes_filt, !is.na(entrez))
intomics_nodes_filt = intomics_nodes_filt %>% mutate(index=1:length(entrez)) %>% dplyr::select(index, entrez)
write_tsv(intomics_nodes_filt, "interactions/intomics/intomics_nodes.txt", col_names=FALSE)

# convert edge nodes to indices
intomics_edges_filt = intomics_edges[,7:8]
write_tsv(intomics_edges, "interactions/intomics/intomics_entrez.txt")
intomics_edges_filt$entrez_in = intomics_nodes_filt$index[match(intomics_edges_filt$entrez_in, intomics_nodes_filt$entrez)]
intomics_edges_filt$entrez_out = intomics_nodes_filt$index[match(intomics_edges_filt$entrez_out, intomics_nodes_filt$entrez)]
write_tsv(intomics_edges, "interactions/intomics/intomics_edges.txt", col_names=FALSE)


# LIGAND RECEPTOR ---------------------------------------------------------

# combine mark's metabase filter with existing data
mark_proxy = read_excel("interactions/mark_proxies/GeneProxies_v3_list.xlsx")
mark_proxy %>% filter(grepl("receptor",`Type of Proxy`,ignore.case=TRUE)) %>% group_by(Database) %>% summarise(n=n())

# filter
mark_proxy_filt = mark_proxy %>% filter(grepl("receptor",`Type of Proxy`,ignore.case=TRUE)) %>% dplyr::select(`GeneID of Interest`, `GeneID of Proxy`)
names(mark_proxy_filt) = c("ligand","receptor")
mark_proxy_filt = distinct(mark_proxy_filt)
write_tsv(mark_proxy_filt, "interactions/ligand_receptor_db/mark_proxy_filt.txt")


# MERGED NETWORK SOURCE ---------------------------------------------------

# merge non-metabase sources
# huri, string, omnipath

networks = as.list(c("huri","string","omnipath"))
net_dat = lapply(networks, function(x) read_tsv(paste0("interactions/",x,"/",x,"_entrez.txt")))
for(i in 1:length(net_dat)) names(net_dat[[i]])[1:2] = c("source","target")
names(net_dat) = networks
net_dat = bind_rows(net_dat, .id="data_source")
write_tsv(net_dat, "interactions/merged_network/net_dat.txt")


