require(dplyr)
require(readr)
require(stringr)
require(sparklyr)

sc <- spark_connect(
  master = "yarn-client",
  spark_home = "/opt/cloudera/parcels/SPARK2/lib/spark2",
  version = "2.3.0"
) 

tbl_change_db(sc, "genetics_portal")
src_tbls(sc)

dat = sc %>% tbl("variant_to_gene_20181128_scored_overall")

gene_ids = dat %>% select(gene_id) %>% distinct %>% collect()
gene_ids = sort(as.character(unlist(gene_ids)))
gene_ids_list = split(gene_ids, ceiling(seq_along(gene_ids)/30))

for(i in 80:length(gene_ids_list)) {
  
  print(i)
  dat_filt = dat %>% filter(gene_id %in% gene_ids_list[[i]])
  dat_filt_collect = dat_filt %>% collect()
  dat_filt_collect = dat_filt_collect %>% arrange(gene_id, variant_id)
  write_tsv(dat_filt_collect, paste0("dat_", str_pad(i,4,pad="0"), ".txt"))
  
}

