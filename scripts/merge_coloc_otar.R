
# load libraries
library(gtx)
library(tidyverse)
library(implyr)
library(glue)
library(pryr)


# COLOCALISATION ----------------------------------------------------------

# set global default params
db = "gene_gwas_use"

# Establish DB connection
dbc <- validate_impala()

# reference the necessary tables
analyses_tbl =  tbl(dbc, glue('{db}.analyses'))
sites_tbl =  tbl(dbc, glue("{db}.sites_ukb_500kv3"))

genes_tbl = tbl(dbc, sql(glue('SELECT * FROM {db}.genes WHERE genetype = "protein_coding"')))

gwas_th_tbl = tbl(dbc, sql(glue('SELECT * FROM {db}.gwas_results_top_hits WHERE signal IS NOT NULL AND chrom IS NOT NULL AND pval_index IS NOT NULL AND analysis ilike "%500kv3%"')))

colocs_raw_tbl <- tbl(dbc, sql(glue('SELECT * FROM {db}.coloc_results WHERE analysis2 ilike "%500kv3%" AND p12 IS NOT NULL'))) 

# filter for reasonable colocs
# reasonable = gwas in cis-window, gwas pval <= 5e-8
colocs_tbl <-   
  colocs_raw_tbl %>% 
  
  # append gene info (hgncid, tss start, genetype)
  # use inner join b/c we only want protein coding colocs
  # tss used to ensure gwas signal is in the cis-window (tss +/- 1mb)
  inner_join(genes_tbl %>%
               select(ensemblid,
                      gene_start=pos_start,
                      genetype,
                      hgncid,
                      chrom),
             by = c("entity"="ensemblid")) %>%
  
  # append gwas top hits info as th_*
  inner_join(gwas_th_tbl %>%
               select(chrom, pos_start, pos_end, pos_index,
                      ref_index, alt_index, pval_index, analysis),
             by = c("analysis2"="analysis", "chrom")) %>%
  
  # reorder and rename cols
  select(everything(),
         th_pos = pos_index, th_pval = pval_index,
         th_start = pos_start, th_end = pos_end,
         th_ref = ref_index, th_alt = alt_index) %>%
  
  # append rsid to each gwas th
  # left join b/c not all chrom_pos_ref_alt have rsid
  left_join(sites_tbl %>%
              select(rs_chrom = "chrom", rs_pos = "pos",
                     rs_ref = "ref", rs_alt = "alt", rs),
            by = c("chrom"  = "rs_chrom", 
                   "th_pos" = "rs_pos",
                   "th_ref" = "rs_ref",
                   "th_alt" = "rs_alt")) %>%
  
  # join gwas trait info - e.g. ncase
  # inner join is okay b/c all analysis ids should be present in both
  inner_join(analyses_tbl %>% select(analysis, ncase, ncohort),
             by = c("analysis2" = "analysis")) %>%
  
  # filter out gwas with low numbers
  filter(ncase >= 200 | ncohort >= 200) %>%
  
  # filter for colocs in cis-window
  filter((gene_start - 1e6 < th_start) &
           (gene_start + 1e6 > th_end) & 
           (th_pval < 5e-8)) %>%
  
  # calc distance b/w each th & gene tss
  # this will be used to help resolve "ties" later
  mutate(th2tss_dist = abs(th_pos - gene_start)) %>%
  
  # calc max_p12 for each gene-gwas & gene-gwas-th to resolve ties later
  group_by(entity, analysis2) %>%
  mutate(gene_gwas_max_p12 = max(p12, na.rm = TRUE)) %>%
  mutate(gene_gwas_min_th2tss_dist = min(th2tss_dist, na.rm = TRUE)) %>%
  
  # for each target indication
  group_by(entity, analysis2, chrom, th_pos, th_ref, th_alt) %>%
  
  # calculate the max h4 (p12) across all tissues
  mutate(gene_gwas_th_max_p12 = max(p12, na.rm = TRUE)) %>%
  ungroup()

colocs <- colocs_tbl %>% collect()

write_tsv(colocs, "v2g/colocs.txt")


# OTAR V2G ----------------------------------------------------------------

# set global default params
db = "genetics_portal"

# Establish DB connection
dbc <- validate_impala()

otar_tbl = tbl(dbc, sql(glue('SELECT * FROM {db}.variant_to_gene_20181128v2')))
otar_tbl %>% head()

