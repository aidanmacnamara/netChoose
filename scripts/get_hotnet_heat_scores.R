require(tidyverse)
require(rhdf5)
require(parallel)

dir_hotnet2 = "/home/am673712/hotnet2"
input_dir = '/home/am673712/network_analysis/hotnet/hotnet/test_all_gwas/input/omnipath'
result_dir = '/home/am673712/network_analysis/hotnet/hotnet/test_all_gwas/results/omnipath'
node_file = "interactions/omnipath/omnipath_nodes.txt"
edge_file = "interactions/omnipath/omnipath_edges.txt"

num_cores = 10
mclapply(as.list(1:712), get_heat, mc.cores=num_cores) 

get_heat <- function(file_ix) {
  
  input_file = paste0(input_dir, "/input_", file_ix, ".tsv")
  output_dir = paste0(result_dir, "/output_", file_ix)
  
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  # get input score and node indices
  
  input_score = read_tsv(input_file, col_names=FALSE)
  nodes = read_tsv(node_file, col_names=FALSE)
  
  # create influence matrix
  
  cmd_1 = paste0("/bin/python2.7.11 ", dir_hotnet2, "/scripts/createInfluenceMatrix.py -e ", edge_file, " -i ", node_file, " -n omnipath -o ", output_dir, "/ppr.h5 -v 1 hotnet2 -b 0.5")
  system(cmd_1)
  
  h5ls(paste0(output_dir, "/ppr.h5"))
  ppr = h5read(paste0(output_dir, "/ppr.h5"),"PPR")
  ppr_nodes = h5read(paste0(output_dir, "/ppr.h5"),"nodes")
  
  # filter and write
  nodes_filt = nodes[nodes$X2 %in% ppr_nodes,]
  nodes_filt$X1 = 1:length(nodes_filt$X1)
  input_score_filt = input_score[input_score$X1 %in% ppr_nodes,]
  
  write_tsv(input_score_filt, paste0(output_dir, "/input_filt.tsv"), col_names=FALSE)
  write_tsv(nodes_filt, paste0(output_dir, "/nodes_filt.tsv"), col_names=FALSE)
  
  cmd_2 = paste0("/bin/python2.7.11 ", dir_hotnet2, "/compute_smoothed_scores.py -mf ", output_dir, "/ppr.h5 -igf ", output_dir, "/nodes_filt.tsv -gsf ", output_dir, "/input_filt.tsv -o ", output_dir, "/heat.txt")
  
  system(cmd_2)
  
}

# check result

heat_score = read_tsv(paste0(output_dir, "/heat.txt"), col_names=FALSE)

comp = data.frame(
  input = input_score$X2,
  heat = heat_score$X2[match(input_score$X1, heat_score$X1)]
) %>% tbl_df()

comp %>% ggplot(aes(input, heat)) + geom_point() + theme_thesis() # + coord_cartesian(xlim=c(9.75,10), ylim=c(4.95,5.1))

