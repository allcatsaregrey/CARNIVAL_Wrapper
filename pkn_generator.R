# Authors: TJ McColl, Matthew McFee
# A simple set of functions to generate prior knowledge networks using the
# OmnipathDB web database interface and then perform CARNIVAL analysis with said 
# PKN. This is mostly adapted form the examples provided by the Saez lab group.

# TODO: Streamline adding network information for nodes that do not exist in
# the OmnipathDB network. Add a function to visualize graph networks in 3D and 
# function to look at subgraphs of generated PKNs in 2D. 

# Library imports
req_pac <- c("dplyr", "ggplot2", "OmnipathR", "igraph", "ggraph", "GGally", 
             "tidyverse", "BiocManager")

bioc_pac <- c("Rgraphviz", "dorothea", "progeny", "viper", "Seurat", "biomaRt",
              "GEOquery", "CARNIVAL")

# Force package install
for (pac in req_pac) {
  if(! pac %in% installed.packages()) {
    install.packages(pac, dependencies = TRUE)
  }
}

# Force package install from Bioconductor
for (pac in bioc_pac) {
  if(! pac %in% installed.packages()) {
    BiocManager::install(pac)
  }
}

# Load necessary libraries
req_pac <- c(req_pac, bioc_pac)
lapply(req_pac, require, character.only = TRUE)

# click on the 'interactions' variable to view the proteins included

init_net <- function(sources){
  
interactions <- import_omnipath_interactions(resources = sources)
user_net <-  interaction_graph(interactions = interactions)

}

net_head <- function(interactions, n) {
  
print_interactions(head(interactions, n))
  
}

proc_data <- function(geo_id, find_weights = FALSE) {
  
  # See dbellViper via data(bcellViper, package = "bcellViper") for sample 
  user_dat <- getGEO(geo_id, GSEMatrix = TRUE)
  
  # Select regulons with high confidence values of "A" and "B"
  regulons <- data(dorothea_hs) %>% filter(confidence %in% c("A", "B"))
  
  # Bulk RNAseq data can be processed this way
  tf_activities <- run_viper(user_dat, regulons, 
                             options =  list(method = "scale", minsize = 4, 
                                             eset.filter = FALSE, cores = 1, 
                                             verbose = FALSE)) %>% 
    as(tf_activities, "data.frame")
  
  # Bulk RNAseq data processing for PROGENy
  # Import data to DESeq2 and variance stabilize
  
  prog_pathways <- list()
  
  if(find_weights == TRUE) {
    dset <- DESeqDataSetFromMatrix(user_dat,
                                colData = as.data.frame(colData(prog_dat)), 
                                design =~ dex)
    dset <- estimateSizeFactors(dset)
    dset <- estimateDispersions(dset)
    gene_expr <- getVarianceStabilizedData(dset)
    
    # Annotate matrix with HGNC symbols
    mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes = getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                values=rownames(gene_expr), mart = mart)
    matched = match(rownames(gene_expr), genes$ensembl_gene_id)
    rownames(gene_expr) = genes$hgnc_symbol[matched]
    
    # Obtain pathway scores
    pathways <- progeny(gene_expr, scale=FALSE)
    # Need to adjust name of control column as necessary
    controls <- prog_dat$dex == "untrt"
    ctl_mean <- apply(pathways[controls,], 2, mean)
    ctl_sd <- apply(pathways[controls,], 2, sd)
    pathways <- t(apply(pathways, 1, function(x) x - ctl_mean))
    prog_pathways <- apply(pathways, 1, function(x) x / ctl_sd)
    
    return(tf_activities, prog_pathways)
    
  }
  
  # Return the data for use in CARNIVAL
  list(tf_activities, prog_pathways)
  
}

remove_nonexist <- function(user_net, sources, targets) {
  
  # L[!(L %in% L1)]
  sources <- sources[sources %in% as_ids(V(user_net))]
  targets <- targets[targets %in% as_ids(V(user_net))]
  
  list(sources, targets)
  
}

gen_pkn <- function(user_net, sources, targets) {
  
  # Remove non-existent nodes before generating the paths
  pruned_src_trg <- remove_nonexist(user_net, sources, targets)
  sources <- pruned_src_trg[[1]]
  targets <- pruned_src_trg[[2]]
  
  collected_path_nodes = list()
  
  for(source in 1:length(sources)){
    
    paths <- shortest_paths(user_net, from = sources[[source]],
                            to = targets,
                            output = 'vpath')
    path_nodes <- lapply(paths$vpath,names) %>% unlist() %>% unique()
    collected_path_nodes[[source]] <- path_nodes
  }
  
  collected_path_nodes <- unlist(collected_path_nodes) %>% unique()
  
  # Generate the prior knowledge network from the paths we generated
  net_nodes <- c(sources, targets, collected_path_nodes) %>% unique()
  generated_pkn <- induced_subgraph(graph = user_net, vids = net_nodes)
  
}

gen_undirected <- function(user_graph, target) {
  
user_graph_undirected <- as.undirected(user_graph, mode = c("mutual"))
user_graph_undirected <- simplify(user_graph_undirected)
cl_results <- cluster_fast_greedy(user_graph_undirected)

# We extract the cluster where a protein of interest is contained
cluster_id <- cl_results$membership[which(cl_results$names == target)]
sub_graph <- induced_subgraph(user_graph_undirected,
                                  V(user_graph_undirected
                                    )$name[which(cl_results$membership == cluster_id)])

}

display_graph <- function(user_graph) {

# Add user customization options for improved visualization  
ggnet2(user_graph, label=TRUE)

}

# Process the network for input into CARNIVAL (ie. convert an igraph object
# to an appropriate data frame)
process_net <- function(user_graph) {
  
  # Convert the graph into a data frame for manipulation
  conv_graph <- as_data_frame(user_graph) %>% 
    select(c("from", "to", "consensus_direction")) %>%
    rename(from = Source, consensus_direction = Effect, to = Target)
}


# Process the network, and user measurements using carnival
carn_pkn <- function(user_net, user_measure, user_weight = NULL) {
  
  result = runCARNIVAL(measObj = user_measure, netObj = user_net,
                       weightObj = user_weight, dir_name = getwd())
  
}

# Custom visualization for CARNIVAL results
alt_vis_carn <- function(result) {
  
  # Prepare the result data for plotting
  user_res <- result$weightedSIF
  attribs <- as.data.frame(result$nodesAttributes)
  user_res <- user_res[ , c(1, 3, 2, 4)] %>% as.data.frame()
  user_res$Effect <- result$nodesAttributes
  user_graph <- graph.data.frame(user_res, directed = TRUE, vertices = NULL)
  
  # Visualize the result
  plot(user_graph, vertex.color = attribs$AvgAct, vertex.frame.color = "Black",
       edge.color = "Black", vertex.size = 50, vertex.label.color = "Black")
  
}


# Move the wrapper to this file since writing an R package would be overkill
# in this case 
runPipeline <- function(net_sources, sources, targets, geo_id) {
  
  # Input the data and generate the necessary PKN
  user_data <- proc_data(geo_id)
  user_net <- init_net(net_sources) 
  user_pkn <- gen_pkn(user_net, sources, targets) %>% process_net()
  
  # Run CARNIVAL and display the resultant network
  result <- carn_pkn(user_net, user_data[[1]], user_data[[2]])
  vis_carn(result) 
  
}