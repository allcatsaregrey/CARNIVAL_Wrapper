# Authors: TJ McColl, Matthew McFee
# A simple set of functions to generate prior knowledge networks using the
# OmnipathDB web database interface and then perform CARNIVAL analysis with said 
# PKN. This is mostly adapted form the examples provided by the Saez lab group.

# Library imports
req_pac <- c("dplyr", "ggplot2", "OmnipathR", "igraph", "ggraph", "GGally", 
             "CARNIVAL", "Rgraphviz")
lapply(req_pac, require, character.only = TRUE)

# click on the 'interactions' variable to view the proteins included

init_net <- function(sources){
  
interactions <- import_omnipath_interactions(filter_databases = sources)
user_net <-  interaction_graph(interactions = interactions)

}

net_head <- function(interactions, n) {
  
print_interactions(head(interactions, n))
  
}

gen_pkn <- function(user_net, sources, targets) {
  
  collected_path_nodes = list()
  
  for(source in 1:length(sources)){
    
    paths <- all_shortest_paths(user_net, from = sources[[source]],
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

# Read in appropriate data 
read_user_data <- function() {
  
  user_measure <- read_csv(file.choose())
  user_input <- read_csv(file.choose())
  
  return(c(user_measure, user_input))
  
}

# Process the network, and user measurements using carnival
carn_pkn <- function(user_net, user_measure, user_input) {
  
  result = runCARNIVAL(inputObj = user_input, measObj = user_measure,
                       netObj = user_net)
  
}

# Visualize the output of the CARNIVAL pipeline with RGraphviz
vis_carn <- function(result) {
  
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