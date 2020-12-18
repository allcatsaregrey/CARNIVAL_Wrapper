# Make sure work directory is appropriately configured 
source("pkn_generator.R")

runPipeline <- function(net_sources, sources, targets) {
  
  # Input the data and generate the necessary PKN
  user_data <- proc_data()
  user_net <- init_net(net_sources) 
  user_pkn <- gen_pkn(user_net, sources, targets) %>% process_net()
  
  # Run CARNIVAL and display the resultant network
  result <- carn_pkn(user_net, user_data[[1]], user_data[[2]])
  vis_carn(result) 
  
}