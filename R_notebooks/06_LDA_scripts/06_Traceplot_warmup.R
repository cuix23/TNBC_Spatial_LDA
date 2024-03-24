library(ggplot2)
library(dplyr)

traceplot_warmup <- function(stan.fit,
                              spe,
                              K = 5,
                              warm_up_iter = NULL,
                              iterations = 2000,
                              chain = 4,
                              # cor_method = c("cor", "cosine"),
                              TissueID_name = "tissue_id") {
  # theta_aligned_with_warmup
  ################### Theta with warmup ###################
  res <- rstan::extract(
    stan.fit,
    permuted = FALSE,
    inc_warmup =TRUE,
    include = TRUE
  )
  # Let's assume 'theta_index' is the position of 'theta' in your Stan model's parameters
  theta_index <- which(colnames(res[, 1, ]) == "theta")
  # Assuming 'res' is a three-dimensional array from rstan::extract() with permuted = FALSE
  # Find indices of all columns that include 'theta' in their names
  theta_indices <- grep("theta", colnames(res[, 1, ]))
  
  theta <- res[, , theta_indices]
  
  
  ################### Alignment matrix ###################
  Chain <- Topic <- topic.dis <- NULL
  # theta is a 3 dimensional array (iterations * Tissues * topic)
  
  # determine the iteration used in posterior sampling (subtract warm up iterations)
  if (is.null(warm_up_iter)) {
    iterUse = iter / 2
  } else {
    iterUse = iter - warm_up_iter
  }
  
  
  
  ## acquire Tissue identity
  dimnames(theta)[[2]] <- unique(spe$tissue_id)
  ## acquire topic number
  dimnames(theta)[[3]] <- c(paste0("Topic_", seq(1,K)))
  
  # array to a dataframe
  ## melt() takes wide-format data and melts into long-format data
  theta_all = reshape2::melt(theta)
  colnames(theta_all) = c("iteration", "Tissue", "Topic", "topic.dis")
  theta_all$Chain = paste0("Chain ", rep(rep(seq(1, chain), each = iterUse), times = K*length(unique(spe$tissue_id))))
  
  theta_all$Topic = factor(theta_all$Topic)
  theta_all$Chain = factor(theta_all$Chain)
  theta_all$Tissue = as.character(theta_all$Tissue)
  #theta_all$Tissue = factor(theta_all$Tissue)
  
  # join the SpatialExperimet object column data to theta_all
  ## colData gets the metadata (clinical data)
  sam = (colData(spe)[, c(1:9, 54)]
         |> unique()
         |> data.frame()
  )
  
  theta_all = dplyr::left_join(
    theta_all,
    sam,
    by = c("Tissue"= TissueID_name)
  )
  
  
  # align matrix
  aligned <- matrix(nrow = K, ncol = chain)
  aligned[, 1] <- seq(1, K)
  corrTop <- numeric()
  for(j in 1:K) { #Topic of first chain
    chains <- lapply(as.list(1:chain), function(x){
      for(top in 1:K){
        corrTop[top] <- cor(
          theta_all |> dplyr::filter(Chain == "Chain 1") |> dplyr::filter(Topic == paste0("Topic_", j)) |> dplyr::select(topic.dis),
          theta_all |> dplyr::filter(Chain == paste0("Chain ",x)) |> dplyr::filter(Topic == paste0("Topic_", top)) |> dplyr::select(topic.dis))
      }
      return(which(corrTop == max(corrTop)))
    })
    
    aligned[j, 2:chain] <- unlist(chains)[2:chain]
  }
  
 
  
  ################### theta_aligned_with_warmup ###################

    # align topics between chains
    # switch samples and Topic dimension in array
    theta <- aperm(theta, c(1,3,2))
    
    theta_chain <- list()
    theta_chain[[1]] <- theta[1:(iter), ,]
    theta_chain[[1]] <- theta_chain[[1]][, aligned[,1],]
    
    theta_aligned <- theta_chain[[1]]
    
    for(ch in 2:chain){
      theta_chain[[ch]] <- theta[((ch-1)*(iter)+1):(ch*(iter)),,]
      theta_chain[[ch]] <- theta_chain[[ch]][,aligned[,ch],]
      
      theta_aligned <- abind::abind(theta_aligned, theta_chain[[ch]], along = 1)
    }

    # theta_aligned <- abind(theta_chain[[1]], theta_chain[[2]], along = 1)
    # for(ch in 3:chain){
    #   theta_aligned <- abind(theta_aligned, theta_chain[[ch]], along = 1)
    # }
    # switch back samples and Topic dimension in array
    theta_aligned <- aperm(theta_aligned, c(1,3,2))
    
    

  
  ################### plot ###################
  
  # Assuming theta_aligned is a 3D array: iterations x tissues x topics
  # Split theta_aligned by chain
  chains <- lapply(1:chain, function(chain) {
    start <- (chain - 1) * iterations + 1
    end <- chain * iterations
    theta_aligned[start:end, , ]
  })
  
  # Prepare data for plotting
  plot_data <- list()
  
  for (chain_index in 1:chain) {
    for (tissue in 1:dim(theta_aligned)[2]) {
      for (topic in 1:dim(theta_aligned)[3]) {
        # Extract data for this tissue-topic combination from the current chain
        temp_data <- data.frame(
          Iteration = 1:iterations,
          Theta = chains[[chain_index]][, tissue, topic],
          Chain = paste("Chain", chain_index),
          TissueTopic = paste("Tissue", tissue, "Topic", topic)  # Combining tissue and topic for unique facetting
        )
        plot_data[[length(plot_data) + 1]] <- temp_data
      }
    }
  }
  
  # Combine all data into a single data frame
  plot_data <- do.call(rbind, plot_data)
  
  # Generate the plot
  p <- ggplot(plot_data, aes(x = Iteration, y = Theta, color = Chain)) +
    geom_line() +
    facet_wrap(~TissueTopic, scales = "free_y",ncol = 5) +  # Using facet_wrap for individual scaling
    labs(x = "Iteration", y = "Theta Value") +
    theme_minimal() +
    # scale_color_manual(values = c("Chain 1" = "blue", "Chain 2" = "red")) +
    theme(strip.text = element_text(size = 8))
  
  return(p)
}
