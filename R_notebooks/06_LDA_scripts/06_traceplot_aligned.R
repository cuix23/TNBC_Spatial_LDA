library(ggplot2)
library(dplyr)

traceplot_aligned <- function(theta_aligned, num_chains, num_iterations_per_chain) {
  # Assuming theta_aligned is a 3D array: iterations x samples x topics
  # Split theta_aligned by chain
  chains <- lapply(1:num_chains, function(chain) {
    start <- (chain - 1) * num_iterations_per_chain + 1
    end <- chain * num_iterations_per_chain
    theta_aligned[start:end, , ]
  })
  
  # Prepare data for plotting
  plot_data <- list()
  
  for (chain_index in 1:num_chains) {
    for (sample in 1:dim(theta_aligned)[2]) {
      for (topic in 1:dim(theta_aligned)[3]) {
        # Extract data for this sample-topic combination from the current chain
        temp_data <- data.frame(
          Iteration = 1:num_iterations_per_chain,
          Theta = chains[[chain_index]][, sample, topic],
          Chain = paste("Chain", chain_index),
          sampleTopic = paste("Sample", sample, "Topic", topic)  # Combining sample and topic for unique facetting
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
    facet_wrap(~SampleTopic, scales = "free_y",ncol = 5) +  # Using facet_wrap for individual scaling
    labs(x = "Iteration", y = "Theta Value") +
    theme_minimal() +
    # scale_color_manual(values = c("Chain 1" = "blue", "Chain 2" = "red")) +
    theme(strip.text = element_text(size = 8))
  
  return(p)
}
