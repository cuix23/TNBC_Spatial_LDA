#' Model assessment for LDA based on the mean abundance for each taxon.
#'
#' @param stan.fit An instance of stanfit.
#' @param ASVsIndexToPlot An integer vector. Use to select ASVs to show the goodness of fit in histograms.
#' @inheritParams alignmentMatrix
#' @return A ggplot2 object. Histogram of data generated from the posterior estimates and the observed data.
#' @importFrom rstan stan
#' @importFrom ggplot2 facet_wrap geom_vline theme_update aes
#' @import phyloseq
#' @importFrom reshape2 melt
#' @export
#'
modelAssessment_mean <- function(
    #spe,
  dtm,
  #stan.fit = NULL,
 res,
  warm_up_iter = NULL,
  iter = 2000,
  cellTypeIndexToPlot = c(1:16)
){
  
  value <- Var2 <- NULL
  
  # determine the iteration used in posterior sampling (subtract warm up iterations)
  if (is.null(warm_up_iter)) {
    iterUse = iter / 2
  } else {
    iterUse = iter - warm_up_iter
  }
  
  x <- dtm
  
  # draws from posterior predictive distribution
  x_sim <- res$x_sim[1:iterUse, , ] # iteration * tissues * topics
  
  # choose only the first chain
  mean_all <- apply(x_sim[1, , ], 2, mean)
  
  # find the mean of each cell types in each iteration acorss sample 
  for (i in 2:iterUse){
    mean_x_sim_i <- apply(
      x_sim[i, ,], 
      2,
      mean)
    mean_all <- rbind(mean_all, mean_x_sim_i)
  }
  
  rownames(mean_all) <- c(
    paste0(
      "x_mean_rep",
      seq(1, iterUse)
    )
  )
  
  colnames(mean_all) <- colnames(x)
  #return(mean_all)
  
  # subset the interested phenotype
  mean_all <- mean_all[, cellTypeIndexToPlot]
  mean_all_long <- reshape2::melt(mean_all)
  
  
  # finding the observed mean value 
  x_mean_cellCount <- data.frame(
    Var1 = rep(
      "x_mean_obs",
      dim(x)[2]
    ),
    Var2 = colnames(x),
    count = apply(x, 2, mean)
  )
  
  # subset the interested phenotype
  x_mean_cellCount <- x_mean_cellCount[cellTypeIndexToPlot, ]
  
  
  # plotting
  p_hist <- ggplot2::ggplot(
    data = mean_all_long
  ) +
    ggplot2::geom_histogram(
      aes(
        x = value,
        group = Var2
      ),
      color = "#0072B2",
      fill = "#0072B2",
      bins = 50) +
    ggplot2::xlab("mean")+
    
    ggplot2::facet_wrap(~Var2, nrow = 4, scales = "free_x") +
    
    ggplot2::geom_vline(
      data = x_mean_cellCount,
      aes(xintercept = count),
      color = "#CC79A7"
    ) +
    ggplot2::theme_update(
      text = element_text(size = 8)
    )
  return(p_hist)
  
}
