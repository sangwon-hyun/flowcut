# Generated from _main.Rmd: do not edit by hand

#' Generate 1d data with 2 clusters from a list (|true_model|)
#' containing true model parameters.
#'
#' @param true_model List containing beta, alpha, mn, prob, numclust.
#' @param nt Particles per time point.
#'
#' @return Cytograms (a |ylist| object)
#' @export
gen_1d <- function(true_model, nt = 1000){

  ## Setup
  stopifnot(true_model$numclust == 2)
  TT = dim(true_model$mn)[1]

  ## Generate cytograms
  ylist = list()
  for(tt in 1:TT){
  
    ## Generate memberships Samples |nt| memberships out of (1:numclust)
    ## according to the cluster probabilities in |prob|.
    nt_by_clust = stats::rmultinom(1, size = nt, true_model$prob[tt,])
    ## draws = sample(1:numclust, size = nt, replace = TRUE, prob = true_model$prob[tt,])
    draws = c(rep(1, nt_by_clust[1]), rep(2, nt_by_clust[2]))
  
    y_onetime = list()
    for(iclust in 1:true_model$numclust){
      ntk = nt_by_clust[iclust]
      membership = rep(iclust, ntk)
      y_onetime[[iclust]] = cbind(MASS::mvrnorm(n = ntk,
                                                mu = true_model$mn[tt,,iclust],
                                                Sigma = true_model$sigma[iclust,,]))
    }
    y = do.call(rbind, y_onetime)
  
    ## Data
    ylist[[tt]] = y
  }
  return(ylist)
}
