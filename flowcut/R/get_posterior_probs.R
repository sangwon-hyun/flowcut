# Generated from _main.Rmd: do not edit by hand

#' Get the membership probabilities of a new dataset |ynew| given GMM model |gmm_model|.
#'
#' @param gmm_model from mclust or me.weighted
#' @param ynew new data to get predicted posterior membership for.
#'
#' @return estimated posterior probabilities from an estimated GMM model
get_posterior_probs <- function(gmm_model, ynew){
  stopifnot(ncol(ynew)==1)
  numclust = length(gmm_model$parameters$mean)
  mn = gmm_model$parameters$mean
  sd = gmm_model$parameters$variance$sigmasq %>% sqrt()
  pro = gmm_model$parameters$pro
  post_prob = lapply(1:numclust, function(iclust){
    weighted_d1 = dnorm(ynew, mean = mn[iclust], sd = sd[iclust]) * pro[iclust]
  })  %>% do.call(cbind, .)
  post_prob = post_prob/rowSums(post_prob)
  return(post_prob)
}
