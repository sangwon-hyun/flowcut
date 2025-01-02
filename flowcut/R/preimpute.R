# Generated from _main.Rmd: do not edit by hand

#' At each time point, we will (1) set aside the censored points, and (2) estimate
#' a 2-component GMM (2) take the censored points, and "impute" them according to those
#' mixtures. (After that, we'll (3) fit flowmix on the resulting imputed data.)
#'
#' This is only designed for 1-dimensional data with two clusters.
#'
#' @param datobj List with ylist (particles), countslist, X, and Cbox.
#'
#' @return datobj with imputed ylist.
#'
#' @export
preimpute <- function(datobj){

  ## Setup
  stopifnot(ncol(datobj$ylist[[1]]) == 1)
  stopifnot("Cbox" %in% names(datobj))

  ## Split the data
  splitres = split_data(datobj)
  inner_ylist = splitres$inner_ylist
  inner_countslist = splitres$inner_countslist
  censored_ylist = splitres$censored_ylist
  censored_direction = splitres$censored_direction
  censor.which.list = splitres$censor.which.list

  ## What happens when the "inner" data is empty-ish at some times? We sample
  ## some points from the rest of the data, and fill the empty-ish cytograms in
  ## so that there are 5 particles.
  inner_ntlist = sapply(inner_ylist, nrow)
  min_num_particles = 5
  times_empty = which(inner_ntlist < min_num_particles)
  collapsed_ylist =
    inner_ylist[-times_empty] %>%
    lapply(function(y){sample(y, size=ceiling(nrow(y)/3), replace=FALSE)}) %>% 
    Reduce(c, .)

  if(length(times_empty) > 0){
    for(tt_empty in times_empty){
      nt = length(inner_ylist[[tt_empty]])
      extra_points = sample(collapsed_ylist, min_num_particles - nt) %>% cbind()
      inner_ylist[[tt_empty]] <- rbind(inner_ylist[[tt_empty]], extra_points)
      inner_countslist[[tt_empty]] <- c(inner_countslist[[tt_empty]], rep(1, length(extra_points)))
    }
  }
  stopifnot(all(sapply(inner_ylist, nrow)>=min_num_particles))

  ## Fit 2-mixture Gaussian at each time point.
  gmm_model_list <- mapply(function(y, counts){

    ## z is just a random initialization of the memberships of the points.
    z = mclust::unmap(sample(1:2, size = length(y), replace=TRUE))

    ## Fit a Gaussian mixture model for weighted data.
    ##res = mclust::me.weighted(as.numeric(y), modelName = "V", weights = counts, z = z)
    res = mclust::Mclust(as.numeric(y), G=2, modelName = "V", verbose = FALSE, warn = FALSE)##, weights = counts, z = z)
    if(is.null(res)){
      res = mclust::Mclust(as.numeric(y), G=2, modelName = "E", verbose = FALSE, warn = FALSE)##, weights = counts, z = z)
      res$parameters$variance$sigmasq = rep(res$parameters$variance$sigmasq,2)
    }
    assertthat::assert_that(!is.null(res))
    return(res)
  }, inner_ylist, inner_countslist, SIMPLIFY = FALSE)

  ## Impute the points
  imputed_censored_ylist = mapply(function(gmm_model, ynew, dirs){

    ## Impute based on 2-cluster model
    if(nrow(ynew) == 0) return(ynew)
    ynew_imputed = ynew
    post_probs = get_posterior_probs(gmm_model, ynew)
    new_mems = post_probs %>% apply(1, function(ps){ sample(1:2, 1, prob=ps) })

    ## Get the gaussian means and SDs
    means = gmm_model$parameters$mean
    sds = gmm_model$parameters$variance$sigmasq %>% sqrt()

    ## Information for imputing each points
    impute_mat = tibble(new_mems = new_mems, dir = dirs)

    ## Impute the points
    for(irow in 1:nrow(impute_mat)){
      iclust = impute_mat$new_mems[irow]
      if(impute_mat$dir[irow] == 1){
        imputed_point = truncnorm::rtruncnorm(1, a = datobj$Cbox[1,2], mean = means[iclust], sd = sds[iclust])
      }
      if(impute_mat$dir[irow] == -1){
        imputed_point = truncnorm::rtruncnorm(1, b = datobj$Cbox[1,1], mean = means[iclust], sd = sds[iclust])
      }
      ynew_imputed[irow] = imputed_point
    }
    return(ynew_imputed)
  }, gmm_model_list, censored_ylist, censored_direction, SIMPLIFY = FALSE)
  stopifnot(all(sapply(imputed_censored_ylist, nrow) == sapply(censored_ylist, nrow)))

  ##  Return the datobj (ylist, countslist, and flowmix).
  TT = length(datobj$ylist)
  new_ylist = datobj$ylist
  for(tt in 1:TT){
    inds = censor.which.list[[tt]]
    if(length(inds) == 0) next
    ynew = datobj$ylist[[tt]]
    ynew[inds] = imputed_censored_ylist[[tt]]
    new_ylist[[tt]] = ynew
  }
  stopifnot(all(sapply(new_ylist, nrow) == sapply(datobj$ylist, nrow)))
  new_datobj = datobj
  new_datobj$ylist = new_ylist
  return(new_datobj)
}
