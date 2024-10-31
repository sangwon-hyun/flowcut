# Generated from _main.Rmd: do not edit by hand

#' Reformatting the results from MCMC to create a flowmix-like object that
#' contains of the posterior means of the model parametr. This allows you to use
## plotting code from the flowmix and flowtrend package.
#' 
#' @param res The result of doing MCMC.
#'
#' @return posterior mean estimates of the cluster coefficients, mean, variance
#'   and probability.
#' 
mcmc_res_to_flowmix <- function(res, last_draws_inds=NULL){

  ## Setup
  X = t(res$dat.info$X)
  TT = nrow(X)
  p = ncol(X)
  dimdat = res$dat.info$ylist[[1]] %>% ncol()
  obj = list()
  Nmc = res$pos.beta %>% dim() %>% tail(1)
  if(is.null(last_draws_inds)) last_draws_inds = 1:Nmc
  if(!is.null(last_draws_inds))stopifnot(all(last_draws_inds<=Nmc))
  ## last_draw_inds = 1:dim(res$pos.beta)[4]

  ## Means
  pos.mn <- list()
  for(kk in 1:numclust){
    last_draws_inds
    pos.mn[[kk]] <- mclapply(last_draws_inds, function(mm){
      res$pos.beta[,,kk,mm,drop=TRUE]%*% rbind(1,datobj$X)
    },
      mc.cores = detectCores())%>%
      abind::abind(.,along=3)
  }
  post.mn.mean <- lapply(pos.mn, function(aa){
    apply(aa,c(1,2),mean)})
  obj$mn = array(NA, dim = c(TT, dimdat,numclust))
  for(iclust in 1:numclust){
    obj$mn[,,iclust] = t(post.mn.mean[[iclust]]) 
  }

  ## probabilities
  ## pos.gamma.mean = res$pos.gamma[,,last_draws_inds] %>% apply(c(1:2), mean)
  pos.SB <- apply(res$pos.gamma[,,last_draws_inds], c(2,3), function(ga)
    1/(1+exp(-t(ga) %*% rbind(1,datobj$X))))
  pos.MN <- apply(pos.SB, c(1,3), flowcut:::SB2MN)  
  post.pi.mean <- apply(pos.MN,c(1,2), mean)
  obj$prob = t(post.pi.mean)

  ## cluster covariances
  obj$sigma = array(NA, dim = c(numclust, dimdat, dimdat))
  obj$numclust = numclust
  res$pos.Sigma %>% dim()
  ## post.Sigma.mean <- apply(res$pos.Sigma,c(1,2,3), mean) %>% aperm(c(3,1,2))
  pos.Sigma.mean = res$pos.Sigma[,,,last_draws_inds] %>% apply(c(1,2,3), mean) %>% aperm(c(3,1,2)) 
  obj$sigma = pos.Sigma.mean

  ## mean regression coefficients
  pos.beta = res$pos.beta[,,,last_draws_inds]
  pos.beta.mean = pos.beta %>% apply(c(1,2,3), mean) 
  obj$beta = pos.beta.mean
  
  ## prob regression coefficients
  pos.gamma.mean = res$pos.gamma[,,last_draws_inds] %>% apply(c(1:2), mean)
  obj$alpha = pos.gamma.mean

  ## Bundle and return
  obj$numclust = numclust
  obj$TT = TT
  class(obj) = "flowmix"
  return(obj)
}
