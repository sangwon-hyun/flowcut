# Generated from _main.Rmd: do not edit by hand

#' Imputing particles to impute the y's based on a model defined by the posterior
#' mean of the last 500 draws. This is separated out from the Gibbs since the result takes up a lot of memory.
#' (This is all code in the Gibbs.fast() function.)
#'
#' 
#' @param ylist Original particles
#' @param countslist Biomass list
#' @param Cbox Bounding box.
#' @param mcres The converted MCMC result, using e.g., |mcmc_res_to_flowmix(res,
#'   last_draws_inds = 1001:1500)|
#' @param n.cores number of CPU cores to use
#' 
#' @export
impute <- function(ylist, countslist, Cbox, mcres, n.cores=1){

  ## if(FALSE){
  ##   ## Temporary
  ##   ylist = datobj$ylist
  ##   countslist = datobj$countslist

  ##   ## TEMPORARY: Add a bounding box |Cbox|
  ##   bounds.lower <- lapply(datobj$ylist,function(xx)
  ##     Rfast::colMins(xx,value = TRUE)) %>%
  ##     do.call(rbind,.) %>% Rfast::colMins(.,value = TRUE)
  ##   bounds.upper <- lapply(datobj$ylist,function(xx)
  ##     Rfast::colMaxs(xx,value = TRUE)) %>%
  ##     do.call(rbind,.) %>% Rfast::colMaxs(.,value = TRUE)
  ##   Cbox <- cbind(bounds.lower,bounds.upper)
  ##   Cbox <- res$dat.info$Cbox
  ##   datobj$Cbox = Cbox

  ##   ## Temporary
  ##   mcres = mcmc_res_to_flowmix(res, last_draws_inds = 1001:1500)
  ## }
  
  ## Get everything
  ntlist = sapply(ylist, nrow)
  TT = length(ylist)
  NN <- sum(ntlist)
  dimdat <- ncol(ylist[[1]])
  countsTotal <- sapply(countslist,sum) %>% sum() 
  alpha.factor <- NN/countsTotal  
  W.list <- lapply(countslist, function(xx) xx*alpha.factor) 

  ## More processing
  Censor.list <- lapply(ylist, function(x) censorIndicator(x,Cbox))  
  censor.01.list <- lapply(Censor.list, function(x) apply(x,1, function(xx){
    sum(abs(xx)) > 0})) ## expensive
  censor.which.list <- lapply(censor.01.list, function(cc) which(cc==TRUE))
  censored.ylist <- mapply(function(yy,c01) {yy[c01==TRUE,,drop=FALSE]},
                                       yy = ylist, c01 = censor.01.list)
  censored.C.list <- mapply(function(c01,cc){cc[c01==TRUE,,drop=FALSE]},
                                        c01 = censor.01.list, cc=Censor.list)
  censored.W.list <- mapply(function(c01,ww){ww[c01==TRUE]},
                                        c01 = censor.01.list, w=W.list)
  ntlist.censor <- sapply(censor.01.list, sum) 

  Z.list = flowmix::gate(mcres, ylist, NULL, 1111, eps_estep = 1E-20)

  censored.Z.list <- mapply(function(c01,zz){zz[c01 == TRUE]},
                            c01 = censor.01.list, zz = Z.list)

  ## Sort of expensive
  samp.region.list <- lapply(1:TT, function(tt){
      c01 <- censor.01.list[[tt]]
      cc <- Censor.list[[tt]]
      if(sum(c01==TRUE)>1){
          apply(cc[c01==TRUE,,drop=FALSE],1,function(xx)
              ## flowcut:::sample.region(xx,Cbox))
              sample.region(xx,Cbox))
      }else if(sum(c01==TRUE)==1){
          matrix(sample.region(cc[c01==TRUE,,drop=FALSE],Cbox),ncol=1)
      }else{
          NULL
      }
  })

  ## Get the censored particles
  Sig.ell = mcres$sigma %>% aperm(c(2,3,1))
  imputed.ylist <- parallel::mclapply(1:TT, function(tt) {
    mu.mat = mcres$mn[tt,,]
    yy = impute.censored(ww = censored.W.list[[tt]], 
                         yy = censored.ylist[[tt]],
                         zz = censored.Z.list[[tt]],
                         cc.info.mat =  censored.C.list[[tt]],
                         bounds.mat = samp.region.list[[tt]],
                         ##mu.list[[tt]]
                         mu.mat = mu.mat, Sigma.ell = Sig.ell,
                         dimdat = dimdat)
    if(dimdat==1 & !is.null(yy)) yy = t(yy)
    return(yy)
  }, mc.cores = n.cores, mc.preschedule = FALSE)             


  ## Finally, ca
  imputed_ylist <-
    parallel::mcmapply(function(yy, c01, imp){replace(yy, c01, imp)}, 
                       yy = ylist, c01 = censor.01.list, imp = imputed.ylist,
                       mc.cores = n.cores, SIMPLIFY = FALSE)            
  return(list(imputed_ylist = imputed_ylist, censor.01.list = censor.01.list,
              memlist = Z.list))
}
