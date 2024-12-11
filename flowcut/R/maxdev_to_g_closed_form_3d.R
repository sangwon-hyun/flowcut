# Generated from _main.Rmd: do not edit by hand

#' Closed form for obtaining the g parameter by simulation. Currently only works for *any* dimdat. 
#' 
#' @param X the design matrix, TT by p, (that doesn't include the intercept) 
#' @param dimdat the dimension of the cytogram space 
#' @param maxdev Maximum deviation of cluster means away from its grand mean. 
#' @param prior_prob prior probability
#' @param viz show the plot of the fitted relationship between the g parameter and the prior probability.
#' 
#' @return the g parameter with desired prior probability on maxdev 
#'
#' @export 
maxdev_to_g_closed_form_3d <- function(X, dimdat, maxdev,
                                        prior_prob = 0.99, viz = FALSE){

  ## stopifnot(dimdat == 3)
  obj = svd(X)
  pp = ncol(X)
  TT = nrow(X)
  ## maxdev = 0.115525

  ## chisq_cutoff = qchisq(p=1-(1-prior_prob^(1/numclust)), df = dimdat)
  chisq_cutoff = qchisq(p=1-(1-prior_prob), df = dimdat)
  
  ## Form the constant c and matrix P, multiplied to form the prior covariance cP.
  P = obj$u %*% t(obj$u)
  ## gg = (maxdev / chisq_cutoff / diag(P)) %>% median()
  gg = (maxdev^2 / chisq_cutoff / diag(P)) %>% median()
  A = gg * obj$v %*% diag(1/(obj$d^2)) %*% t(obj$v) 
  gP = (X %*% A %*% t(X)) ## prior covariance; not used now. 

  if(viz){
  
    ## Generate beta and form w = X beta_k
    nsim = 1000
    many_betas = matrixsampling::rmatrixnormal(n = nsim,
                                               M=matrix(0, pp, dimdat), U=A,
                                               V = diag(rep(1,dimdat)))
    wlist = lapply(1:nsim, function(ii)X %*% many_betas[,,ii])
    wlist_rownorms = wlist %>% lapply(function(w) apply(w, 1, function(a)sqrt(sum(a*a))))
    
    ## Make a plot of all |w|_2, w = X beta_k
    wlist_rownorms %>% do.call(rbind, .) %>% t() %>%  matplot(pch = 16, col = rgb(0,0,0,0.1),
                                                              cex = .5,
                                                                   ylab = "2norm of w[t]", xlab = "t")

    ## Draw the lines where the 95'th quantiles are.
    qs = wlist_rownorms %>% do.call(rbind, .) %>% t() %>% apply(1, function(one_w){
      quantile(one_w, probs=prior_prob)})
    lines(qs, col = 'red')
    abline(h = maxdev, lwd = 2, lty = 2, col = 'red')
    abline(h = median(qs), lwd = 2, lty = 2, col = 'yellow')
    legend("topleft", col = 'red', lwd=1, lty=1, legend = c(paste0(prior_prob," quantile at each t")),
           bty = "n")
    return(NULL)
  } else {
    return(gg)
  }
}
