# Generated from _main.Rmd: do not edit by hand

#' Closed form for obtaining the g parameter by simulation. Currently only works for dimdat=1
#' 
#' @param X the design matrix, TT by p, (that doesn't include the intercept) 
#' @param dimdat the dimension of the cytogram space 
#' @param maxdev Maximum deviation of cluster means away from its grand mean. 
#' @param prior_prob prior probability
#'
#' 
#' @return the g parameter with desired prior probability on maxdev 
#'
#' @export 
maxdev_to_g_closed_form_1d <- function(X, dimdat, maxdev, prior_prob = 0.99){

  stopifnot(dimdat == 1)
  obj = svd(X)
  pp = ncol(X)
  TT = nrow(X)
  ## maxdev = 0.115525

  ## z_cutoff = 2.58 ## for 99% for one cluster
  ## z_cutoff = qnorm(1-(1-prior_prob^(1/numclust))/2)
  z_cutoff = qnorm(1-(1-prior_prob)/2)
  
  ## Form the constant c and matrix P, multiplied to form the prior covariance cP.
  P = obj$u %*% t(obj$u)
  c = (maxdev^2 / z_cutoff^2 / diag(P)) %>% median
  A = c * obj$v %*% diag(1/(obj$d^2)) %*% t(obj$v) 
  cP = (X %*% A %*% t(X)) ##%>% diag() %>% plot()

  return(c)
}
