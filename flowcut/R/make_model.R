# Generated from _main.Rmd: do not edit by hand

#' From an original set of model parameters (|true_model|),
#' generate synthetic 2-cluster 1-dimensional data with equal probabilities.
#'
#' @param isignal 0 to 10, which generates the means.
#' @param orig_model Original model of class flowmix or flowcut; a list that contains alpha, beta and TT.
#'
#' @return A list with beta, mn, alpha, prob, X, sigma, TT, numclust.
#' @export
make_model <- function(orig_model, isignal){

  ## Setup
  stopifnot(isignal %in% 0:10)
  new_model = orig_model
  new_model$numclust = 2

  ## Calculate + renormalize the probabilities
  link = cbind(1, orig_model$X) %*% t(orig_model$alpha)
  new_model$prob = exp(link) / rowSums(exp(link))
  new_model$prob  %>% matplot(type = 'l', lty = 1)

  ## ## Probabilities
  ## alphamat = orig_model$alpha
  ## alphamat[,-1] = 0
  ## alphamat[,1] = 1
  ## new_model$alpha = alphamat
  ## ## new_model$prob = matrix(1/2, nrow = orig_model$TT, ncol = 2)
  ## cbind(1, new_model$X) %*% alphamat

  ## Take the two intercepts
  intp_high = orig_model$beta %>% .[[1]]%>% .["intp",]  
  intp_low = orig_model$beta %>% .[[2]]%>% .["intp",] 
  increment = (intp_high - intp_low)/10
  
  ## Bring the larger mean down.
  new_model$beta[[1]]["intp",] = intp_low + increment * isignal
  new_model$mn = array(NA, dim = c(orig_model$TT, 1, 2))
  new_model$mn[,,1] = (cbind(1,new_model$X))  %*%  (new_model$beta[[1]])
  new_model$mn[,,2] = (cbind(1,new_model$X))  %*%  (new_model$beta[[2]])
  
  ## Optional: plot the means
  if(FALSE){
    new_model$mn[,1,] %>% matplot(type = 'l', lty = 1)
    new_model$prob %>% matplot(type = 'l', lty = 1)
  }
  
  return(new_model)
}
