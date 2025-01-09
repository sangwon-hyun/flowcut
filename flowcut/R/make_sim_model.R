# Generated from _main.Rmd: do not edit by hand

#' From an original set of model parameters (|true_model|),
#' generate synthetic 2-cluster 1-dimensional data with equal probabilities.
#'
#' @param isignal 0 to 10, which generates the means.
#' @param orig_model Original model of class flowmix or flowcut; a list that
#'   contains alpha, beta and TT.
#' @param shrink_alpha If TRUE, "shrink" the alpha coefficients to 40% their
#'   size. Defaults to FALSE.
#' @param coef_is_dense If TRUE, add some Gaussian noise on the non-intercept
#'   coefficients to make them dense (instead of sparse).
#'
#' @return A list with beta, mn, alpha, prob, X, sigma, TT, numclust.
#' @export
make_sim_model <- function(orig_model, isignal, shrink_alpha = FALSE, coef_is_dense = FALSE){

  ## Setup
  stopifnot(isignal %in% 0:10)
  new_model = orig_model
  new_model$numclust = 2
  ## stopifnot(class(orig_model) %in% c("flowmix")) ##, "flowcut"

  ## If dense model is called for, make the coefficients all nonzero
  if(coef_is_dense){
    dense_model = new_model
    beta1 = dense_model$beta[[1]][-1]
    beta2 = dense_model$beta[[2]][-1]
    sd1 = sd(beta1[which(beta1!=0)])/3
    sd2 = sd(beta2[which(beta2!=0)])/3
    p = nrow(dense_model$beta[[1]]) - 1
    set.seed(4891)
    dense_model$beta[[1]][-1] = dense_model$beta[[1]][-1] + rnorm(n = p, mean=0, sd=sd1)
    dense_model$beta[[2]][-1] = dense_model$beta[[2]][-1] + rnorm(n = p, mean=0, sd=sd2)
    alpha = dense_model$alpha[,-1]
    sd_alpha = sd(alpha[which(alpha!=0)])/3
    dense_model$alpha[,-1] = dense_model$alpha[,-1] + rnorm(n=2*p, mean=0, sd = sd_alpha)
    new_model = dense_model
  }

  ## Calculate + renormalize the probabilities
  link = cbind(1, orig_model$X) %*% t(orig_model$alpha)
  new_model$prob = exp(link) / rowSums(exp(link))
  new_model$prob  %>% matplot(type = 'l', lty = 1)

  ## Take the two intercepts
  intp_high = orig_model$beta %>% .[[1]]%>% .["intp",]  
  intp_low = orig_model$beta %>% .[[2]]%>% .["intp",] 
  increment = (intp_high - intp_low)/10
  
  ## Bring the larger mean down.
  new_model$beta[[1]]["intp",] = intp_low + increment * isignal
  new_model$mn = array(NA, dim = c(orig_model$TT, 1, 2))
  new_model$mn[,,1] = (cbind(1,new_model$X))  %*%  (new_model$beta[[1]])
  new_model$mn[,,2] = (cbind(1,new_model$X))  %*%  (new_model$beta[[2]])


  ## Shrink the probabilities to be closer to each other.
  if(shrink_alpha){
    Xta1 = (new_model$X) %*% ((new_model$alpha[1,-1]) * 0.40)
    Xta2 = (new_model$X) %*% ((new_model$alpha[2,-1]) * 0.40)
    ##lines(exp(Xta1)/(exp(Xta1) + exp(Xta2)), ylim = c(0,1))
    ##lines(exp(Xta2)/(exp(Xta1) + exp(Xta2)), col = 'red')
    new_model$prob[,1] = exp(Xta1)/(exp(Xta1) + exp(Xta2))
    new_model$prob[,2] = exp(Xta2)/(exp(Xta1) + exp(Xta2))
  }

  return(new_model)
}
