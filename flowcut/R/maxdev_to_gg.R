# Generated from _main.Rmd: do not edit by hand

#' Obtain the g parameter by simulation. 
#' 
#' @param X the design matrix, TT by p, (that doesn't include the intercept) 
#' @param dimdat the dimension of the cytogram space 
#' @param maxdev Maximum deviation of cluster means away from its grand mean. 
#' @param numclust the number of experts 
#' @param ggvec the vector of g parameter values to calculate the prior probability by Monte Carlo samples  
#' @param Nmc Monte Carlo simulation sample size, with default value being 1e4. 
#' @param n.cores the number of CPU cores to be used for parallelization 
#' @param viz show the plot of the fitted relationship between the g parameter and the prior probability.
#'
#' 
#' @return the g parameter with desired prior probability on maxdev 
#'
#' @export 
maxdev_to_gg <- function(X, dimdat, maxdev, numclust, ggvec,
                           Nmc = 1e4 , prior.prob = 0.99, 
                         viz = FALSE, n.cores = 1, verbose = FALSE){
    ## X is T x p 
    ## Basic setup
    p = ncol(X) ## 39
    d = dimdat
    
  # Helper function
  # @param tX is the transpose of X
    ball.deviance <- function(gg, rr, tX, 
                              Nmc=5000,
                              nu0=d,
                              nu1=p+1,
                              S0=diag(d),
                              S1=diag(p+1)){
        inv.XTX  <- solve(tcrossprod(tX))
        max.deviance <- rep(NA,Nmc) 
        
        Sig.ell <- matrixsampling::rinvwishart(Nmc,nu0+d,S0)
        beta.ell <- apply(Sig.ell,3,function(xx)
            as.matrix(matrixNormal::rmatnorm(M = matrix(0,d,p), 
                                             U = xx, 
                                             V = inv.XTX*gg, 
                                             tol = .Machine$double.eps^0.5)),
            simplify = FALSE)
        xb <- lapply(beta.ell, function(bb) bb%*% tX) 
                                        # length of list: Nmc 
                                        # each element: d x T 
        max.deviance <- lapply(xb, function(mat){
            apply(mat,2,function(cols) crossprod(cols)) %>% max()
        }) %>% unlist()
        prob <- mean(max.deviance <= rr^2)  
        return(list(gg = gg, rr = rr, prob = prob))
    }
    
    gglist <- as.list(ggvec)##as.list(1:40/100)
    n.cores <- min(n.cores, length(gglist))
    plist <- parallel::mclapply(1:length(gglist), function(igg){
      if(verbose) printprogress(igg, length(gglist))
      gg = gglist[[igg]]
      ball.deviance(gg, maxdev, t(X), Nmc = Nmc)$prob
    },
    mc.cores = n.cores)

    
    ## Make linear interpolation at a fine grid
    res = stats::approx(x = unlist(gglist), y = unlist(plist), method="linear", n = 100)
    newx =  res$x
    newy = res$y
    
    ## Get closest point
    imin =  which.min(abs(newy - prior.prob ^(1/numclust))) 
    
  ## Make some plots to confirm
    if(viz){
        plot(gglist,plist,type = "l",
             ylab = "Prob outside of \n radius r", xlab = "Value of |gg|")
        graphics::abline(h=0.05,lwd=2,col = "red")
        graphics::abline(h=0.01,lwd=2,col = "red")
        ## lines(y=newy, x=newx,lwd=2,col="blue")
        ## abline(v=newx[imin])
    }
    
    return(newx[imin])
  }


