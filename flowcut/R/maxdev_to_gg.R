# Generated from _main.Rmd: do not edit by hand

#' X is a (T x p) covariate matrix (that doesn't include the intercept)
#' 
#' @param X (T x p) matrix.
#' @param dimdat Dimension of cytogram data.
#' @param maxdev Maximum deviation of cluster means away from its grand mean.
#' @param viz Whether to produce a plot.
#' @param n.cores Number of cores.
#'
#' 
maxdev_to_gg <- function(X, dimdat, maxdev, viz = FALSE, n.cores = 1){

  ## Basic setup
  p = ncol(X) ## 39
  d = dimdat

  # Helper function
  # @param tX is the transpose of X
  ball.deviance <- function(gg, rr, tX, Nmc=5000,
                            nu0=d, nu1=p+1,
                            S0=diag(d), S1=diag(p+1),
                            simple= TRUE){
    inv.XTX  <- solve(tcrossprod(tX))
    ms.deviance <- rep(NA,Nmc) 

    Sig.ell <- matrixsampling::rinvwishart(Nmc,nu0+d,S0)
    beta.ell <- apply(Sig.ell,3,function(xx)
        as.matrix(matrixNormal::rmatnorm(M = matrix(0,d,p), 
                 U = xx, 
                 V = inv.XTX*gg, 
                 tol = .Machine$double.eps^0.5)),
        simplify = FALSE)
    xb <- lapply(beta.ell, function(bb) bb%*% tX)  
    mean.vec <- lapply(xb, function(bb)
        apply(bb,1,mean))
    ms.deviance <- mapply(function(xx,mm){
        mean(apply(xx,2,function(cols) crossprod(cols-mm)))
    }, xx=xb,mm=mean.vec)
    prob <- mean(ms.deviance > rr^2)
    if(simple){
      return(list(gg=gg,rr=rr,prob=prob))
    }else{
      return(list(gg = gg, rr = rr,
                  prob=prob,
                  msd = ms.deviance)) 
    }
  }

  gglist <- as.list((1:10)/10)##as.list(1:40/100)
  plist <- parallel::mclapply(gglist, function(gg){
    print('here')
      ball.deviance(gg, maxdev, t(X), Nmc = 1e4)$prob},
      mc.cores = min(n.cores, length(gglist)))

  ## Make linear interpolation at a fine grid
  res = stats::approx(x = unlist(gglist), y = unlist(plist), method="linear", n = 100)
  newx =  res$x
  newy = res$y

  ## Get closest point
  imin =  which.min(abs(newy - 0.05))

  # Make some plots to confirm
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

## hist(ball.deviance(0.01,0.5,t(X))$msd, breaks = "FD")
## hist(ball.deviance(0.1,0.5,X)$msd, breaks = "FD")
## hist(ball.deviance(0.2,0.5,X)$msd, breaks = "FD")
## hist(ball.deviance(0.5,0.5,X)$msd, breaks = "FD")
## hist(ball.deviance(1,0.5,X)$msd, breaks = "FD")
## hist(ball.deviance(2,0.5,X)$msd, breaks = "FD")
