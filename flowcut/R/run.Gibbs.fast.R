# Generated from _main.Rmd: do not edit by hand

#' Runs a gibbs sampler to estimate the posterior distribution of the flowcut 
#' model.
#'
#' @param ylist data.
#' @param countlist weights (counts) for data in |ylist|.
#' @param X covariates; a (p x T) matrix (TODO: change code so that X is T x p).
#' @param Nmc Number of MCMC iterations
#' @param Nburn Number of burn-in iterations.
#' @param Cbox Censored box.
#' @param user.prior User-supplied prior. Otherwise, prior defaults to ___.
#' @param gg Size of Normal prior.
#' @param verbose Whether to be loud.
#' @param warm.start If supplied, restart the MCMC at these values.
#' @param tol NOT USED.
#' @param n_cores Number of cores for multiple cores.
#' 
#'
#' @return 
run.Gibbs.fast <- function(ylist, countslist, X, 
                           numclust,
                           Nmc = 3000,
                           Nburn = 500,
                           Cbox = NULL,  
                           user.prior = NULL,
                           gg=1,
                           verbose = FALSE,
                           warm.start = NULL,
                           tol = 1/1e8,
                           n.cores = 1){

  ## Basic setup
  TT <- length(ylist)
  p <- dim(X)[1]
  dimdat = ncol(ylist[[1]])
  ntlist = sapply(ylist, nrow)
  NN <- sum(ntlist)
  tt.impute <- min(20,floor(Nburn/5))
  ylist.raw <- ylist ## store raw ylist 
  ## n.cores <- detectCores()
  n.cores = min(n.cores, TT)

  ##### pre computed quantities
  X.list <- as.list(as.data.frame(X))
  Xp <- rbind(1,X) 
  Xp.list <- as.list(as.data.frame(Xp))
  W.sq.list <- lapply(countslist, function(xx) sqrt(xx))
  mt <- lapply(countslist, sum)%>%unlist()
  MM <- sum(mt) 
      
  XtXtT <- lapply(X.list,function(x) x%*%t(x))
  XtXtTp <- lapply(Xp.list, function(x) x%*%t(x))
  ggXtXtTp <- lapply(Xp.list, function(x) x%*%t(x)/gg)
  XTX <- X%*%t(X)
  XTXp <- Xp%*%t(Xp)
  inv.XTX <- Rfast::spdinv(XTX)
  inv.XTXp <- Rfast::spdinv(XTXp)
  ##    Y.grand.mean <- Rfast::colmeans(do.call(rbind,ylist))

  ## Build censored box
  if(!is.null(Cbox)){
      Censor.list <- mclapply(ylist, function(x) censorIndicator(x,Cbox),
                              mc.cores = n.cores)
      censor.01.list <- mclapply(Censor.list, function(x)
          apply(x,1, function(xx)
              sum(abs(xx))>0), mc.cores = n.cores)
      censored.ylist <- mcmapply(function(yy,c01) {yy[c01==TRUE,,drop=FALSE]},
                                  yy = ylist, c01 = censor.01.list,
                                  mc.cores = n.cores)
      censored.C.list <- mcmapply(function(c01,cc){cc[c01==TRUE,,drop=FALSE]},
                                  c01 = censor.01.list, cc=Censor.list,
                                  mc.cores = n.cores)
      censored.countslist <- mcmapply(function(c01,ww){ww[c01==TRUE]},
                                  c01 = censor.01.list, w=countslist,
                                  mc.cores = n.cores)  
      ntlist.censor <- sapply(1:TT,function(tt) sum(censor.01.list[[tt]])) 
      samp.region.list <- lapply(1:TT, function(tt){
          c01 <- censor.01.list[[tt]]
          cc <- Censor.list[[tt]]
          if(verbose){
              print(paste("Time=",tt,", censored obs.= ",sum(c01==TRUE),
                          ", censored ratio =", round(sum(c01==TRUE)/ntlist[tt],2), sep=" ")) 
          }
          if(sum(c01==TRUE)>1){
              apply(cc[c01==TRUE,,drop=FALSE],1,function(xx) sample.region(xx,Cbox))
          }else if(sum(c01==TRUE)==1){
              matrix(sample.region(cc[c01==TRUE,,drop=FALSE],Cbox),ncol=1)
          }else{
              NULL
          }}) 
  }

##### prior specifications
    if(is.null(user.prior)){
        nu0=dimdat
        nu1=p+1 
        S0=diag(dimdat)
        S1=diag(p+1)
        ## inv.Omega <- solve(S1)
    }

##### initialize
    if(!is.null(warm.start)){
        print("continue with a previous draw") 
        Z.list <- warm.start$Z.list
        ylist <- warm.start$ylist
        ## Omega <- warm.start$Omega 
        beta.ell <-  warm.start$beta 
        gamma.ell <- warm.start$gamma
        Sig.ell <- warm.start$Sigma  
        
        Nburn <- 0 ## no need of burn-in 
        tt.impute <- min(25,floor(Nburn/5)) 
    }else{
        Z.list <- lapply(1:TT,function(tt) sample(1:numclust,ntlist[tt],replace = TRUE)) 
        ## Omega <- rinvwishart(1,nu1+p+1,S1)[,,1]
        beta.ell <-  array(rnorm(numclust*(p+1)*dimdat), c(dimdat,p+1,numclust))
        gamma.ell <- matrix(rnorm((p+1)*(numclust-1)),nrow=p+1,ncol=numclust-1) 
        Sig.ell <- rinvwishart(numclust,nu0+dimdat,S0)## %>% as.matrix()
    }

##    Omega.ell <- rinvwishart(numclust-1,nu1+p+1,S1)
##    beta.ell <-  array(rnorm(numclust*p*dimdat), c(dimdat,p,numclust))
##    beta0.ell <- matrix(0,nrow=dimdat,ncol=numclust)

    SX.ell <- array(0, c(p,p,numclust))
    Sy.ell <- matrix(0,nrow = dimdat, ncol = TT)
    Sxy.ell <- array(0, c(p,dimdat,numclust))
    
##### store
    ## burn.beta0 <- array(0,c(dimdat,numclust,Nburn))
    ## burn.beta <- array(0,c(dimdat,p,numclust,Nburn))
    burn.beta <- array(0,c(dimdat,p+1,numclust,Nburn))
    burn.Sigma <- array(0,c(dimdat,dimdat,numclust,Nburn)) 
    burn.gamma <- array(0,c(p+1,numclust-1,Nburn))
    ## burn.Omega <- array(0,c(p+1,p+1,Nburn))
    burn.avgloglik <- rep(NA, Nburn)
    
    ## pos.beta0 <- array(0,c(dimdat,numclust,Nmc))
    ##pos.beta <- array(0,c(dimdat,p,numclust,Nmc))
    pos.beta <- array(0,c(dimdat,p+1,numclust,Nmc))
    pos.Sigma <- array(0,c(dimdat,dimdat,numclust,Nmc)) 
    pos.gamma <- array(0,c(p+1,numclust-1,Nmc))
    ##pos.Omega <- array(0,c(p+1,p+1,Nmc))
    pos.avgloglik <- rep(NA, Nmc) 
    ## if(is.null(Cbox)==FALSE){
    ##     pos.imputedimdat.Y <- array(0,c(dimdat,sum(nt.censor),Nmc))
    ## } 

    if(verbose){
        print("The Gibbs sampler starts now.")
        list.iter <- c(1, Nburn+1, 1:floor((Nmc+ Nburn-1)/100)*100, Nmc + Nburn)
        ptm <- NULL 
    }
    
##### posterior sampling starts here
    for(jj in 1:(Nburn+Nmc)){
        cat("MCMC iteration:", jj, fill = TRUE) 
        if(verbose){
            if(jj ==1 & Nburn > 0){
                print("Burn-in period starts.")
                print(Sys.time())                
            }
            if(jj == tt.impute + 1){
                ptm <- proc.time()
            }
            if(jj ==Nburn+1 & Nburn > 0){
                ## burn.in.last.draw <- list(Z.list=Z.list,ylist=ylist)
                ## save(burn.in.last.draw, file = "Burn-in-last-draw.Rdata")
                ## try(dev.off(),silent = TRUE)
                ## par(mfrow = c(dimdat,numclust))
                ## par(mar = c(numclust,dimdat,numclust,dimdat)/1.5)
                ## try(for(dd in 1:dimdat){
                ##         for(kk in 1:numclust){
                ##             title <- paste("beta0 traceplot", "dim = ", dd, ", clust = ", kk, sep = " ")
                ##             plot(floor(Nburn/2):Nburn,  
                ##                  burn.beta0[dd,kk,floor(Nburn/2):Nburn],
                ##                  type="l", main=title, lwd = 2)
                ##         }    
                ##     },silent = TRUE)
                plot(burn.avgloglik[(tt.impute+1):Nburn],
                     type="l",lwd=2)
                print(paste("Burn-in sample size: ", Nburn, sep = " "))
                print("Burn-in period ends.")
                print("Burn-in time cost per iteration, in seconds")
                burn.tc <- (proc.time()-ptm)/(Nburn - tt.impute) 
                print(round(burn.tc,2))

                print(Sys.time())
                
                print("Collecting posterior samples......")
                print(paste("Expected time to draw", Nmc, "posterior samples: ", 
                            round(burn.tc[3]*Nmc/60/60,2), " hours", sep=" "))
            }
            if(jj %in% list.iter & jj > Nburn){
                progress(jj - Nburn,Nmc)
            }
            if(jj %in% list.iter & jj < Nburn+1){
                progress(jj,Nburn)
            }
        }

################################################ 
        ## experts' estimation ###
################################################ 
        
        mt.ell <- mcmapply(function(ww,zz){
            sapply(1:numclust, function(kk) sum(ww[zz==kk]))},
            ww = countslist, zz = Z.list, SIMPLIFY = TRUE,
            mc.cores = n.cores)
        m.ell <-  Rfast::rowsums(mt.ell)

        ## nt.ell  <- do.call(cbind,mclapply(Z.list, function(zz){
        ##     sapply(1:numclust, function(kk) as.numeric(sum(zz == kk)))},
        ##     mc.cores = min(n.cores, T))) 
        ## n.ell <- Rfast::rowsums(nt.ell)
        
        for(ell in 1:numclust){ 
            ## mm0t <- do.call(rbind,mapply(function(xx,yy,zz){
            ##     Rfast::colsums(Rfast::eachrow(as.matrix(yy[zz==ell,]), beta.ell[,,ell] %*% xx,"-"))},
            ##     xx=X.list, yy = ylist, zz = Z.list, SIMPLIFY = FALSE))
            ## mm0 <- Rfast::colsums(mm0t)/n.ell[ell] 
            ## beta0.ell[,ell] <- mvrnorm(1, mm0, Sig.ell[,,ell]/n.ell[ell])

            ## SX.ell <- XTXp/gg + Reduce('+',Map(`*`, XtXtTp, nt.ell[ell,]))  
            ## inv.SX.ell <- Rfast::spdinv(SX.ell) ## (p+1) x (p+1)
            ## Sy.ell <- do.call(rbind, mcmapply(function(yy,zz){
            ##     Rfast::colsums(as.matrix(yy[zz==ell,]))},
            ##     yy=ylist,zz=Z.list, SIMPLIFY = FALSE,
            ##     mc.cores = min(n.cores, T)))
            

            SX.ell <- XTXp/gg + Reduce('+',Map(`*`, XtXtTp, mt.ell[ell,]))  
            inv.SX.ell <- Rfast::spdinv(SX.ell) ## (p+1) x (p+1)
            Sy.ell <- mapply(function(ww, yy, zz){
                t(as.matrix(ww[zz==ell])) %*% yy[zz==ell,,drop=FALSE]},
                ww = countslist, yy = ylist, zz = Z.list, SIMPLIFY = TRUE)
          if(dimdat == 1) Sy.ell = rbind(Sy.ell)
          ## if(!is.numeric(Sy.ell[1])) browser()
            Sxy.ell <- Xp %*% t(Sy.ell)  ## (p+1) x d
            beta.ell[,,ell] <- rmatnorm(M = t(Sxy.ell) %*% inv.SX.ell, 
                                        U = Sig.ell[,,ell], 
                                        V = inv.SX.ell, 
                                        tol = .Machine$double.eps^0.95)
            
            sse <- mcmapply(function(ww,xx,yy,zz) {
                wcrossprod.fast(Rfast::eachrow(as.matrix(yy[zz==ell,,drop=FALSE]),
                                               beta.ell[,,ell]%*%xx, '-'),
                                ww[zz==ell], weighting = TRUE)},
                xx = Xp.list, yy = ylist, zz=Z.list, ww = W.sq.list, 
                SIMPLIFY = FALSE,
                mc.cores = n.cores) 
          if(dimdat == 1) Sn.ell = sum(unlist(sse))
          if(dimdat > 1) Sn.ell <- Reduce('+',sse)
            Sig.ell[,,ell] <- rinvwishart(1,nu0+dimdat + m.ell[ell], S0 + Sn.ell)[,,1]
        }
        
        ################################################ 
        ## expert assignment ###
################################################
        ## mt.ell <- mcmapply(function(ww,zz){
        ##     sapply(1:numclust, function(kk) sum(ww[zz==kk]))},
        ##     ww = countslist, zz = Z.list, SIMPLIFY = TRUE,
        ##     mc.cores = n.cores)
        ## m.ell <-  Rfast::rowsums(mt.ell)

        ## XpGamma.abs <- abs(t(gamma.ell) %*% Xp) ## numclust-1 x T 
        ## mt.cumsum <- Rfast::colCumSums(mt.ell) 
        ## Mt.ell <- rbind(nt,-sweep(mt.cumsum[-numclust,], 2, mt.cumsum[numclust,])) 
        ## omega.tell <- matrix(mcmapply(pgdraw, round(Mt.ell[-numclust,]), XpGamma.abs,
        ##                               mc.cores = n.cores),
        ##                      nrow=numclust-1, ncol = TT)
        ## kappa <- mt.ell[-numclust,] - Mt.ell[-numclust,]/2


        XpGamma.abs <- abs(t(gamma.ell) %*% Xp) ## numclust-1 x TT  
        ## nt.cumsum <- Rfast::colCumSums(as.matrix(nt.ell))  
        ## Nt.ell <- rbind(nt,-sweep(nt.cumsum[-numclust,], 2, nt.cumsum[numclust,])) 
        ## omega.tell <- matrix(mcmapply(pgdraw, Nt.ell[-numclust,], XpGamma.abs,
        ##                               mc.cores = n.cores),
        ##                      nrow=numclust-1, ncol = TT)
        ## kappa <- nt.ell[-numclust,] - Nt.ell[-numclust,]/2

        mt.cumsum <- Rfast::colCumSums(as.matrix(mt.ell))  
        Mt.ell <- rbind(mt,-sweep(mt.cumsum[-numclust,,drop=FALSE], 2, mt.cumsum[numclust,,drop=FALSE])) 
        omega.tell <- matrix(mcmapply(pgdraw, round(Mt.ell[-numclust,]), XpGamma.abs,
                                      mc.cores = n.cores),
                             nrow=numclust-1, ncol = TT)
        kappa <- mt.ell[-numclust,,drop=FALSE] - Mt.ell[-numclust,,drop=FALSE]/2

        ##        inv.Omega <- Rfast::spdinv(Omega)
        for(ell in 1:(numclust-1)){
            V.omell <- Rfast::spdinv(Reduce('+', Map('*',XtXtTp, omega.tell[ell,]))) ##+diag(p+1)/TT ) Mod2 
            m.omell <- V.omell%*% Reduce('+', Map('*',Xp.list,kappa[ell,]))
            gamma.ell[,ell] <- Rfast::rmvnorm(1,m.omell, V.omell)
        }        


        XpGamma <- t(gamma.ell) %*% Xp ## K-1 x TT 
        pi.sb <- 1/(1+exp(-XpGamma))
        pi.mn <- apply(pi.sb,2,SB2MN)
        logpi.list <- mclapply(1:TT,function(t) log(pi.mn[,t]),
                               mc.cores = min(n.cores,TT)) 
        
        chol.Sig.ell <- apply(Sig.ell,3, chol)
        if(dimdat == 1) chol.Sig.ell = rbind(chol.Sig.ell)
        chol.Sig.list <- lapply(1:numclust,function(kk) matrix(chol.Sig.ell[,kk],nrow = dimdat))
        mu.list <- lapply(1:TT, function(tt){
          one_mu = sapply(1:numclust, function(kk){
            beta.ell[,,kk,drop=FALSE] %*% Xp[,tt,drop=FALSE]
          })
          if(dimdat==1) one_mu = rbind(one_mu)
          return(one_mu)
        })
#            beta0.ell[,kk] + beta.ell[,,kk] %*% X[,tt])) 
        
        logPiZ <- mcmapply(function(ww, xx, yy, mm, pp){
            sapply(1:numclust, function(kk)
                mvnfast::dmvn(yy, mm[,kk], chol.Sig.list[[kk]],
                              log=TRUE, isChol = TRUE) + pp[kk])}, ## mod here 
            ww = countslist, xx =X.list, yy = ylist, mm = mu.list, pp=logpi.list,
            mc.cores = n.cores, SIMPLIFY = FALSE)

        Z.list <- mclapply(logPiZ, function(pp) apply(pp,1, function(lpi)
            sample(1:numclust,1,prob=softmax(lpi))),
            mc.cores = n.cores)

        ## mt.ell <- mcmapply(function(ww,zz){
        ##     sapply(1:numclust, function(kk) sum(ww[zz==kk]))},
        ##     ww = countslist, zz = Z.list, SIMPLIFY = TRUE,
        ##     mc.cores = n.cores)
        ## m.ell <-  Rfast::rowsums(mt.ell)


        ## print(sort(round(n.ell/NN,3)))
        ## print(min(n.ell))

################################################ 
            ## censored data imputation
################################################
         if((jj> tt.impute | !is.null(warm.start)) & !is.null(Cbox)){
##        if(is.null(Cbox)==FALSE){
            censored.Z.list <- mcmapply(function(c01,zz){zz[c01==TRUE]},
                                        c01 = censor.01.list, zz=Z.list,
                                        mc.cores = n.cores)  
            imputed.ylist <- mclapply(1:TT, function(tt) {
            ## imputed.ylist <- lapply(1:TT, function(tt){
              yy =
                impute.censored(ww = censored.countslist[[tt]], 
                                yy = censored.ylist[[tt]],
                                zz = censored.Z.list[[tt]],
                                cc.info.mat =  censored.C.list[[tt]],
                                bounds.mat = samp.region.list[[tt]],
                                mu.mat = mu.list[[tt]], Sigma.ell = Sig.ell,
                                dimdat = dimdat)
            ## })             
            },mc.cores = n.cores, mc.preschedule = FALSE)             
            for(tt in 1:TT){
                ylist[[tt]][censor.01.list[[tt]]==TRUE,] <- imputed.ylist[[tt]]
            }            
        }

        loglik <- loglik_eval(mu.list, chol.Sig.list, 
                              countslist, X.list, ylist, Z.list,
                              as.list(ntlist),
                              simple = TRUE) 
        ## print(sort(round(n.ell/NN,3)))
        ## print(sort(round(m.ell/MM,3)))
        if(jj %% 10 == 0) print(paste("avg loglikelihood: ", round(loglik,2), sep=" "))


################################################ 
            ## collecting posterior samples 
################################################ 
        
        if(jj >Nburn){
###            pos.beta0[,,jj-Nburn] <- beta0.ell
            pos.beta[,,,jj-Nburn] <- beta.ell
            pos.Sigma[,,,jj-Nburn] <- Sig.ell
            pos.gamma[,,jj-Nburn] <- gamma.ell 
##            pos.Omega[,,,jj-Nburn] <- Omega.ell
##            pos.Omega[,,jj-Nburn] <- Omega
            pos.avgloglik[jj-Nburn] <- loglik
            ## if(is.null(Cbox)==FALSE){
            ##     pos.imputed.Y[,,jj-Nburn] <- matrix(unlist(imputed.ylist),nrow = d) 
            ## }
        } else{
###            burn.beta0[,,jj] <- beta0.ell
            burn.beta[,,,jj] <- beta.ell
            burn.Sigma[,,,jj] <- Sig.ell
            burn.gamma[,,jj] <- gamma.ell 
##            burn.Omega[,,,jj] <- Omega.ell   
##            burn.Omega[,,jj] <- Omega 
            burn.avgloglik[jj] <- loglik
        }    
    }

    if(verbose ){
        print("All work is done.")
        print(paste("Stored MCMC sample size: ", Nmc, sep = " ")) 
        print(Sys.time())                
        print("total time cost, in min")
        print(round((proc.time()-ptm)/60,1))
        print("total time cost, in hours")
        print(round((proc.time()-ptm)/60/60,2))
        print("time cost per iteration, in seconds") 
        print((proc.time()-ptm)/(Nburn+Nmc - tt.impute))
    }

    if(is.null(Cbox)){
        ret <- list( ## burn.beta0=burn.beta0,
            burn.beta=burn.beta,burn.Sigma=burn.Sigma,
            burn.gamma=burn.gamma, ## burn.Omega=burn.Omega,
            ## pos.beta0=pos.beta0,
            pos.beta=pos.beta,pos.Sigma=pos.Sigma,
            pos.gamma=pos.gamma, ## pos.Omega=pos.Omega,
            burn.avgloglik=burn.avgloglik,
            pos.avgloglik=pos.avgloglik)
    }else{
        ret <- list( ## burn.beta0 = burn.beta0,
            burn.beta = burn.beta,
            burn.Sigma = burn.Sigma,
            burn.gamma = burn.gamma,
            ## burn.Omega = burn.Omega,
            ## pos.beta0 = pos.beta0,
            pos.beta = pos.beta,
            pos.Sigma = pos.Sigma,
            pos.gamma = pos.gamma,
            ## pos.Omega = pos.Omega,
            pos.imputed.Y = do.call(rbind,imputed.ylist),
            burn.avgloglik = burn.avgloglik,
            pos.avgloglik = pos.avgloglik,
            raw.ylist = ylist.raw, 
            last.ylist = ylist,
            last.Z.list = Z.list,
            last.gamma = gamma.ell,
            ## last.Omega = Omega,
            last.beta = beta.ell, 
            last.Sigma = Sig.ell )
    }
    return(ret)
}



## benchmark("ss"={sapply(1:TT, function(t) sapply(1:(numclust-1), function(kk)
##     pgdraw(Nt.ell[kk,t],abs(XpGamma[kk,t]))))},
##     "mcmapp"={matrix(mcmapply(pgdraw, Nt.ell[-numclust,], abs(XpGamma),
##                               mc.cores = n.cores),
##                      nrow=numclust-1, ncol = TT)},
##     "mapp"={matrix(mapply(pgdraw, Nt.ell[-numclust,], abs(XpGamma)),
##                    nrow=numclust-1, ncol = TT)},
##     replications = 10,
##     columns = c("test", "replications", "elapsed",
##               "relative", "user.self", "sys.self"))
