library("tmvnsim") ## for tmvnsim
library("MASS") ## for mvrnorm
library("scatterplot3d") # for 3d plot 
library("RColorBrewer") ## distinct colors 
library("matrixsampling") ## r-inverse Wishart
library("matrixNormal") ## matrix normal random elements 
library("pgdraw") ## draw PG r.v.'s 
library("parallel") ## mc-mapply
library("svMisc")
library("mclust") ## for warm start 
library("plotly")
library("ggcorrplot")

library("Rfast")
library("Rfast2")

n.cores <- detectCores(); print(paste("n.cores = ",n.cores,sep = ""))

run.Gibbs.fast <- function(Nmc=3000, Nburn=500,
                           Y.list, W.list, X, 
                           nt, K, d,
                           Cbox=NULL,  
                           user.prior = NULL, gg=1,
                           verbose = FALSE,
                           plot.cytogram = TRUE, 
                           warm.start = NULL, tol = 1/1e8){
    T <- length(Y.list)
    p <- dim(X)[1]
    NN <- sum(nt)
    tt.impute <- min(20,floor(Nburn/5))
    Y.list.raw <- Y.list ## store raw Y.list 
##### pre computed     
    X.list <- as.list(as.data.frame(X))
    Xp <- rbind(1,X) 
    Xp.list <- as.list(as.data.frame(Xp))
    W.sq.list <- lapply(W.list, function(xx) sqrt(xx))
    mt <- lapply(W.list, sum)%>%unlist()
    MM <- sum(mt) 
        
    XtXtT <- lapply(X.list,function(x) x%*%t(x))
    XtXtTp <- lapply(Xp.list, function(x) x%*%t(x))
    ggXtXtTp <- lapply(Xp.list, function(x) x%*%t(x)/gg)
    XTX <- X%*%t(X)
    XTXp <- Xp%*%t(Xp)
    inv.XTX <- Rfast::spdinv(XTX)
    inv.XTXp <- Rfast::spdinv(XTXp)
##    Y.grand.mean <- Rfast::colmeans(do.call(rbind,Y.list))


    if(!is.null(Cbox)){
        Censor.list <- mclapply(Y.list, function(x) censorIndicator(x,Cbox),
                                mc.cores = min(n.cores, T))
        censor.01.list <- mclapply(Censor.list, function(x)
            apply(x,1, function(xx)
                sum(abs(xx))>0), mc.cores = min(n.cores, T))
        censor.1.list <- mclapply(Censor.list, function(x)
            apply(x,1, function(xx)
                sum(abs(xx))==1), mc.cores = min(n.cores, T)) 
        censor.23.list <- mclapply(Censor.list, function(x)
            apply(x,1, function(xx)
                sum(abs(xx))>1), mc.cores = min(n.cores, T)) 
        censored.Y.list <- mcmapply(function(yy,c01) {yy[c01==TRUE,]},
                                    yy = Y.list, c01 = censor.01.list,
                                    mc.cores = min(n.cores, T))
        censored.C.list <- mcmapply(function(c01,cc){cc[c01==TRUE,]},
                                    c01 = censor.01.list, cc=Censor.list,
                                    mc.cores = min(n.cores, T))
        censored.W.list <- mcmapply(function(c01,ww){ww[c01==TRUE]},
                                    c01 = censor.01.list, w=W.list,
                                    mc.cores = min(n.cores, T))  
        nt.censor <- sapply(1:T,function(tt) sum(censor.01.list[[tt]])) 
        samp.region.list <- lapply(1:T, function(tt){
            c01 <- censor.01.list[[tt]]
            cc <- Censor.list[[tt]]
            if(verbose){
                print(paste("Time=",tt,", censored obs.= ",sum(c01==TRUE),
                            ", censored ratio =", round(sum(c01==TRUE)/nt[tt],2), sep=" ")) 
            }
            if(sum(c01==TRUE)>1){
                apply(cc[c01==TRUE,],1,function(xx) sample.region(xx,Cbox))
            }else if(sum(c01==TRUE)==1){
                matrix(sample.region(cc[c01==TRUE,],Cbox),ncol=1)
            }else{
                NULL
            }}) 
    }

##### prior specifications
    if(is.null(user.prior)){
        nu0=d
        nu1=p+1 
        S0=diag(d)
        S1=diag(p+1)
        ## inv.Omega <- solve(S1)
    }

##### initialize
    if(!is.null(warm.start)){
        print("continue with a previous draw") 
        Z.list <- warm.start$Z.list
        Y.list <- warm.start$Y.list
        ## Omega <- warm.start$Omega 
        beta.ell <-  warm.start$beta 
        gamma.ell <- warm.start$gamma
        Sig.ell <- warm.start$Sigma  
        
        Nburn <- 0 ## no need of burn-in 
        tt.impute <- min(25,floor(Nburn/5)) 
    }else{
        Z.list <- lapply(1:T,function(tt) sample(1:K,nt[tt],replace = TRUE)) 
        ## Omega <- rinvwishart(1,nu1+p+1,S1)[,,1]
        beta.ell <-  array(rnorm(K*(p+1)*d), c(d,p+1,K))
        gamma.ell <- matrix(rnorm((p+1)*(K-1)),nrow=p+1,ncol=K-1) 
        Sig.ell <- rinvwishart(K,nu0+d,S0)
    }

##    Omega.ell <- rinvwishart(K-1,nu1+p+1,S1)
##    beta.ell <-  array(rnorm(K*p*d), c(d,p,K))
##    beta0.ell <- matrix(0,nrow=d,ncol=K)

    SX.ell <- array(0, c(p,p,K))
    Sy.ell <- matrix(0,nrow = d, ncol = T)
    Sxy.ell <- array(0, c(p,d,K))
    
##### store
    ## burn.beta0 <- array(0,c(d,K,Nburn))
    ## burn.beta <- array(0,c(d,p,K,Nburn))
    burn.beta <- array(0,c(d,p+1,K,Nburn))
    burn.Sigma <- array(0,c(d,d,K,Nburn)) 
    burn.gamma <- array(0,c(p+1,K-1,Nburn))
    ## burn.Omega <- array(0,c(p+1,p+1,Nburn))
    burn.avgloglik <- rep(NA, Nburn)
    
    ## pos.beta0 <- array(0,c(d,K,Nmc))
    ##pos.beta <- array(0,c(d,p,K,Nmc))
    pos.beta <- array(0,c(d,p+1,K,Nmc))
    pos.Sigma <- array(0,c(d,d,K,Nmc)) 
    pos.gamma <- array(0,c(p+1,K-1,Nmc))
    ##pos.Omega <- array(0,c(p+1,p+1,Nmc))
    pos.avgloglik <- rep(NA, Nmc) 
    ## if(is.null(Cbox)==FALSE){
    ##     pos.imputed.Y <- array(0,c(d,sum(nt.censor),Nmc))
    ## } 

    if(verbose){
        print("The Gibbs sampler starts now.")
        list.iter <- c(1, Nburn+1, 1:floor((Nmc+ Nburn-1)/100)*100, Nmc + Nburn)
        ptm <- NULL 
    }
    
##### posterior sampling starts here
    for(jj in 1:(Nburn+Nmc)){
        print(jj) 
        if(verbose){
            if(jj ==1 & Nburn > 0){
                print("Burn-in period starts.")
                print(Sys.time())                
            }
            if(jj == tt.impute + 1){
                ptm <- proc.time()
            }
            if(jj ==Nburn+1 & Nburn > 0){
                ## burn.in.last.draw <- list(Z.list=Z.list,Y.list=Y.list)
                ## save(burn.in.last.draw, file = "Burn-in-last-draw.Rdata")
                ## try(dev.off(),silent = TRUE)
                ## par(mfrow = c(d,K))
                ## par(mar = c(K,d,K,d)/1.5)
                ## try(for(dd in 1:d){
                ##         for(kk in 1:K){
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
            sapply(1:K, function(kk) sum(ww[zz==kk]))},
            ww = W.list, zz = Z.list, SIMPLIFY = TRUE,
            mc.cores = min(n.cores, T))
        m.ell <-  Rfast::rowsums(mt.ell)

        ## nt.ell  <- do.call(cbind,mclapply(Z.list, function(zz){
        ##     sapply(1:K, function(kk) as.numeric(sum(zz == kk)))},
        ##     mc.cores = min(n.cores, T))) 
        ## n.ell <- Rfast::rowsums(nt.ell)
        
        for(ell in 1:K){ 
            ## mm0t <- do.call(rbind,mapply(function(xx,yy,zz){
            ##     Rfast::colsums(Rfast::eachrow(as.matrix(yy[zz==ell,]), beta.ell[,,ell] %*% xx,"-"))},
            ##     xx=X.list, yy = Y.list, zz = Z.list, SIMPLIFY = FALSE))
            ## mm0 <- Rfast::colsums(mm0t)/n.ell[ell] 
            ## beta0.ell[,ell] <- mvrnorm(1, mm0, Sig.ell[,,ell]/n.ell[ell])

            ## SX.ell <- XTXp/gg + Reduce('+',Map(`*`, XtXtTp, nt.ell[ell,]))  
            ## inv.SX.ell <- Rfast::spdinv(SX.ell) ## (p+1) x (p+1)
            ## Sy.ell <- do.call(rbind, mcmapply(function(yy,zz){
            ##     Rfast::colsums(as.matrix(yy[zz==ell,]))},
            ##     yy=Y.list,zz=Z.list, SIMPLIFY = FALSE,
            ##     mc.cores = min(n.cores, T)))
            

            SX.ell <- XTXp/gg + Reduce('+',Map(`*`, XtXtTp, mt.ell[ell,]))  
            inv.SX.ell <- Rfast::spdinv(SX.ell) ## (p+1) x (p+1)
            Sy.ell <- mcmapply(function(ww, yy, zz){
                t(as.matrix(ww[zz==ell]))%*%yy[zz==ell,]},
                ww = W.list, yy = Y.list, zz = Z.list, SIMPLIFY = TRUE,
                mc.cores = min(n.cores, T)) 
            Sxy.ell <- Xp %*% t(Sy.ell)  ## (p+1) x d
            beta.ell[,,ell] <- rmatnorm(M = t(Sxy.ell) %*% inv.SX.ell, 
                                        U = Sig.ell[,,ell], 
                                        V = inv.SX.ell, 
                                        tol = .Machine$double.eps^0.95)
            
            sse <- mcmapply(function(ww,xx,yy,zz) {
                wcrossprod.fast(Rfast::eachrow(as.matrix(yy[zz==ell,]),
                                               beta.ell[,,ell]%*%xx, '-'),
                                ww[zz==ell], weighting = TRUE)},
                xx = Xp.list, yy = Y.list, zz=Z.list, ww = W.sq.list, 
                SIMPLIFY = FALSE,
                mc.cores = min(n.cores, T)) 
            ## sse <- mcmapply(function(ww,xx,yy,zz) {
            ##     wcrossprod.fast(Rfast::eachrow(as.matrix(yy[zz==ell,]),
            ##                                    beta.ell[,,ell]%*%xx, '-'),
            ##                     ww[zz==ell], weighting = FALSE)},
            ##     xx = Xp.list, yy = Y.list, zz=Z.list, ww = W.sq.list, 
            ##     SIMPLIFY = FALSE,
            ##     mc.cores = min(n.cores, T)) 
            Sn.ell <- Reduce('+',sse) 
            Sig.ell[,,ell] <- rinvwishart(1,nu0+d + m.ell[ell], S0 + Sn.ell)[,,1]
        }
        
################################################ 
        ## expert assignment ###
################################################
        ## mt.ell <- mcmapply(function(ww,zz){
        ##     sapply(1:K, function(kk) sum(ww[zz==kk]))},
        ##     ww = W.list, zz = Z.list, SIMPLIFY = TRUE,
        ##     mc.cores = min(n.cores, T))
        ## m.ell <-  Rfast::rowsums(mt.ell)

        ## XpGamma.abs <- abs(t(gamma.ell) %*% Xp) ## K-1 x T 
        ## mt.cumsum <- Rfast::colCumSums(mt.ell) 
        ## Mt.ell <- rbind(nt,-sweep(mt.cumsum[-K,], 2, mt.cumsum[K,])) 
        ## omega.tell <- matrix(mcmapply(pgdraw, round(Mt.ell[-K,]), XpGamma.abs,
        ##                               mc.cores = min(n.cores, T)),
        ##                      nrow=K-1, ncol = T)
        ## kappa <- mt.ell[-K,] - Mt.ell[-K,]/2


        XpGamma.abs <- abs(t(gamma.ell) %*% Xp) ## K-1 x T  
        ## nt.cumsum <- Rfast::colCumSums(as.matrix(nt.ell))  
        ## Nt.ell <- rbind(nt,-sweep(nt.cumsum[-K,], 2, nt.cumsum[K,])) 
        ## omega.tell <- matrix(mcmapply(pgdraw, Nt.ell[-K,], XpGamma.abs,
        ##                               mc.cores = min(n.cores, T)),
        ##                      nrow=K-1, ncol = T)
        ## kappa <- nt.ell[-K,] - Nt.ell[-K,]/2

        mt.cumsum <- Rfast::colCumSums(as.matrix(mt.ell))  
        Mt.ell <- rbind(mt,-sweep(mt.cumsum[-K,], 2, mt.cumsum[K,])) 
        omega.tell <- matrix(mcmapply(pgdraw, round(Mt.ell[-K,]), XpGamma.abs,
                                      mc.cores = min(n.cores, T)),
                             nrow=K-1, ncol = T)
        kappa <- mt.ell[-K,] - Mt.ell[-K,]/2

        ##        inv.Omega <- Rfast::spdinv(Omega)
        for(ell in 1:(K-1)){
            V.omell <- Rfast::spdinv(Reduce('+', Map('*',XtXtTp, omega.tell[ell,]))) ##+diag(p+1)/T ) Mod2 
            m.omell <- V.omell%*% Reduce('+', Map('*',Xp.list,kappa[ell,]))
            gamma.ell[,ell] <- Rfast::rmvnorm(1,m.omell, V.omell)
        }        


        XpGamma <- t(gamma.ell) %*% Xp ## K-1 x T 
        pi.sb <- 1/(1+exp(-XpGamma))
        pi.mn <- apply(pi.sb,2,SB2MN)
        logpi.list <- mclapply(1:T,function(t) log(pi.mn[,t]),
                               mc.cores = min(n.cores,T)) 
        
        chol.Sig.ell <- apply(Sig.ell,3, chol)
        chol.Sig.list <- lapply(1:K,function(kk) matrix(chol.Sig.ell[,kk],nrow = d))
        mu.list <- lapply(1:T, function(tt) sapply(1:K, function(kk)
            beta.ell[,,kk] %*% Xp[,tt]))
#            beta0.ell[,kk] + beta.ell[,,kk] %*% X[,tt])) 
        
        logPiZ <- mcmapply(function(ww, xx, yy, mm, pp){
            sapply(1:K, function(kk)
                mvnfast::dmvn(yy, mm[,kk], chol.Sig.list[[kk]],
                              log=TRUE, isChol = TRUE) + pp[kk])}, ## mod here 
            ww = W.list, xx =X.list, yy = Y.list, mm = mu.list, pp=logpi.list,
            mc.cores = min(n.cores, T), SIMPLIFY = FALSE)

        Z.list <- mclapply(logPiZ, function(pp) apply(pp,1, function(lpi)
            sample(1:K,1,prob=softmax(lpi))),
            mc.cores = min(n.cores, T))

        ## mt.ell <- mcmapply(function(ww,zz){
        ##     sapply(1:K, function(kk) sum(ww[zz==kk]))},
        ##     ww = W.list, zz = Z.list, SIMPLIFY = TRUE,
        ##     mc.cores = min(n.cores, T))
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
                                        mc.cores = min(n.cores, T))  
            imputed.Y.list <- mclapply(1:T, function(tt) {
                impute.censored(ww = censored.W.list[[tt]], 
                                yy = censored.Y.list[[tt]],
                                zz = censored.Z.list[[tt]],
                                cc.info.mat =  censored.C.list[[tt]],
                                bounds.mat = samp.region.list[[tt]],
                                mu.mat = mu.list[[tt]], Sigma.ell = Sig.ell)},
                mc.cores = min(n.cores, T), mc.preschedule = FALSE)             
            for(tt in 1:T){
                Y.list[[tt]][censor.01.list[[tt]]==TRUE,] <- imputed.Y.list[[tt]]
            }            
        }

        loglik <- loglik_eval(mu.list, chol.Sig.list, 
                              W.list, X.list, Y.list, Z.list,
                              as.list(nt),
                              simple = TRUE) 
        ## print(sort(round(n.ell/NN,3)))
        print(sort(round(m.ell/MM,3)))
        print(paste("avg loglikelihood: ", round(loglik,2), sep=" "))

################################################ 
            ## graphics: clustering results  
################################################ 

        if(plot.cytogram){
            plt.tt <- jj %% (T-6) 
            cols <- brewer.pal(K,"Paired") ## K need to be <=12
            pchs <- c(16,17)
            dim.names <- colnames(Y.list[[1]])
            dims <- c(1,sample(2:3,1))  
            par(mfrow = c(2,2)) 
            for(plt.time in (plt.tt + c(0,6))){
                plt.dta <- Y.list.raw[[plt.time]]
                title.imputed <- paste("Time=",plt.time,
                                       ", K=",K,
                                       ", iter=", jj, sep=" ")
                plot(plt.dta[,dims], 
                     pch = 16,
                     cex = .5,
                     xlab = dim.names[dims[1]], 
                     ylab = dim.names[dims[2]], 
                     col="steelblue",
                     main=title.imputed)
            }
            for(plt.time in (plt.tt + c(0,6))){
                plt.dta.imputed <- Y.list[[plt.time]]
                plt.zz <- Z.list[[plt.time]]
                plt.cc1 <- censor.1.list[[plt.time]]
                plt.cc23 <- censor.23.list[[plt.time]]
                plt.cc <- plt.cc1+plt.cc23 
                title.imputed <- paste("Time=",plt.time,
                                       ", K=",K,
                                       ", iter=", jj, sep=" ")
                ## scatterplot3d(plt.dta.imputed,
                ##               pch = pchs[plt.cc23+1],
                ##               color=cols[plt.zz],
                ##               main=title.imputed, 
                ##               grid=FALSE, box=FALSE,
                ##               angle = plt.angle)
                plot(plt.dta.imputed[,dims], 
                     pch = pchs[plt.cc+1],
                     cex = .5,
                     xlab = dim.names[dims[1]], 
                     ylab = dim.names[dims[2]], 
                     col=cols[plt.zz],
                     main=title.imputed)
                points(t(mu.list[[plt.time]][dims,]),
                       pch=8, cex=1, 
                       col=cols)
                ## points(t(mu.list[[plt.time]][dims,]),
                ##        pch=8, cex=0.25, 
                ##        col="black")
            }
        }


        

       ##  if(jj %in% list.iter & verbose){
       ##      print(sort(round(n.ell/NN,3)))
       ##      print(min(n.ell))
       ##      print(sort(round(m.ell/MM,3)))
       ##      print(min(m.ell))
       ##      print(paste("avg loglikelihood: ", round(loglik,2), sep=" "))
       ## }

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
            ##     pos.imputed.Y[,,jj-Nburn] <- matrix(unlist(imputed.Y.list),nrow = d) 
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
            pos.imputed.Y = do.call(rbind,imputed.Y.list),
            burn.avgloglik = burn.avgloglik,
            pos.avgloglik = pos.avgloglik,
            raw.Y.list = Y.list.raw, 
            last.Y.list = Y.list,
            last.Z.list = Z.list,
            last.gamma = gamma.ell,
            ## last.Omega = Omega,
            last.beta = beta.ell, 
            last.Sigma = Sig.ell )
    }
    return(ret)
}



## benchmark("ss"={sapply(1:T, function(t) sapply(1:(K-1), function(kk)
##     pgdraw(Nt.ell[kk,t],abs(XpGamma[kk,t]))))},
##     "mcmapp"={matrix(mcmapply(pgdraw, Nt.ell[-K,], abs(XpGamma),
##                               mc.cores = min(n.cores, T)),
##                      nrow=K-1, ncol = T)},
##     "mapp"={matrix(mapply(pgdraw, Nt.ell[-K,], abs(XpGamma)),
##                    nrow=K-1, ncol = T)},
##     replications = 10,
##     columns = c("test", "replications", "elapsed",
##               "relative", "user.self", "sys.self"))

