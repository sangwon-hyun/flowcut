## Load the code for running MCMC
## codedir =  "~/Dropbox/Apps/Overleaf/censor-flowmix/BMEcensor/code/"
codedir =  "~/Dropbox/Apps/Overleaf/censor-flowmix/BMEcensor/code/"
source(file.path(codedir, 'BME-functions.r'))

## Load dataset
datadir = "~/Dropbox/research-new/censored-flowmix/code/data"
dta <- readRDS(file.path(datadir, "pre-censor-datobj.RDS"))

## Assign data to objects
Y.list <- dta$ylist
nt <- unlist(lapply(Y.list,function(xx) dim(xx)[1]))
T <- length(Y.list)
d <- 3
W.list = dta$countslist

## Define censor limits
bounds.lower <- rowMins(matrix(unlist(
    lapply(dta$ylist, function(xx) colMins(xx,value = TRUE))),
    nrow = d), value = TRUE) 
bounds.upper <- rowMaxs(matrix(unlist(
    lapply(dta$ylist, function(xx) colMaxs(xx,value = TRUE))),
    nrow = d), value = TRUE)

## Modify censorship slightly
Cbox <- cbind(signif(log(bounds.lower),5),log(bounds.upper))
more.censor <- function(ymat,bd){
    cbind(sapply(ymat[,1], function(xx) ifelse(xx<bd, bd, xx)),
          ymat[,c(2,3)])
}
Y.list <- mclapply(dta$ylist,function(yy) more.censor(log(yy),Cbox[1,1]),
                   mc.cores = min(n.cores, T))
dim(Y.list[[1]])

rowMins(matrix(unlist(
    lapply(Y.list, function(xx) colMins(xx,value = TRUE))),
    nrow = d), value = TRUE) 


X <- dta$X
p <- dim(X)[2] 
head(X) 
X <- t(X)
dim(X) ####  p x T


clean.dta <- list(Y.list = Y.list, W.list = dta$countslist, X=X,
            nt=nt, K=10,T=T ,Cbox=Cbox) 
dim(dta$Y.list[[1]])

Nmc <- 1e3*5
Nburn <- 500
K  <-  10

Gibbs.res0 <- run.Gibbs.fast(Nmc=Nmc, Nburn=Nburn,
                             Y.list = clean.dta$Y.list,
                             W.list = clean.dta$W.list,
                             nt=clean.dta$nt,
                             gg=0.1, 
                             X=clean.dta$X,
                             K = clean.dta$K,
                             d =d,
                             Cbox=clean.dta$Cbox,
                             verbose = TRUE, warm.start = NULL)

