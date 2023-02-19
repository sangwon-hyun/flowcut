using Distributions 
using LinearAlgebra
using RCall 
using JLD2
using Random 
using StatsBase


include("utils.jl")
fitdata = load("fit_test.jld2")
simdata = load("test.jld2")
y = simdata["y"]
x = simdata["x"]

@rput y x
R"""
png("observed.png")
plot(x=rep(1,100), y=y[[1]], xlim=c(1,30), ylim=c(min(unlist(x)),max(unlist(x))), ylab="y", xlab="t")
for(i in 2:30)
{    
    points(x=rep(i,length(y[[i]])), y=y[[i]])
    # hist(y[[i]])
}
dev.off() 

png("actual.png")
plot(x=rep(1,100), y=x[[1]], xlim=c(1,30), ylim=c(min(unlist(x)),max(unlist(x))), ylab="y", xlab="t")
for(i in 2:30)
{    
    points(x=rep(i,length(x[[i]])), y=x[[i]])
}
dev.off() 
"""

pos = fitdata["pos"]

nsam = length(pos["sigma1"])
nburn = div(nsam, 4)
nthin = div(nsam-nburn,2000)
keep_index = [nburn+1:nthin:nsam;]
nkeep = length(keep_index)

sigma1_save = pos["sigma1"][keep_index]
sigma2_save = pos["sigma2"][keep_index]
mu1_save = pos["mu1"][keep_index,:]
mu2_save = pos["mu2"][keep_index,:]
phi1_save = pos["phi1"][keep_index,:]
w1_save = exp.(phi1_save) ./ (exp.(phi1_save) .+ 1)

@rput sigma1_save sigma2_save 
R"""
par(mfrow=c(2,2))
hist(sigma1_save)
plot(sigma1_save, type="l")
hist(sigma2_save)
plot(sigma2_save, type="l")
par(mfrow=c(1,1))
"""

@rput mu1_save mu2_save 
R"""
mu1_mean = apply(mu1_save, 2, mean) 
mu1_quan = apply(mu1_save, 2, quantile, prob=c(0.025, 0.975))
mu2_mean = apply(mu2_save, 2, mean) 
mu2_quan = apply(mu2_save, 2, quantile, prob=c(0.025, 0.975)) 


df_mu = data.frame(x=1:30, y1=mu1_mean, l1=mu1_quan[1,], h1=mu1_quan[2,],
                           y2=mu2_mean, l2=mu2_quan[1,], h2=mu2_quan[2,])
library(ggplot2)

p = ggplot(df_mu) + geom_line(aes(x=x,y=y1), color="red") 
p = p + geom_ribbon(aes(x=x,ymin=l1,ymax=h1), fill="red", alpha=0.5)
p = p + geom_line(aes(x=x,y=y2), color="blue")
p = p + geom_ribbon(aes(x=x,ymin=l2,ymax=h2), fill="blue", alpha=0.5)
p = p + ylab("mu") + xlab("t") + theme_bw(base_size=25)
ggsave("mu.png", p)
"""

_y_save = pos["_y"][keep_index]
y_impute_mean = mean(_y_save)
@rput y_impute_mean
R"""
png("imputed.png")
plot(x=rep(1,100), y= y_impute_mean[[1]], xlim=c(1,30), ylim=c(min(unlist(y_impute_mean)),max(unlist(y_impute_mean))), col="red", ylab="y", xlab="t")
for(i in 2:30)
{    
    points(x=rep(i,length(y_impute_mean[[i]])), y=y_impute_mean[[i]], col="red")
    # hist(y[[i]])
}
dev.off() 
"""