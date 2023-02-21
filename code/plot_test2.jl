using Distributions 
using LinearAlgebra
using RCall 
using JLD2
using Random 
using StatsBase
using ProgressBars


include("utils.jl")
simdata = load("../data/test2.jld2")
y = simdata["y"]
x = simdata["x"]

fitdata = load("../output/fit_test2.jld2")
@rput y x
R"""
# png("observed.png")
plot(x=rep(1,100), y=y[[1]], xlim=c(1,30), ylim=c(min(unlist(x)),max(unlist(x))), ylab="y", xlab="t")
for(i in 2:30)
{    
    points(x=rep(i,length(y[[i]])), y=y[[i]])
    # hist(y[[i]])
}
# dev.off() 

# png("actual.png")
plot(x=rep(1,100), y=x[[1]], xlim=c(1,30), ylim=c(min(unlist(x)),max(unlist(x))), ylab="y", xlab="t")
for(i in 2:30)
{    
    points(x=rep(i,length(x[[i]])), y=x[[i]])
}
# dev.off() 
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

l1_save = pos["l1"][keep_index]
l2_save = pos["l2"][keep_index]
l3_save = pos["l3"][keep_index]
e1_save = pos["e1"][keep_index]
e2_save = pos["e2"][keep_index]
e3_save = pos["e3"][keep_index]
@rput l1_save l2_save l3_save
@rput e1_save e2_save e3_save
R"""
par(mfrow=c(2,3))
plot(l1_save, type="l")
plot(l2_save, type="l")
plot(l3_save, type="l")
plot(e1_save, type="l")
plot(e2_save, type="l")
plot(e3_save, type="l")
par(mfrow=c(1,1))
"""

T = 30
mu1_truth = 0.5 .* sin.(1:T) .+ 1 
mu2_truth = 0.5 .* cos.(1:T) .- 1
@rput mu1_truth mu2_truth
@rput mu1_save mu2_save 
R"""
mu1_mean = apply(mu1_save, 2, mean) 
mu1_quan = apply(mu1_save, 2, quantile, prob=c(0.025, 0.975))
mu2_mean = apply(mu2_save, 2, mean) 
mu2_quan = apply(mu2_save, 2, quantile, prob=c(0.025, 0.975)) 


df_mu = data.frame(x=1:30, y1=mu1_mean, l1=mu1_quan[1,], h1=mu1_quan[2,], t1=mu1_truth, 
                           y2=mu2_mean, l2=mu2_quan[1,], h2=mu2_quan[2,], t2=mu2_truth)
library(ggplot2)

p = ggplot(df_mu) + geom_line(aes(x=x,y=y1), linetype="dashed", color="red") 
p = p + geom_line(aes(x=x,y=t1), color="red")
p = p + geom_ribbon(aes(x=x,ymin=l1,ymax=h1), fill="red", alpha=0.5)
p = p + geom_line(aes(x=x,y=y2), linetype="dashed", color="blue")
p = p + geom_line(aes(x=x,y=t2), color="blue")
p = p + geom_ribbon(aes(x=x,ymin=l2,ymax=h2), fill="blue", alpha=0.5)
p = p + ylab("mu") + xlab("t") + theme_bw(base_size=25)
ggsave("test2_mu.png", p)
p
"""

T_grids = range(1, 30, length=300)
nrep = 100
mu1_pos = zeros(nkeep, length(T_grids), nrep)
for i in ProgressBar(1:nkeep)
    Kvv = get_dist_matrix(1:T, 1:T, l1_save[i], e1_save[i])
    Kwv = get_dist_matrix(T_grids, 1:T, l1_save[i], e1_save[i])
    Kww = get_dist_matrix(T_grids, T_grids, l1_save[i], e1_save[i])

    # K_sigma = Kvv # + sigma1_save[i] * Diagonal(ones(T))
    # K_sigma_inv = svd2inv(K_sigma)
    # mean = Kwv * K_sigma_inv * mu1_save[i,:]
    # Cov = Kww - Kwv * K_sigma_inv * Kwv'
    Kvv_inv = svd2inv(Kvv)
    mean = Kwv * Kvv_inv * mu1_save[i,:]
    Cov = Kww - Kwv * Kvv_inv * Kwv'
    if !isposdef(Cov)
        Cov = (Cov + Cov') / 2
        Cov = Cov + Matrix(Diagonal(10e-7 .* ones(length(T_grids)))) 
    end 
    mu1_pos[i,:,:] = rand(MvNormal(mean, Cov), nrep)
end

mu2_pos = zeros(nkeep, length(T_grids), nrep)
for i in ProgressBar(1:nkeep)
    Kvv = get_dist_matrix(1:T, 1:T, l2_save[i], e2_save[i])
    Kwv = get_dist_matrix(T_grids, 1:T, l2_save[i], e2_save[i])
    Kww = get_dist_matrix(T_grids, T_grids, l2_save[i], e2_save[i])

    Kvv_inv = svd2inv(Kvv)
    mean = Kwv * Kvv_inv * mu2_save[i,:]
    Cov = Kww - Kwv * Kvv_inv * Kwv'
    if !isposdef(Cov)
        Cov = (Cov + Cov') / 2
        Cov = Cov + Matrix(Diagonal(10e-7 .* ones(length(T_grids)))) 
    end 
    mu2_pos[i,:,:] = rand(MvNormal(mean, Cov), nrep)
end

mu1_truth1 = 0.5 .* sin.(T_grids) .+ 1 
mu2_truth1 = 0.5 .* cos.(T_grids) .- 1 
@rput mu1_truth1  mu2_truth1
@rput mu1_pos  mu2_pos  
@rput T_grids
R"""
mu1_pos_mean = apply(mu1_pos, 2, mean)
mu1_pos_quan = apply(mu1_pos, 2, quantile, prob=c(0.025, 0.975))
mu2_pos_mean = apply(mu2_pos, 2, mean)
mu2_pos_quan = apply(mu2_pos, 2, quantile, prob=c(0.025, 0.975))

mu_pos_df = data.frame(x=T_grids, l1=mu1_pos_quan[1,], h1=mu1_pos_quan[2,], m1=mu1_pos_mean, t1=mu1_truth1, 
                                  l2=mu2_pos_quan[1,], h2=mu2_pos_quan[2,], m2=mu2_pos_mean, t2=mu2_truth1)
p_pos_mu = ggplot(mu_pos_df) 
p_pos_mu = p_pos_mu + geom_ribbon(aes(x=x,ymin=l1,ymax=h1), fill="red", alpha=0.5)
p_pos_mu = p_pos_mu + geom_line(aes(x=x,y=m1), color="red", linetype="dashed")
p_pos_mu = p_pos_mu + geom_line(aes(x=x,y=t1), color="red")
p_pos_mu = p_pos_mu + geom_ribbon(aes(x=x,ymin=l2,ymax=h2), fill="blue", alpha=0.5)
p_pos_mu = p_pos_mu + geom_line(aes(x=x,y=m2), color="blue", linetype="dashed")
p_pos_mu = p_pos_mu + geom_line(aes(x=x,y=t2), color="blue")
p_pos_mu = p_pos_mu + xlab("t") + ylab("y") + theme_bw(base_size=25)
ggsave("test2_mu_pos.png", p_pos_mu)
p_pos_mu
"""

@rput w1_save
R"""
w1_mean = apply(w1_save, 2, mean) 
w1_quan = apply(w1_save, 2 ,quantile, prob=c(0.025, 0.975))
df_w1 = data.frame(x=1:30, y1=w1_mean, l1=w1_quan[1,], h1=w1_quan[2,], t1=rep(0.5,30))
pw = ggplot(df_w1) 
pw = pw + geom_line(aes(x=x, y=y1), color="blue", linetype="dashed") 
pw = pw + geom_line(aes(x=x, y=t1), color="blue")
pw = pw + geom_ribbon(aes(x=x,ymin=l1,ymax=h1), fill="blue", alpha=0.5)
ggsave("test2_w1.png", pw)
pw
"""

_y_save = pos["_y"][keep_index]
Z = pos["Z"][keep_index]
# y_impute_mean = mean(_y_save)
# @rput y_impute_mean
# R"""
# # png("imputed.png")
# plot(x=rep(1,100), y= y_impute_mean[[1]], xlim=c(1,30), ylim=c(min(unlist(y_impute_mean)),max(unlist(y_impute_mean))), col="red", ylab="y", xlab="t")
# for(i in 2:30)
# {    
#     points(x=rep(i,length(y_impute_mean[[i]])), y=y_impute_mean[[i]], col="red")
#     # hist(y[[i]])
# }
# # dev.off() 
# """