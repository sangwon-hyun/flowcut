# Yunzhe
# Version 2
# Date Mar. 2
# Julia 1.7.2 

using Distributions
using StatsBase
using LinearAlgebra
using TOML 
using JLD2
using ProgressMeter
using Random 
Random.seed!(100)

include("utils.jl")
include("impute_data.jl")
include("update_l_e_k.jl")
include("update_s_p_k.jl")
include("update_mu_k.jl")
include("update_phi_k.jl")
include("update_sigma_k.jl")
include("update_Z.jl")


function MCMC(config_file) 

    config = TOML.parsefile(config_file) 

    nsam = config["nsam"]
    dat = load(config["datafile"])
    hyper = config["hyper"]
    K = hyper["K"]

    # dat is a list of list 
    # yᵢᵗ for i=1,...,nₜ, and t=1,...,T 

    T = length(dat["y"])
    nt = length.(dat["y"])
    dat["T"] = T
    dat["nt"] = nt 

    Z = []
    for t in 1:T
        tmp = sample(1:K, nt[t])
        push!(Z, tmp) 
    end 

    cur = Dict() 
    for k in 1:K 
        if k < K
            cur["phi"*string(k)] = zeros(T)
            cur["s"*string(k)] = 2
            cur["p"*string(k)] = 2
        end
        cur["mu"*string(k)] = k .* ones(T)
        cur["sigma"*string(k)] = 1
        cur["l"*string(k)] = 2
        cur["e"*string(k)] = 2
        cur["_y"] = dat["y"]
        cur["Z"] = Z
    end 

    pos = Dict() 
    for k in 1:K
        if k < K 
            pos["phi"*string(k)] = zeros(nsam, T)
            pos["s"*string(k)] = zeros(nsam) 
            pos["p"*string(k)] = zeros(nsam)
        end
        pos["mu"*string(k)] = zeros(nsam, T)
        pos["sigma"*string(k)] = zeros(nsam)
        pos["l"*string(k)] = zeros(nsam) 
        pos["e"*string(k)] = zeros(nsam)
        pos["_y"] = []
        pos["Z"] = [] 
    end

    @showprogress for i in (1:nsam)

        cur["_y"] = impute_data(cur, dat) 
        push!(pos["_y"], cur["_y"])

        cur["Z"] = update_Z(cur, hyper, dat)
        push!(pos["Z"], cur["Z"])

        for k in 1:K 
            if k < K
                pos["phi"*string(k)][i,:] = cur["phi"*string(k)] = update_phi_k(cur, dat, k) 
                cur["s"*string(k)], cur["p"*string(k)] = update_s_p_k(cur, hyper, dat, k) 
                pos["s"*string(k)][i] = cur["s"*string(k)]
                pos["p"*string(k)][i] = cur["p"*string(k)]
            end
            pos["mu"*string(k)][i,:] = cur["mu"*string(k)] = update_mu_k(cur, dat, k) 
            pos["sigma"*string(k)][i] = cur["sigma"*string(k)] = update_sigma_k(cur, hyper, dat, k)
            cur["l"*string(k)], cur["e"*string(k)] = update_l_e_k(cur, hyper, dat, k)
            pos["l"*string(k)][i] = cur["l"*string(k)]
            pos["e"*string(k)][i] = cur["e"*string(k)]
        end
        
    end

    result = Dict("pos" => pos,
                  "hyper" => hyper)

    savefile = config["savepath"] 
    save(savefile, result) 
end