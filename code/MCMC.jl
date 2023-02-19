using Distributions
using StatsBase
using LinearAlgebra
using TOML 
using JLD2
using ProgressBars
using Random 
Random.seed!(100)

include("utils.jl")
include("impute_data.jl")
include("update_l1_e1.jl")
include("update_l2_e2.jl")
include("update_l3_e3.jl")
include("update_mu1.jl")
include("update_mu2.jl")
include("update_phi1.jl")
include("update_sigma1.jl")
include("update_sigma2.jl")
include("updat_Z.jl")


function MCMC(config_file) 

    config = TOML.parsefile(config_file) 

    nsam = config["nsam"]
    dat = load(config["data_path"] * config["datafile"])
    hyper = config["hyper"]

    # dat is a list of list 
    # yᵢᵗ for i=1,...,nₜ, and t=1,...,T 

    T = length(dat["y"])
    nt = length.(dat["y"])
    dat["T"] = T
    dat["nt"] = nt 

    Z = []
    for t in 1:T
        tmp = zeros(Int64, nt[t])
        for i in 1:nt[t]
            tmp[i] = sample([1,2])
        end 
        push!(Z, tmp) 
    end 

    cur = Dict(
        "mu1" => ones(T),
        "mu2" => -1 .* ones(T),
        "phi1" => zeros(T),
        "Z" => Z, 
        "sigma1" => 1, 
        "sigma2" => 1, 
        "l1" => 1,
        "e1" => 1,
        "l2" => 1,
        "e2" => 1,
        "l3" => 1,
        "e3" => 1,
        "_y" => dat["y"] # imputated data 
    )

    pos = Dict(
        "mu1" => zeros(nsam, T),
        "mu2" => zeros(nsam, T),
        "phi1" => zeros(nsam, T), 
        "sigma1" => zeros(nsam), 
        "sigma2" => zeros(nsam),
        "l1" => zeros(nsam), 
        "l2" => zeros(nsam),
        "l3" => zeros(nsam),
        "e1" => zeros(nsam),
        "e2" => zeros(nsam),
        "e3" => zeros(nsam),
        "_y" => []
    )

    for i in ProgressBar(1:nsam)
        cur["_y"] = impute_data(cur, dat) 
        push!(pos["_y"], cur["_y"])
        cur["Z"] = update_Z(cur, dat)
        pos["mu1"][i,:] = cur["mu1"] = update_mu1(cur, dat) 
        pos["mu2"][i,:] = cur["mu2"] = update_mu2(cur, dat)
        pos["phi1"][i,:] = cur["phi1"] = update_phi1(cur, dat) 
        pos["sigma1"][i] = cur["sigma1"] = update_sigma1(cur, hyper, dat)
        pos["sigma2"][i] = cur["sigma2"] = update_sigma2(cur, hyper, dat)
        cur["l1"], cur["e1"] = update_l1_e1(cur, hyper, dat) 
        pos["l1"][i] = cur["l1"]
        pos["e1"][i] = cur["e1"]
        cur["l2"], cur["e2"] = update_l2_e2(cur, hyper, dat) 
        pos["l2"][i] = cur["l2"]
        pos["e2"][i] = cur["e2"]
        cur["l3"], cur["e3"] = update_l3_e3(cur, hyper, dat) 
        pos["l3"][i] = cur["l3"]
        pos["e3"][i] = cur["e3"]
    end

    result = Dict("pos" => pos,
                  "hyper" => hyper)

    savefile = config["save_path"] 
    save(savefile, result) 
end