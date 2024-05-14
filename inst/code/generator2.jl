using JLD2 
using StatsBase 
using Distributions 
using Random 
Random.seed!(100)

T = 30
nt = fill(100, 30)
x = []
y = []
L, R = -2.5, 2.5
mu1 = 0.5 .* sin.(1:T) .+ 1 
mu2 = 0.5 .* cos.(1:T) .- 1
for t in 1:T
    tmp1 = rand(Normal(mu1[t], 1), 50)
    tmp2 = rand(Normal(mu2[t], 1), 50)
    tmp = vcat(tmp1, tmp2)
    push!(x, tmp)

    tmp1[findall(tmp1 .> R)] .= R
    tmp1[findall(tmp1 .< L)] .= L
    tmp2[findall(tmp2 .> R)] .= R
    tmp2[findall(tmp2 .< L)] .= L
    tmp = vcat(tmp1, tmp2)
    push!(y, tmp)
end 

dat = Dict()
dat["y"] = y
dat["x"] = x
dat["L"] = L
dat["R"] = R 

save("../data/test2.jld2", dat)