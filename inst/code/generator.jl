using JLD2 
using StatsBase 
using Distributions 
using Random 
Random.seed!(100)

T = 30 
nt = fill(100, 30)
x = []
y = []
L, R = -3.5, 3.5
nc = 0 
for t in 1:T
    tmp1 = rand(Normal(2, 1), 50)
    tmp2 = rand(Normal(-2,1), 50)
    tmp = vcat(tmp1, tmp2)
    push!(x, tmp)

    nc += length(findall(tmp1 .> R)) + length(findall(tmp2 .< L))
    tmp1[findall(tmp1 .> R)] .= R
    tmp2[findall(tmp2 .< L)] .= L
    tmp = vcat(tmp1, tmp2)
    push!(y, tmp)
end 
print(nc / 3000)

dat = Dict()
dat["y"] = y
dat["x"] = x
dat["L"] = L
dat["R"] = R 

save("../data/test.jld2", dat)