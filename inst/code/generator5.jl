using JLD2 
using StatsBase 
using Distributions 
using Random 
Random.seed!(100)

T = 30 
nt = fill(100, 30)
x = []
y = []
L, R = -5.3, 5.3
nc = 0 
for t in 1:T
    tmp1 = rand(Normal(4, 1), 30)
    tmp2 = rand(Normal(0, 1), 30)
    tmp3 = rand(Normal(-4,1), 30)
    tmp = vcat(tmp1, tmp2, tmp3)
    push!(x, tmp)

    nc += length(findall(tmp1 .> R)) 
    nc += length(findall(tmp1 .< L)) 
    nc += length(findall(tmp2 .> R)) 
    nc += length(findall(tmp2 .< L))
    nc += length(findall(tmp3 .> R)) 
    nc += length(findall(tmp3 .< L)) 

    tmp1[findall(tmp1 .> R)] .= R
    tmp1[findall(tmp1 .< L)] .= L

    tmp2[findall(tmp2 .> R)] .= R
    tmp2[findall(tmp2 .< L)] .= L

    tmp3[findall(tmp3 .> R)] .= R
    tmp3[findall(tmp3 .< L)] .= L

    tmp = vcat(tmp1, tmp2, tmp3)
    push!(y, tmp)
end 
print(nc / (90*30))

dat = Dict()
dat["y"] = y
dat["x"] = x
dat["L"] = L
dat["R"] = R 

save("../data/test5.jld2", dat)