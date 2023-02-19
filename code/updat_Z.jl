function update_Z(cur, dat)

    T = dat["T"]
    nt = dat["nt"]
    phi1 = cur["phi1"]
    mu1 = cur["mu1"]
    mu2 = cur["mu2"]
    sigma1 = cur["sigma1"]
    sigma2 = cur["sigma2"]
    _y = cur["_y"]

    Z = []
    for t in 1:T 
        tmp = zeros(Int64, nt[t])
        p1 = exp(phi1[t]) / (exp(phi1[t]) + 1)
        p2 = 1 / (exp(phi1[t]) + 1)
        for i in 1:nt[t]
            p1new = p1 * pdf(Normal(mu1[t], sqrt(sigma1)), _y[t][i])
            p2new = p2 * pdf(Normal(mu2[t], sqrt(sigma2)), _y[t][i])
            p = p1new / (p1new + p2new)
            tmp[i] = sample([1,2], weights([p,1-p]), 1)[1]
        end
        push!(Z, tmp)
    end 

    return Z 
end 