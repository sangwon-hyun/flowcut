function impute_data(cur, dat)

    y = dat["y"]
    T = dat["T"]
    nt = dat["nt"]
    L, R = dat["L"], dat["R"] # left and right boundaries  

    Z = cur["Z"]
    mu1 = cur["mu1"]
    mu2 = cur["mu2"]
    sigma1 = cur["sigma1"]
    sigma2 = cur["sigma2"]

    _y = cur["_y"]

    for t in 1:T
        for i in 1:nt[t]
            if Z[t][i] == 1
                d = Normal(mu1[t], sigma1)
            else
                d = Normal(mu2[t], sigma2)
            end 
            
            if L == y[t][i]
                _y[t][i] = rand(truncated(d, upper=L), 1)[1]
            elseif R == y[t][i] 
                _y[t][i] = rand(truncated(d, lower=R), 1)[1]
            else
                _y[t][i] = y[t][i]
            end
        end
    end 

    return _y
end 