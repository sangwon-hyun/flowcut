function impute_data(cur, dat)

    y = dat["y"]
    T = dat["T"]
    nt = dat["nt"]
    L, R = dat["L"], dat["R"] # left and right boundaries  

    Z = cur["Z"]
    _y = deepcopy(cur["_y"])

    for t in 1:T
        for i in 1:nt[t]
            k = Z[t][i]
            d = Normal(cur["mu"*string(k)][t], sqrt(cur["sigma"*string(k)]))
            
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