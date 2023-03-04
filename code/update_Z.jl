function update_Z(cur, hyper, dat)

    T = dat["T"]
    nt = dat["nt"]
    K = hyper["K"]

    _y = cur["_y"]

    Z = []
    for t in 1:T 

        tmp = zeros(Int64, nt[t])
        p = zeros(K)
        for k in 1:(K-1) 
            p[k] = exp(cur["phi"*string(k)][t])
        end
        p[K] = 1

        for i in 1:nt[t]
            unnormalized_p_new = zeros(K)
            for k in 1:K
                unnormalized_p_new[k] = p[k] * pdf(
                    Normal(
                        cur["mu"*string(k)][t], 
                        sqrt(cur["sigma"*string(k)])
                        ), 
                    _y[t][i])
            end
            p_new = unnormalized_p_new ./ sum(unnormalized_p_new)
            tmp[i] = sample(1:K, weights(p_new), 1)[1]
        end
        push!(Z, tmp)
    end 

    return Z 
end 