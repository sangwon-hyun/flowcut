function update_sigma_k(cur, hyper, dat, k)
    """
    Update σₖ² 
    """
    T = dat["T"]
    nt = dat["nt"]
    Z = cur["Z"]
    _y = cur["_y"]
    
    a = hyper["a"*string(k)]
    b = hyper["b"*string(k)]

    mu = cur["mu"*string(k)]

    N = zeros(T)
    y_mu = 0 
    for t in 1:T
        for i in 1:nt[t]
            if Z[t][i] == k
                N[t] += 1
                y_mu += (_y[t][i] - mu[t])^2
            end
        end
    end 

    a_new = a+sum(N)/2 
    b_new = b+y_mu/2

    sigmanew = rand(InverseGamma(a_new, b_new), 1)[1]
    return sigmanew 
end 