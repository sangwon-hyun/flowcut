function update_sigma2(cur, hyper, dat)
    """
    Update σ₂² 
    """
    T = dat["T"]
    nt = dat["nt"]
    Z = cur["Z"]
    _y = cur["_y"]
    
    a2 = hyper["a2"]
    b2 = hyper["b2"]
    mu2 = cur["mu2"]

    N = zeros(T)
    y_mu = 0
    for t in 1:T
        for i in 1:nt[t]
            if Z[t][i] == 2
                N[t] += 1
                y_mu += (_y[t][i] - mu2[t])^2
            end
        end 
    end 

    a2new = a2+sum(N)/2 
    b2new = b2+y_mu/2

    sigma2new = rand(InverseGamma(a2new, b2new), 1)[1]
    return sigma2new 
end 