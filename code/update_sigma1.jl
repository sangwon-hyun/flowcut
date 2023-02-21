function update_sigma1(cur, hyper, dat)
    """
    Update σ₁² 
    """
    T = dat["T"]
    nt = dat["nt"]
    Z = cur["Z"]
    _y = cur["_y"]
    
    a1 = hyper["a1"]
    b1 = hyper["b1"]
    mu1 = cur["mu1"]

    N = zeros(T)
    y_mu = 0 
    for t in 1:T
        for i in 1:nt[t]
            if Z[t][i] == 1
                N[t] += 1
                y_mu += (_y[t][i] - mu1[t])^2
            end
        end
    end 

    a1new = a1+sum(N)/2 
    b1new = b1+y_mu/2

    sigma1new = rand(InverseGamma(a1new, b1new), 1)[1]
    return sigma1new 
end 