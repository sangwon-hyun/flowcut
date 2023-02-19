function update_sigma2(cur, hyper, dat)
    T = dat["T"]
    Z = cur["Z"]
    _y = cur["_y"]
    
    a2 = hyper["a2"]
    b2 = hyper["b2"]
    mu2 = cur["mu2"]

    N = zeros(T)
    ybar2 = zeros(T) 
    for t in 1:T
        index2 = findall(Z[t] .== 2)
        N[t] = length(index2)
        if N[t] == 0
            ybar2[t] = 0
        else
            ybar2[t] = mean(_y[t][index2])
        end
    end 

    SigmaTilde2 = Matrix(Diagonal(1 ./ N))

    a2new = a2 + T/2 

    tmp = (mu2 - ybar2)' * SigmaTilde2 * (mu2 - ybar2)
    b2new = 1 / (1 / b2 + tmp[1])

    sigma2new = rand(InverseGamma(a2new, b2new), 1)[1]
    return sigma2new 
end 