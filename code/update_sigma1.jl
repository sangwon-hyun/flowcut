function update_sigma1(cur, hyper, dat)
    T = dat["T"]
    Z = cur["Z"]
    _y = cur["_y"]
    
    a1 = hyper["a1"]
    b1 = hyper["b1"]
    mu1 = cur["mu1"]

    N = zeros(T)
    ybar1 = zeros(T) 
    for t in 1:T
        index1 = findall(Z[t] .== 1)
        N[t] = length(index1)
        if N[t] == 0
            ybar1[t] = 0
        else
            ybar1[t] = mean(_y[t][index1])
        end 
    end 

    SigmaTilde1 = Matrix(Diagonal(1 ./ N))

    a1new = a1 + T/2 
    tmp = (mu1 - ybar1)' * SigmaTilde1 * (mu1 - ybar1)
    b1new = 1 / (1/b1 + tmp[1])

    if isnan(b1new)
        print("mu1: ", mu1, "\n")
        print("ybar1: ", ybar1, "\n")
    end

    sigma1new = rand(InverseGamma(a1new, b1new), 1)[1]
    return sigma1new 
end 