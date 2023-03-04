function update_mu_k(cur, dat, k)

    T = dat["T"]
    Z = cur["Z"]

    sigma = cur["sigma"*string(k)]
    _y = cur["_y"]

    N = zeros(T)
    ybar = zeros(T) 
    for t in 1:T
        index = findall(Z[t] .== k)
        N[t] = length(index)
        if N[t] > 0
            ybar[t] = mean(_y[t][index])
        end
    end 

    SigmaBar_inv = Matrix(Diagonal(N ./ sigma))
    Sigma = get_dist_matrix(1:T, 1:T, cur["l"*string(k)], cur["e"*string(k)])
    Sigma_inv = svd2inv(Sigma) 

    SigmaNew = svd2inv(Sigma_inv + SigmaBar_inv)
    SigmaNew = (SigmaNew + SigmaNew') / 2
    SigmaNew = stablizeMatrix(SigmaNew)
    muNew = SigmaNew * SigmaBar_inv * ybar

    mu_new = rand(MvNormal(muNew, SigmaNew), 1)

    return mu_new 
end 