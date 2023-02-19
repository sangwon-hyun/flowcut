function update_mu2(cur, dat)

    T = dat["T"]
    Z = cur["Z"]
    sigma2 = cur["sigma2"]
    # mu2 = cur["mu2"]
    _y = cur["_y"]

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
    SigmaBar2_inv = (1 / sigma2) .* Matrix(Diagonal(N))

    Sigma2 = get_dist_matrix(1:T, cur["l1"], cur["e1"])

    SigmaNew = svd2inv(svd2inv(Sigma2) + SigmaBar2_inv)
    muNew = SigmaNew * SigmaBar2_inv * ybar2

    mu2_new = rand(MvNormal(muNew, SigmaNew), 1)

    return mu2_new 
end 