function update_mu1(cur, dat)

    T = dat["T"]
    Z = cur["Z"]
    sigma1 = cur["sigma1"]
    _y = cur["_y"]

    N = zeros(T)
    ybar1 = zeros(T) 
    for t in 1:T
        index1 = findall(Z[t] .== 1)
        N[t] = length(index1)
        if N[t] > 0
            ybar1[t] = mean(_y[t][index1])
        end
    end 

    SigmaBar1_inv = Matrix(Diagonal(N ./ sigma1))
    Sigma1 = get_dist_matrix(1:T, 1:T, cur["l1"], cur["e1"])
    Sigma1_inv = svd2inv(Sigma1) 

    SigmaNew = svd2inv(Sigma1_inv + SigmaBar1_inv)
    SigmaNew = (SigmaNew + SigmaNew') / 2
    SigmaNew = stablizeMatrix(SigmaNew)
    muNew = SigmaNew * SigmaBar1_inv * ybar1

    mu1_new = rand(MvNormal(muNew, SigmaNew), 1)

    return mu1_new 
end 