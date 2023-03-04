function phi_loglikelihood(index, _phi, phi_t, _phi_t, Zt, Sigma, nt, k)
    """
    Calculate log-likelihood of phi1[index]
    """

    Sigma11 = Sigma[index, index]
    Sigma12 = Sigma[index, 1:end .!= index]
    Sigma21 = Sigma[1:end .!= index, index]
    Sigma22 = Sigma[1:end .!= index, 1:end .!= index] 
    Sigma22inv = stablizeMatrix(svd2inv(Sigma22))

    muBar = Sigma12' * Sigma22inv * _phi_t
    SigmaBar = Sigma11 - Sigma12' * Sigma22inv * Sigma21 

    ll = logpdf(Normal(muBar, sqrt(SigmaBar)), phi_t)
    for i in 1:nt
        if Zt[i] == k
            ll = ll + phi_t
        else
            ll = ll 
        end
        ll = ll - log(exp(phi_t) + sum(exp.(_phi)) + 1)
    end

    return ll 
end 

function update_phi_k(cur, hyper, dat, k)
    """
    Update phik vector 
    """

    T  = dat["T"]
    nt = dat["nt"]
    K = hyper["K"]

    Z  = cur["Z"]
    s = cur["s"*string(k)]
    p = cur["p"*string(k)]

    Sigma = get_dist_matrix(1:T, 1:T, s, p)
    Sigma = stablizeMatrix(Sigma)
    
    phi_cur = cur["phi"*string(k)]

    for t in 1:T
        phi_rest = []
        for l in 1:(K-1)
            if l != k
                phi_rest = vcat(phi_rest, cur["phi"*string(l)][t])
            end
        end

        phi_pro = phi_cur[t] + rand(Normal(0,0.1), 1)[1]

        llcur = phi_loglikelihood(t, phi_rest, phi_cur[t], phi_cur[1:end .!= t], Z[t], Sigma, nt[t], k)
        llpro = phi_loglikelihood(t, phi_rest, phi_pro,    phi_cur[1:end .!= t], Z[t], Sigma, nt[t], k)

        if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
            phi_cur[t] = phi_pro
        end
    end

    return phi_cur 
end 