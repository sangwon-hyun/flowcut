function phi1_loglikelihood(index, phi1t, _phi1t, Zt, Sigma, nt)
    """
    Calculate log-likelihood of phi1[index]
    """

    Sigma11 = Sigma[index, index]
    Sigma12 = Sigma[index, 1:end .!= index]
    Sigma21 = Sigma[1:end .!= index, index]
    Sigma22 = Sigma[1:end .!= index, 1:end .!= index] 
    Sigma22inv = stablizeMatrix(svd2inv(Sigma22))

    muBar = Sigma12' * Sigma22inv * _phi1t
    SigmaBar = Sigma11 - Sigma12' * Sigma22inv * Sigma21 

    ll = logpdf(Normal(muBar, sqrt(SigmaBar)), phi1t)
    for i in 1:nt
        if Zt[i] == 1
            ll = ll + phi1t - log(exp(phi1t) + 1)
        else
            ll = ll - log(exp(phi1t) + 1)
        end
    end

    return ll 
end 

function update_phi1(cur, dat)
    """
    Update phi1 vector 
    """

    T  = dat["T"]
    nt = dat["nt"]

    Z  = cur["Z"]
    l3 = cur["l3"]
    e3 = cur["e3"]

    Sigma3 = get_dist_matrix(1:T, 1:T, l3, e3)
    Sigma3 = stablizeMatrix(Sigma3)
    
    phi1cur = cur["phi1"]

    for t in 1:T
        phi1pro = phi1cur[t] + rand(Normal(0,0.1), 1)[1]

        llcur = phi1_loglikelihood(t, phi1cur[t], phi1cur[1:end .!= t], Z[t], Sigma3, nt[t])
        llpro = phi1_loglikelihood(t, phi1pro,    phi1cur[1:end .!= t], Z[t], Sigma3, nt[t])

        if log(rand(Uniform(0,1),1)[1]) < (llpro - llcur)
            phi1cur[t] = phi1pro
        end
    end

    return phi1cur 

end 